from __future__ import annotations

from operator import truediv
from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils import stats
from fast_array_utils.numba import njit
from fast_array_utils.stats import mean_var

from .. import logging as logg
from .._compat import CSBase, CSCBase, CSRBase, DaskArray, warn
from .._utils import axis_mul_or_truediv, dematrix, view_to_actual
from ..get import _get_obs_rep, _set_obs_rep

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


def _compute_nnz_median(counts: np.ndarray | DaskArray) -> np.floating:
    """Given a 1D array of counts, compute the median of the non-zero counts."""
    if isinstance(counts, DaskArray):
        counts = counts.compute()
    counts_greater_than_zero = counts[counts > 0]
    median = np.median(counts_greater_than_zero)
    return median


@njit
def _normalize_csr(
    mat: CSRBase,
    *,
    rows,
    columns,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
):
    """For sparse CSR matrix, compute the normalization factors."""
    counts_per_cell = np.zeros(rows, dtype=mat.data.dtype)
    for i in numba.prange(rows):
        count = 0.0
        for j in range(mat.indptr[i], mat.indptr[i + 1]):
            count += mat.data[j]
        counts_per_cell[i] = count
    if exclude_highly_expressed:
        counts_per_cols_t = np.zeros((n_threads, columns), dtype=np.int32)
        counts_per_cols = np.zeros(columns, dtype=np.int32)

        for i in numba.prange(n_threads):
            for r in range(i, rows, n_threads):
                for j in range(mat.indptr[r], mat.indptr[r + 1]):
                    if mat.data[j] > max_fraction * counts_per_cell[r]:
                        minor_index = mat.indices[j]
                        counts_per_cols_t[i, minor_index] += 1
        for c in numba.prange(columns):
            counts_per_cols[c] = counts_per_cols_t[:, c].sum()

        for i in numba.prange(rows):
            count = 0.0
            for j in range(mat.indptr[i], mat.indptr[i + 1]):
                if counts_per_cols[mat.indices[j]] == 0:
                    count += mat.data[j]
            counts_per_cell[i] = count

    return counts_per_cell, counts_per_cols


def _normalize_total_helper(
    x: np.ndarray | CSBase | DaskArray,
    *,
    exclude_highly_expressed: bool,
    max_fraction: float,
    target_sum: float | None,
) -> tuple[np.ndarray | CSBase | DaskArray, np.ndarray, np.ndarray | None]:
    """Calculate the normalized data, counts per cell, and gene subset.

    Parameters
    ----------
    See `normalize_total` for details.

    Returns
    -------
    X
        The normalized data matrix.
    counts_per_cell
        The normalization factors used for each cell (counts / target_sum).
    gene_subset
        If `exclude_highly_expressed=True`, a boolean mask indicating which genes
        were not considered highly expressed. Otherwise, `None`.
    """
    gene_subset = None
    counts_per_cell = None
    if isinstance(x, CSRBase):
        n_threads = numba.get_num_threads()
        counts_per_cell, counts_per_cols = _normalize_csr(
            x,
            rows=x.shape[0],
            columns=x.shape[1],
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            n_threads=n_threads,
        )
        if target_sum is None:
            target_sum = np.median(counts_per_cell)
        if exclude_highly_expressed:
            gene_subset = ~np.where(counts_per_cols)[0]
    else:
        counts_per_cell = stats.sum(x, axis=1)
        if exclude_highly_expressed:
            # at least one cell as more than max_fraction of counts per cell
            hi_exp = dematrix(x > counts_per_cell[:, None] * max_fraction)
            gene_subset = stats.sum(hi_exp, axis=0) == 0

            counts_per_cell = stats.sum(x[:, gene_subset], axis=1)
        if target_sum is None:
            target_sum = _compute_nnz_median(counts_per_cell)

    counts_per_cell = counts_per_cell / target_sum
    out = x if isinstance(x, np.ndarray | CSBase) else None
    x = axis_mul_or_truediv(
        x, counts_per_cell, op=truediv, out=out, allow_divide_by_zero=False, axis=0
    )
    return x, counts_per_cell, gene_subset


def normalize_total(  # noqa: PLR0912
    adata: AnnData,
    *,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: str | None = None,
    layer: str | None = None,
    obsm: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | dict[str, np.ndarray] | None:
    """Normalize counts per cell.

    Normalize each cell by total counts over all genes,
    so that every cell has the same total count after normalization.
    If choosing `target_sum=1e6`, this is CPM normalization.

    If `exclude_highly_expressed=True`, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes :cite:p:`Weinreb2017`.

    Similar functions are used, for example, by Seurat :cite:p:`Satija2015`, Cell Ranger
    :cite:p:`Zheng2017` or SPRING :cite:p:`Weinreb2017`.

    .. array-support:: pp.normalize_total

    .. note::
        When used with a :class:`~dask.array.Array` in `adata.X`, this function will have to
        call functions that trigger `.compute()` on the :class:`~dask.array.Array` if `exclude_highly_expressed`
        is `True`, `layer_norm` is not `None`, or if `key_added` is not `None`.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    target_sum
        If `None`, after normalization, each observation (cell) has a total
        count equal to the median of total counts for observations (cells)
        before normalization.
    exclude_highly_expressed
        Exclude (very) highly expressed genes for the computation of the
        normalization factor (size factor) for each cell. A gene is considered
        highly expressed, if it has more than `max_fraction` of the total counts
        in at least one cell. The not-excluded genes will sum up to
        `target_sum`.  Providing this argument when `adata.X` is a :class:`~dask.array.Array`
        will incur blocking `.compute()` calls on the array.
    max_fraction
        If `exclude_highly_expressed=True`, consider cells as highly expressed
        that have more counts than `max_fraction` of the original total counts
        in at least one cell.
    key_added
        Name of the field in `adata.obs` where the normalization factor is
        stored.
    layer
        Layer to normalize instead of `X`.
    obsm
        Array to normalize instead of `X`.
    inplace
        Whether to update `adata` or return dictionary with normalized copies of
        `adata.X` and `adata.layers`.
    copy
        Whether to modify copied input object. Not compatible with inplace=False.

    Returns
    -------
    Returns dictionary with normalized copies of `adata.X` and `adata.layers`
    or updates `adata` with normalized version of the original
    `adata.X` and `adata.layers`, depending on `inplace`.

    Example
    -------
    >>> import sys
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> sc.settings.verbosity = "info"
    >>> sc.settings.logfile = sys.stdout  # for doctests
    >>> np.set_printoptions(precision=2) and None
    >>> adata = AnnData(
    ...     np.array(
    ...         [
    ...             [3, 3, 3, 6, 6],
    ...             [1, 1, 1, 2, 2],
    ...             [1, 22, 1, 2, 2],
    ...         ],
    ...         dtype="float32",
    ...     )
    ... )
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)["X"]
    normalizing counts per cell
        finished (0:00:00)
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(
    ...     adata,
    ...     target_sum=1,
    ...     exclude_highly_expressed=True,
    ...     max_fraction=0.2,
    ...     inplace=False,
    ... )["X"]
    normalizing counts per cell
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
        finished (0:00:00)
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)

    """
    if copy:
        if not inplace:
            msg = "`copy=True` cannot be used with `inplace=False`."
            raise ValueError(msg)
        adata = adata.copy()

    if max_fraction < 0 or max_fraction > 1:
        msg = "Choose max_fraction between 0 and 1."
        raise ValueError(msg)

    view_to_actual(adata)

    x = _get_obs_rep(adata, layer=layer, obsm=obsm)
    if isinstance(x, CSCBase):
        x = x.tocsr()
    if not inplace:
        x = x.copy()
    if issubclass(x.dtype.type, int | np.integer):
        x = x.astype(np.float32)  # TODO: Check if float64 should be used

    start = logg.info("normalizing counts per cell")

    x, counts_per_cell, gene_subset = _normalize_total_helper(
        x,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        target_sum=target_sum,
    )

    if exclude_highly_expressed:
        logg.info(
            "The following highly-expressed genes are not considered during normalization factor computation:\n"
            f"{adata.var_names[~gene_subset].tolist()}"
        )

    cell_subset = counts_per_cell > 0
    if not isinstance(cell_subset, DaskArray) and not np.all(cell_subset):
        warn("Some cells have zero counts", UserWarning)

    dat = dict(
        X=x,
        norm_factor=counts_per_cell,
    )
    if inplace:
        if key_added is not None:
            adata.obs[key_added] = dat["norm_factor"]
        _set_obs_rep(adata, dat["X"], layer=layer, obsm=obsm)

    logg.info(
        "    finished ({time_passed})",
        time=start,
    )
    if key_added is not None:
        logg.debug(
            f"and added {key_added!r}, counts per cell before normalization (adata.obs)"
        )

    if copy:
        return adata
    elif not inplace:
        return dat
    return None


# -----------------------------------------------------------------------------
# Shifted Centered Log-Ratio (PFlogPF) Core Implementation
#
# The mathematical logic and sparse matrix optimization ("offset trick")
# implemented below are adapted from the Pachter Lab reference repository:
# https://github.com/pachterlab/bhgp_2022
#
# Reference Publication:
# Booeshaghi et al. (2022, 2026) "Depth normalization for single-cell genomics
# count data" bioRxiv. doi: 10.1101/2022.05.06.490859
#
# Copyright (c) 2022, Pachter Lab. All rights reserved.
# Licensed under the BSD 2-Clause License.
# -----------------------------------------------------------------------------


def _estimate_overdispersion(x: np.ndarray | CSBase | DaskArray) -> float:
    r"""Estimate the negative-binomial overdispersion :math:`α` from raw counts.

    Fits :math:`\mathrm{Var}_g = μ_g + α \cdot μ_g^2` across genes, where
    :math:`μ_g` and :math:`\mathrm{Var}_g` are the per-gene mean and (population)
    variance over cells. The model is linear in :math:`α`, so the ordinary
    least-squares solution is closed form

    .. math::
        α = \frac{\sum_g (\mathrm{Var}_g - μ_g) \, μ_g^2}{\sum_g μ_g^4},

    which is exactly the minimizer a non-linear `curve_fit` would converge to,
    but without the dependency. :func:`~fast_array_utils.stats.mean_var` is
    dispatched for dense, sparse and dask input alike, so the only dask-specific
    step is computing the two final scalar sums.
    """
    mu, var = mean_var(x, axis=0, correction=0)
    mu2 = mu**2
    numerator = np.sum((var - mu) * mu2)
    denominator = np.sum(mu2 * mu2)
    if isinstance(x, DaskArray):
        import dask

        numerator, denominator = dask.compute(numerator, denominator)
    if denominator == 0.0:
        msg = (
            "Cannot estimate overdispersion: every gene has zero mean. "
            "Pass `alpha` or `target_sum` explicitly."
        )
        raise ValueError(msg)
    return float(numerator / denominator)


def _clr_log_center(u: np.ndarray | CSBase) -> np.ndarray:
    """Apply log1p then per-cell mean-centering to PF-scaled counts; return dense.

    Operates on a full-width block of cells: every gene must be present so the
    per-cell mean spans the whole row. Used directly for in-memory input, and
    per row-chunk under :meth:`~dask.array.Array.map_blocks` (after the genes
    have been rechunked into a single chunk).
    """
    n_genes = u.shape[1]
    if isinstance(u, CSBase):
        # Sparse "offset trick": store log1p only at the nonzero entries. The
        # structural zeros already hold log1p(0) = 0, so the per-cell mean is just
        # the stored row sum / n_genes, and no densification is needed until the
        # final centering step.
        u.data = np.log1p(u.data)
        per_cell_mean = np.asarray(u.sum(axis=1)).ravel() / n_genes
        dense_log = u.toarray()
    else:
        u = np.asarray(u)
        dense_log = np.log1p(u)
        per_cell_mean = dense_log.mean(axis=1)

    # Additive centering (the CLR step): map each cell onto the zero-sum
    # (Aitchison) hyperplane. This is what fills the zeros and densifies the matrix.
    return dense_log - per_cell_mean[:, None]


def _normalize_clr_helper(
    x: np.ndarray | CSBase | DaskArray,
    *,
    target_sum: float | None,
    alpha: float | Literal["auto"] | None,
) -> tuple[np.ndarray | DaskArray, np.ndarray | DaskArray]:
    """Compute the shifted CLR (PFlog1pPF) transform backend logic.

    This function acts as the internal processing engine for `normalize_clr`,
    adapted from the reference implementations provided by the Pachter Lab
    (Booeshaghi et al. 2022). It applies initial depth-scaling followed by
    an optimized sparse log-shift operation before mapping coordinates onto
    the zero-sum Aitchison hyperplane.

    Returns the dense, per-cell-centered matrix and the raw per-cell depths
    (the latter is used by the caller to warn about empty cells).

    See `normalize_clr` for the meaning of the parameters.
    """
    # Keep the depths lazy for dask; `.ravel()` would otherwise materialize them.
    cell_depths = stats.sum(x, axis=1)
    if not isinstance(x, DaskArray):
        cell_depths = np.asarray(cell_depths).ravel()

    # Resolve `target_sum` (the paper's PF target constant K, Booeshaghi et al.
    # 2022): the depth the first proportional-fitting step scales every cell to.
    if alpha is not None:
        if alpha == "auto":
            # Estimate the overdispersion from the raw counts (closed-form OLS,
            # as in the reference implementation).
            alpha = _estimate_overdispersion(x)
        if not alpha > 0:
            msg = (
                f"`alpha` must be positive to derive K = 4 * alpha * s, got {alpha}. "
                "The data may be underdispersed; pass `target_sum` explicitly instead."
            )
            raise ValueError(msg)
        # Delta-method calibration: K = 4 * alpha * s, where s is the mean cell
        # depth. This sets the count-scale pseudocount to the variance-stabilizing
        # value y0 = 1 / (4 * alpha), and overrides any value passed as `target_sum`.
        target_sum = 4.0 * alpha * float(cell_depths.mean())
    elif target_sum is None:
        # `float(...)` triggers a single blocking `.compute()` for dask input.
        target_sum = float(cell_depths.mean())
    # else: use the `target_sum` passed by the caller.

    # Proportional fitting: rescale each cell so its counts sum to `target_sum` (K).
    # Dividing by `cell_depths / target_sum` (with allow_divide_by_zero=False)
    # leaves empty cells as all-zero rows instead of producing invalid infinities.
    u = axis_mul_or_truediv(
        x, cell_depths / target_sum, op=truediv, axis=0, allow_divide_by_zero=False
    )

    if isinstance(x, DaskArray):
        # Per-cell centering needs the whole row, so collapse the genes into a
        # single chunk; the transform is then embarrassingly parallel over the
        # row-chunks. The output is always dense, hence the dense numpy `meta`.
        u = u.rechunk({1: -1})
        dense_log = u.map_blocks(
            _clr_log_center, dtype=np.float64, meta=np.array([], dtype=np.float64)
        )
    else:
        dense_log = _clr_log_center(u)

    return dense_log, cell_depths


def normalize_clr(
    adata: AnnData,
    *,
    target_sum: float | None = None,
    alpha: float | Literal["auto"] | None = None,
    layer: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | dict[str, np.ndarray] | None:
    r"""Normalize counts with the shifted centered log-ratio (PFlog1pPF) transform.

    Computes the shifted centered log-ratio (CLR) transform

    .. math::
        T(x)_i = \log(u_i + 1) - \frac{1}{D} \sum_{j=1}^D \log(u_j + 1),

    where :math:`u_i = K \, x_i / \sum_j x_j` are the depth-normalized counts
    (proportional fitting to a target depth :math:`K`) and :math:`D` is the
    number of genes. Equivalently this is proportional fitting, then
    :obj:`~numpy.log1p`, then per-cell mean-centering in log space (the
    centered-log-ratio step). This count transform is simultaneously
    variance-stabilizing, depth-invariant, and rank-preserving :cite:p:`Booeshaghi2022`.

    Because the per-cell centering fills the zero entries, the output is always
    dense, and each cell sums to exactly zero.

    .. note::
        When used with a :class:`~dask.array.Array` in `adata.X`, deriving the
        proportional-fitting target :math:`K` from the data requires a global
        reduction and therefore triggers a blocking `.compute()`: once for the
        default mean-depth target, and once more for `alpha` (including
        ``alpha="auto"``). Only the scalar reduction is materialized, not the
        matrix. Passing `target_sum` explicitly avoids this and keeps the whole
        operation lazy.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    target_sum
        Target depth :math:`K` for the first proportional-fitting step. If `None`
        (and `alpha` is not given), the empirical mean cell depth is used. Unlike
        :func:`~scanpy.pp.normalize_total`, this is only an *intermediate* target:
        cells sum to `target_sum` after proportional fitting, but the subsequent
        log and centering steps make each cell sum to exactly zero in the output.
    alpha
        Negative-binomial overdispersion of the dataset (``var = μ + α·μ²``). When
        given, it overrides `target_sum` and sets :math:`K = 4 \cdot α \cdot s` by
        the delta method, where :math:`s` is the mean cell depth, calibrating the
        count-scale pseudocount to the variance-stabilizing value ``y0 = 1/(4·α)``
        :cite:p:`Booeshaghi2022`. Pass ``"auto"`` to estimate :math:`α` from the
        data (closed-form least squares of ``var = μ + α·μ²`` across genes). A
        :class:`ValueError` is raised if the estimated or supplied :math:`α` is
        not positive (e.g. underdispersed data); pass `target_sum` instead.
    layer
        Layer to normalize instead of `X`.
    inplace
        Whether to update `adata` or return a dictionary with the normalized
        matrix.
    copy
        Whether to modify a copied input object. Not compatible with
        `inplace=False`.

    Returns
    -------
    Returns a dictionary with the normalized copy of `adata.X` (key `"X"`) or
    updates `adata` with the normalized matrix, depending on `inplace`.

    Example
    -------
    >>> import numpy as np
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]], dtype="float32"))
    >>> sc.pp.normalize_clr(adata)
    >>> np.allclose(adata.X.sum(axis=1), 0, atol=1e-6)
    True
    """
    if copy:
        if not inplace:
            msg = "`copy=True` cannot be used with `inplace=False`."
            raise ValueError(msg)
        adata = adata.copy()

    view_to_actual(adata)

    x = _get_obs_rep(adata, layer=layer)
    if isinstance(x, CSCBase):
        x = x.tocsr()
    if issubclass(x.dtype.type, int | np.integer):
        x = x.astype(np.float64)

    start = logg.info("normalizing counts per cell via shifted CLR")

    x, cell_depths = _normalize_clr_helper(x, target_sum=target_sum, alpha=alpha)

    if not isinstance(cell_depths, DaskArray) and not np.all(cell_depths > 0):
        warn("Some cells have zero counts", UserWarning)

    dat = dict(X=x)
    if inplace:
        _set_obs_rep(adata, dat["X"], layer=layer)

    logg.info(
        "    finished ({time_passed})",
        time=start,
    )

    if copy:
        return adata
    elif not inplace:
        return dat
    return None
