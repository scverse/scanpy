from __future__ import annotations

from operator import truediv
from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils import stats
from fast_array_utils.numba import njit
from fast_array_utils.stats import mean_var
from scipy import sparse

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
            "Pass a positive `alpha` explicitly."
        )
        raise ValueError(msg)
    alpha = float(numerator / denominator)
    if not alpha > 0:
        msg = (
            f"Estimated overdispersion is non-positive (alpha = {alpha}); "
            "pass a positive `alpha` explicitly."
        )
        raise ValueError(msg)
    return alpha


def _log1p_sparse_block(x: np.ndarray | CSBase) -> np.ndarray | CSBase:
    """Apply log1p to a dense or sparse block while preserving sparse zeros."""
    if isinstance(x, CSBase):
        x = x.copy()
        x.data = np.log1p(x.data)
        return x
    return np.log1p(x)


def _densify_shifted_clr(
    log_values: np.ndarray | CSBase | DaskArray, row_center: np.ndarray | DaskArray
) -> np.ndarray | DaskArray:
    dense_log = log_values.toarray() if isinstance(log_values, CSBase) else log_values
    return dense_log - row_center[:, None]


def _normalize_clr_helper(  # noqa: PLR0912
    x: np.ndarray | CSBase | DaskArray,
    *,
    target: Literal["auto", "mean", "median"] | float,
    alpha: float | None,
) -> tuple[
    np.ndarray | CSBase | DaskArray,
    np.ndarray | DaskArray,
    np.ndarray | DaskArray,
    dict[str, object],
]:
    """Compute sparse PFlog / shifted-CLR values and row centers."""
    # Keep the depths lazy for dask; `.ravel()` would otherwise materialize them.
    cell_depths = stats.sum(x, axis=1)
    if not isinstance(x, DaskArray):
        cell_depths = np.asarray(cell_depths).ravel()

    mean_depth = (
        float(cell_depths.mean().compute())
        if isinstance(cell_depths, DaskArray)
        else float(cell_depths.mean())
    )
    if alpha is not None:
        if not alpha > 0:
            msg = (
                f"`alpha` must be positive to compute PFlog, got {alpha}. "
                "The data may be underdispersed."
            )
            raise ValueError(msg)
        scale: float | np.ndarray | DaskArray = 4.0 * float(alpha)
        target_sum = 4.0 * float(alpha) * mean_depth
        resolved_target = "alpha"
    elif target == "auto":
        alpha = _estimate_overdispersion(x)
        scale = 4.0 * alpha
        target_sum = 4.0 * alpha * mean_depth
        resolved_target = "auto"
    else:
        alpha = None
        if target == "mean":
            target_sum = mean_depth
        elif target == "median":
            depths = (
                cell_depths.compute()
                if isinstance(cell_depths, DaskArray)
                else cell_depths
            )
            target_sum = float(np.median(depths))
        elif isinstance(target, bool) or not isinstance(target, int | float):
            msg = "`target` must be 'auto', 'mean', 'median', or a positive number."
            raise ValueError(msg)
        else:
            target_sum = float(target)
        if not target_sum > 0:
            msg = f"`target` must resolve to a positive depth, got {target_sum}."
            raise ValueError(msg)
        resolved_target = str(target)
        x = axis_mul_or_truediv(
            x, cell_depths / target_sum, op=truediv, axis=0, allow_divide_by_zero=False
        )

    if isinstance(x, DaskArray):
        if alpha is not None:
            x = x * scale
        log_values = x.map_blocks(
            _log1p_sparse_block, dtype=np.float64, meta=x._meta.astype(np.float64)
        )
    else:
        if alpha is not None:
            x = x * scale
        log_values = _log1p_sparse_block(x)

    row_center = stats.sum(log_values, axis=1) / x.shape[1]
    if not isinstance(row_center, DaskArray):
        row_center = np.asarray(row_center).ravel()
    report = dict(
        target=resolved_target,
        alpha=None if alpha is None else float(alpha),
        k=float(target_sum),
        mean_depth=mean_depth,
    )
    return log_values, row_center, cell_depths, report


def normalize_clr(
    adata: AnnData,
    *,
    target: Literal["auto", "mean", "median"] | float = "auto",
    alpha: float | None = None,
    key_added: str = "pflog",
    densify: bool = False,
    layer: str | None = None,
    inplace: bool = True,
    copy: bool = False,
) -> AnnData | dict[str, np.ndarray | CSBase | DaskArray | dict[str, object]] | None:
    r"""Normalize counts with the shifted centered log-ratio (PFlog) transform.

    With `target="auto"` (or explicit `alpha`), computes PFlog as

    .. math::
        T(x)_i = \log(1 + 4 α x_i)
            - \frac{1}{D} \sum_{j=1}^D \log(1 + 4 α x_j),

    which is equivalent to centering
    :math:`\log(x_i + 1 / (4 α))` because the constant :math:`\log(4 α)`
    cancels during CLR centering. For memory efficiency, the sparse log values
    are stored separately from the per-cell centering vector.

    .. note::
        `target="auto"` is the recommended PFlog default. Fixed targets
        (``"mean"``, ``"median"``, or a positive number) apply shifted CLR
        after rescaling each cell to the same total count.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    target
        ``"auto"`` estimates the negative-binomial overdispersion and applies
        PFlog. ``"mean"``, ``"median"``, or a positive numeric value set a
        fixed total count used before applying the shifted CLR transform.
    alpha
        Negative-binomial overdispersion of the dataset (``var = μ + α·μ²``).
        When given, it overrides `target` and applies PFlog with sparse values
        ``log1p(4 * alpha * x)``.
    key_added
        Layer key used to store the sparse log values. The row-centering vector
        is stored in ``adata.obs[f"{key_added}_center"]``.
    densify
        If `True`, also materialize the dense centered matrix into `X` or the
        selected `layer`. This is intended for small data or generic downstream
        tools that cannot consume the sparse-plus-row-center representation.
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
    Returns a dictionary with normalized values and metadata or updates `adata`,
    depending on `inplace`.

    Example
    -------
    >>> import numpy as np
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> adata = AnnData(np.array([[1, 2, 30], [4, 50, 6]], dtype="float32"))
    >>> sc.pp.normalize_clr(adata, alpha=0.5)
    >>> "pflog" in adata.layers and "pflog_center" in adata.obs
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
    if not inplace:
        x = x.copy()
    if issubclass(x.dtype.type, int | np.integer):
        x = x.astype(np.float64)

    start = logg.info("normalizing counts per cell via PFlog")

    log_values, row_center, cell_depths, report = _normalize_clr_helper(
        x, target=target, alpha=alpha
    )

    if not isinstance(cell_depths, DaskArray) and not np.all(cell_depths > 0):
        warn("Some cells have zero counts", UserWarning)

    row_center_obs = (
        np.asarray(row_center.compute()).ravel()
        if isinstance(row_center, DaskArray)
        else np.asarray(row_center).ravel()
    )
    if not isinstance(log_values, DaskArray | CSBase):
        log_values = sparse.csr_matrix(log_values)  # noqa: TID251

    dense_x = _densify_shifted_clr(log_values, row_center) if densify else None
    row_center_key = f"{key_added}_center"
    metadata = dict(
        report,
        method="PFlog",
        encoding_type="shifted_clr",
        row_center_key=row_center_key,
        layer=layer,
        densify=densify,
    )
    dat = dict(
        X=dense_x if densify else log_values,
        row_center=row_center_obs,
        metadata=metadata,
    )
    if inplace:
        adata.layers[key_added] = log_values
        adata.obs[row_center_key] = row_center_obs
        adata.uns[key_added] = metadata
        if densify:
            _set_obs_rep(adata, dense_x, layer=layer)

    logg.info(
        "    finished ({time_passed})",
        time=start,
    )

    if copy:
        return adata
    elif not inplace:
        return dat
    return None
