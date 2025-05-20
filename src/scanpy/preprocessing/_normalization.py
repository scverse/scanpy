from __future__ import annotations

from operator import truediv
from typing import TYPE_CHECKING
from warnings import warn

import numba
import numpy as np

from .. import logging as logg
from .._compat import CSBase, CSCBase, CSRBase, DaskArray, njit, old_positionals
from .._utils import axis_mul_or_truediv, axis_sum, dematrix, view_to_actual
from ..get import _get_obs_rep, _set_obs_rep

try:
    import dask
    import dask.array as da
except ImportError:
    da = None
    dask = None

if TYPE_CHECKING:
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
    indptr,
    indices,
    data,
    *,
    rows,
    columns,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    n_threads: int = 10,
):
    """For sparse CSR matrix, compute the normalization factors."""
    counts_per_cell = np.zeros(rows, dtype=data.dtype)
    for i in numba.prange(rows):
        count = 0.0
        for j in range(indptr[i], indptr[i + 1]):
            count += data[j]
        counts_per_cell[i] = count
    if exclude_highly_expressed:
        counts_per_cols_t = np.zeros((n_threads, columns), dtype=np.int32)
        counts_per_cols = np.zeros(columns, dtype=np.int32)

        for i in numba.prange(n_threads):
            for r in range(i, rows, n_threads):
                for j in range(indptr[r], indptr[r + 1]):
                    if data[j] > max_fraction * counts_per_cell[r]:
                        minor_index = indices[j]
                        counts_per_cols_t[i, minor_index] += 1
        for c in numba.prange(columns):
            counts_per_cols[c] = counts_per_cols_t[:, c].sum()

        for i in numba.prange(rows):
            count = 0.0
            for j in range(indptr[i], indptr[i + 1]):
                if counts_per_cols[indices[j]] == 0:
                    count += data[j]
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
            x.indptr,
            x.indices,
            x.data,
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
        counts_per_cell = axis_sum(x, axis=1)
        if exclude_highly_expressed:
            # at least one cell as more than max_fraction of counts per cell
            hi_exp = dematrix(x > counts_per_cell[:, None] * max_fraction)
            gene_subset = axis_sum(hi_exp, axis=0) == 0

            counts_per_cell = axis_sum(x[:, gene_subset], axis=1)
        if target_sum is None:
            target_sum = _compute_nnz_median(counts_per_cell)

    counts_per_cell = counts_per_cell / target_sum
    out = x if isinstance(x, np.ndarray | CSBase) else None
    X = axis_mul_or_truediv(
        x, counts_per_cell, op=truediv, out=out, allow_divide_by_zero=False, axis=0
    )
    return X, counts_per_cell, gene_subset


@old_positionals(
    "target_sum",
    "exclude_highly_expressed",
    "max_fraction",
    "key_added",
    "layer",
    "inplace",
    "copy",
)
def normalize_total(  # noqa: PLR0912
    adata: AnnData,
    *,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: str | None = None,
    layer: str | None = None,
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

    .. note::
        When used with a :class:`~dask.array.Array` in `adata.X`, this function will have to
        call functions that trigger `.compute()` on the :class:`~dask.array.Array` if `exclude_highly_expressed`
        is `True`, `layer_norm` is not `None`, or if `key_added` is not `None`.

    Params
    ------
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
        Layer to normalize instead of `X`. If `None`, `X` is normalized.
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
    >>> np.set_printoptions(precision=2)
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

    x = _get_obs_rep(adata, layer=layer)
    if x is None:
        msg = f"Layer {layer!r} not found in adata."
        raise ValueError(msg)
    if isinstance(x, CSCBase):
        x = x.tocsr()
    if not inplace:
        x = x.copy()
    if issubclass(x.dtype.type, int | np.integer):
        x = x.astype(np.float32)  # TODO: Check if float64 should be used

    start = logg.info("normalizing counts per cell")

    X, counts_per_cell, gene_subset = _normalize_total_helper(
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
        warn("Some cells have zero counts", UserWarning, stacklevel=2)

    dat = dict(
        X=X,
        norm_factor=counts_per_cell,
    )
    if inplace:
        if key_added is not None:
            adata.obs[key_added] = dat["norm_factor"]
        _set_obs_rep(adata, dat["X"], layer=layer)

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
