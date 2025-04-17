from __future__ import annotations

from operator import truediv
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
from fast_array_utils import stats

from .. import logging as logg
from .._compat import CSBase, DaskArray, old_positionals
from .._utils import axis_mul_or_truediv, view_to_actual
from ..get import _get_obs_rep, _set_obs_rep

try:
    import dask
    import dask.array as da
except ImportError:
    da = None
    dask = None

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from anndata import AnnData


def _compute_nnz_median(counts: np.ndarray | DaskArray) -> np.floating:
    """Given a 1D array of counts, compute the median of the non-zero counts."""
    if isinstance(counts, DaskArray):
        counts = counts.compute()
    counts_greater_than_zero = counts[counts > 0]
    median = np.median(counts_greater_than_zero)
    return median


def _normalize_data(
    X: np.ndarray | CSBase,
    counts: np.ndarray,
    after: float | None = None,
    *,
    copy: bool = False,
) -> np.ndarray:
    if counts.ndim != 1:
        msg = "counts must be a 1D array"
        raise ValueError(msg)
    X = X.copy() if copy else X
    if issubclass(X.dtype.type, int | np.integer):
        X = X.astype(np.float32)  # TODO: Check if float64 should be used
    if after is None:
        after = _compute_nnz_median(counts)
    counts = counts / after
    out = X if isinstance(X, np.ndarray | CSBase) else None
    return axis_mul_or_truediv(
        X, counts, op=truediv, out=out, allow_divide_by_zero=False, axis=0
    )


@old_positionals(
    "target_sum",
    "exclude_highly_expressed",
    "max_fraction",
    "key_added",
    "layer",
    "layers",
    "layer_norm",
    "inplace",
    "copy",
)
def normalize_total(  # noqa: PLR0912, PLR0915
    adata: AnnData,
    *,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: str | None = None,
    layer: str | None = None,
    layers: Literal["all"] | Iterable[str] | None = None,
    layer_norm: str | None = None,
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
    normalizing counts per cell. The following highly-expressed genes are not considered during normalization factor computation:
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

    # Deprecated features
    if layers is not None:
        warn(
            "The `layers` argument is deprecated. Instead, specify individual "
            "layers to normalize with `layer`.",
            FutureWarning,
            stacklevel=2,
        )
    if layer_norm is not None:
        warn(
            "The `layer_norm` argument is deprecated. Specify the target size "
            "factor directly with `target_sum`.",
            FutureWarning,
            stacklevel=2,
        )

    if layers == "all":
        layers = adata.layers.keys()
    elif isinstance(layers, str):
        msg = f"`layers` needs to be a list of strings or 'all', not {layers!r}"
        raise ValueError(msg)

    view_to_actual(adata)

    x = _get_obs_rep(adata, layer=layer)

    gene_subset = None
    msg = "normalizing counts per cell"

    counts_per_cell = stats.sum(x, axis=1)
    assert counts_per_cell.ndim == 1
    if exclude_highly_expressed:
        # at least one cell as more than max_fraction of counts per cell

        x = x > counts_per_cell[:, None] * max_fraction
        if isinstance(x, np.matrix):
            x = x.A
        elif isinstance(x, DaskArray) and isinstance(x._meta, np.matrix):
            x = x.map_blocks(np.asarray, meta=np.array([], dtype=x.dtype))
        gene_subset = stats.sum(x, axis=0) == 0
        assert gene_subset.ndim == 1

        msg += (
            ". The following highly-expressed genes are not considered during "
            f"normalization factor computation:\n{adata.var_names[~gene_subset].tolist()}"
        )
        counts_per_cell = stats.sum(x[:, gene_subset], axis=1)
        assert counts_per_cell.ndim == 1

    start = logg.info(msg)

    cell_subset = counts_per_cell > 0
    if not isinstance(cell_subset, DaskArray) and not np.all(cell_subset):
        warn("Some cells have zero counts", UserWarning, stacklevel=2)

    if inplace:
        if key_added is not None:
            adata.obs[key_added] = counts_per_cell
        _set_obs_rep(
            adata, _normalize_data(x, counts_per_cell, target_sum), layer=layer
        )
    else:
        # not recarray because need to support sparse
        dat = dict(
            X=_normalize_data(x, counts_per_cell, target_sum, copy=True),
            norm_factor=counts_per_cell,
        )

    # Deprecated features
    if layer_norm == "after":
        after = target_sum
    elif layer_norm == "X":
        after = np.median(counts_per_cell[cell_subset])
    elif layer_norm is None:
        after = None
    else:
        msg = 'layer_norm should be "after", "X" or None'
        raise ValueError(msg)

    for layer_to_norm in layers if layers is not None else ():
        res = normalize_total(
            adata, layer=layer_to_norm, target_sum=after, inplace=inplace
        )
        if not inplace:
            dat[layer_to_norm] = res["X"]

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
