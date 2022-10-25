from typing import Optional, Union, Iterable, Dict
from warnings import warn

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs

try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        pass


from scanpy import logging as logg
from scanpy._compat import Literal
from scanpy._utils import view_to_actual

from scanpy.get import _get_obs_rep, _set_obs_rep


def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    if issubclass(X.dtype.type, (int, np.integer)):
        X = X.astype(np.float32)  # TODO: Check if float64 should be used
    if isinstance(counts, DaskArray):
        counts_greater_than_zero = counts[counts > 0].compute_chunk_sizes()
    else:
        counts_greater_than_zero = counts[counts > 0]

    after = np.median(counts_greater_than_zero, axis=0) if after is None else after
    counts += counts == 0
    counts = counts / after
    if issparse(X):
        sparsefuncs.inplace_row_scale(X, 1 / counts)
    elif isinstance(counts, np.ndarray):
        np.divide(X, counts[:, None], out=X)
    else:
        X = np.divide(X, counts[:, None])  # dask does not support kwarg "out"
    return X


def normalize_total(
    adata: AnnData,
    target_sum: Optional[float] = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    layers: Union[Literal['all'], Iterable[str]] = None,
    layer_norm: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
) -> Optional[Dict[str, np.ndarray]]:
    """\
    Normalize counts per cell.

    Normalize each cell by total counts over all genes,
    so that every cell has the same total count after normalization.
    If choosing `target_sum=1e6`, this is CPM normalization.

    If `exclude_highly_expressed=True`, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes [Weinreb17]_.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

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
        `target_sum`.
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
    --------
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> sc.settings.verbosity = 2
    >>> np.set_printoptions(precision=2)
    >>> adata = AnnData(np.array([
    ...    [3, 3, 3, 6, 6],
    ...    [1, 1, 1, 2, 2],
    ...    [1, 22, 1, 2, 2],
    ... ]))
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)['X']
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(
    ...     adata, target_sum=1, exclude_highly_expressed=True,
    ...     max_fraction=0.2, inplace=False
    ... )['X']
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)
    """
    if copy:
        if not inplace:
            raise ValueError("`copy=True` cannot be used with `inplace=False`.")
        adata = adata.copy()

    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('Choose max_fraction between 0 and 1.')

    # Deprecated features
    if layers is not None:
        warn(
            FutureWarning(
                "The `layers` argument is deprecated. Instead, specify individual "
                "layers to normalize with `layer`."
            )
        )
    if layer_norm is not None:
        warn(
            FutureWarning(
                "The `layer_norm` argument is deprecated. Specify the target size "
                "factor directly with `target_sum`."
            )
        )

    if layers == 'all':
        layers = adata.layers.keys()
    elif isinstance(layers, str):
        raise ValueError(
            f"`layers` needs to be a list of strings or 'all', not {layers!r}"
        )

    view_to_actual(adata)

    X = _get_obs_rep(adata, layer=layer)

    gene_subset = None
    msg = 'normalizing counts per cell'
    if exclude_highly_expressed:
        counts_per_cell = X.sum(1)  # original counts per cell
        counts_per_cell = np.ravel(counts_per_cell)

        # at least one cell as more than max_fraction of counts per cell

        gene_subset = (X > counts_per_cell[:, None] * max_fraction).sum(0)
        gene_subset = np.ravel(gene_subset) == 0

        msg += (
            ' The following highly-expressed genes are not considered during '
            f'normalization factor computation:\n{adata.var_names[~gene_subset].tolist()}'
        )
        counts_per_cell = X[:, gene_subset].sum(1)
    else:
        counts_per_cell = X.sum(1)
    start = logg.info(msg)
    counts_per_cell = np.ravel(counts_per_cell)

    cell_subset = counts_per_cell > 0
    if not np.all(cell_subset):
        warn(UserWarning('Some cells have zero counts'))

    if inplace:
        if key_added is not None:
            adata.obs[key_added] = counts_per_cell
        _set_obs_rep(
            adata, _normalize_data(X, counts_per_cell, target_sum), layer=layer
        )
    else:
        # not recarray because need to support sparse
        dat = dict(
            X=_normalize_data(X, counts_per_cell, target_sum, copy=True),
            norm_factor=counts_per_cell,
        )

    # Deprecated features
    if layer_norm == 'after':
        after = target_sum
    elif layer_norm == 'X':
        after = np.median(counts_per_cell[cell_subset])
    elif layer_norm is None:
        after = None
    else:
        raise ValueError('layer_norm should be "after", "X" or None')

    for layer_to_norm in layers if layers is not None else ():
        res = normalize_total(
            adata, layer=layer_to_norm, target_sum=after, inplace=inplace
        )
        if not inplace:
            dat[layer_to_norm] = res["X"]

    logg.info(
        '    finished ({time_passed})',
        time=start,
    )
    if key_added is not None:
        logg.debug(
            f'and added {key_added!r}, counts per cell before normalization (adata.obs)'
        )

    if copy:
        return adata
    elif not inplace:
        return dat
