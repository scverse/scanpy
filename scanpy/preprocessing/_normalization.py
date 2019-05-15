import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from .. import logging as logg


def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    after = np.median(counts[counts>0]) if after is None else after
    counts += (counts == 0)
    counts /= after
    if issparse(X):
        sparsefuncs.inplace_row_scale(X, 1/counts)
    else:
        X /= counts[:, None]
    return X if copy else None


def normalize_total(
        adata,
        target_sum=None,
        exclude_highly_expressed=False,
        max_fraction=0.05,
        key_added=None,
        layers=None,
        layer_norm=None,
        inplace=True):
    """\
    Normalize counts per cell.

    If choosing `target_sum=1e6`, this is CPM normalization.

    If `exclude_highly_expressed=True`, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes [Weinreb17]_.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Params
    ------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    target_sum : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each observation (cell) has a total count
        equal to the median of total counts for observations (cells)
        before normalization.
    exclude_highly_expressed : `bool`, optional (default: `False`)
        Exclude (very) highly expressed genes for the computation of the
        normalization factor (size factor) for each cell. A gene is considered
        highly expressed, if it has more than `max_fraction` of the total counts
        in at least one cell. The not-excluded genes will sum up to
        `target_sum`.
    max_fraction : `float`, optional (default: 0.05)
        If `exclude_highly_expressed=True`, consider cells as highly expressed
        that have more counts than `max_fraction` of the original total counts
        in at least one cell.
    key_added : `str`, optional (default: `None`)
        Name of the field in `adata.obs` where the normalization factor is
        stored.
    layers : `str` or list of `str`, optional (default: `None`)
        List of layers to normalize. Set to `'all'` to normalize all layers.
    layer_norm : `str` or `None`, optional (default: `None`)
        Specifies how to normalize layers:

        * If `None`, after normalization, for each layer in *layers* each cell\
        has a total count equal to the median of the *counts_per_cell* before\
        normalization of the layer.

        * If `'after'`, for each layer in *layers* each cell has\
        a total count equal to `target_sum`.

        * If `'X'`, for each layer in *layers* each cell has a total count equal\
        to the median of total counts for observations (cells) of `adata.X`\
        before normalization.

    inplace : `bool`, optional (default: `True`)
        Whether to update `adata` or return dictionary with normalized copies of
        `adata.X` and `adata.layers`.\

    Returns
    -------
    Returns dictionary with normalized copies of `adata.X` and `adata.layers`
    or updates `adata` with normalized version of the original
    `adata.X` and `adata.layers`, depending on `inplace`.

    Example
    --------
    >>> sc.settings.verbosity = 2
    >>> np.set_printoptions(precision=2)
    >>> adata = sc.AnnData(np.array([[3, 3, 3, 6, 6], [1, 1, 1, 2, 2], [1, 22, 1, 2, 2]]))
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)['X']
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=True, max_fraction=0.2, inplace=False)['X']
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)
    """
    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('Choose max_fraction between 0 and 1.')

    X = adata.X

    if not inplace:
        dat = {}  # not recarray because need to support sparse

    gene_subset = None
    if exclude_highly_expressed:
        counts_per_cell = X.sum(1)  # original counts per cell
        counts_per_cell = np.ravel(counts_per_cell)

        # at least one cell as more than max_fraction of counts per cell
        gene_subset = (X>counts_per_cell[:, None]*max_fraction).sum(0)
        gene_subset = (np.ravel(gene_subset) == 0)
        logg.info(
            'The following highly-expressed genes are not considered during normalization factor computation:\n{}'
            .format(adata.var_names[~gene_subset].tolist()))

    # counts per cell for subset, if max_fraction!=1
    X = X if gene_subset is None else adata[:, gene_subset].X
    counts_per_cell = X.sum(1)
    # get rid of adata view
    counts_per_cell = np.ravel(counts_per_cell).copy()

    if key_added is not None:
        adata.obs[key_added] = counts_per_cell

    cell_subset = counts_per_cell>0
    if not np.all(cell_subset):
        logg.warn('Some cells have total count of genes equal to zero')

    if layer_norm == 'after':
        after = target_sum
    elif layer_norm == 'X':
        after = np.median(counts_per_cell[cell_subset])
    elif layer_norm is None:
        after = None
    else:
        raise ValueError('layer_norm should be "after", "X" or None')
    del cell_subset

    if inplace:
        if hasattr(adata.X, '__itruediv__'):
            _normalize_data(adata.X, counts_per_cell, target_sum)
        else:
            adata.X = _normalize_data(adata.X, counts_per_cell, target_sum, copy=True)
    else:
        dat['X'] = _normalize_data(adata.X, counts_per_cell, target_sum, copy=True)

    layers = adata.layers.keys() if layers == 'all' else layers
    if layers is not None:
        for layer in layers:
            L = adata.layers[layer]
            counts = np.ravel(L.sum(1))
            if inplace:
                if hasattr(L, '__itruediv__'):
                    _normalize_data(L, counts, after)
                else:
                    adata.layers[layer] = _normalize_data(L, counts, after, copy=True)
            else:
                dat[layer] = _normalize_data(L, counts, after, copy=True)

    logg.msg('    finished', t=True, end=': ')
    logg.msg('normalized adata.X')
    if key_added is not None:
        logg.msg('and added \'{}\', counts per cell before normalization (adata.obs)'
            .format(key_added))

    return dat if not inplace else None
