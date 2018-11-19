import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from .. import logging as logg

def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    after = np.median(counts) if after is None else after
    counts += (counts == 0)
    counts /= after
    if issparse(X):
        X = sparsefuncs.inplace_row_scale(X, 1/counts)
    else:
        X /= counts[:, None]
    return X if copy else None

def normalize_quantile(data, counts_per_cell_after=None, counts_per_cell=None,
                       quantile=1, min_counts=1, key_n_counts=None, copy=False,
                       layers=[], use_rep=None):
    """Normalize total counts per cell.

    Normalize each cell by total counts over genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    quantile : `float`, optional (default: 1)
        Only use genes are less than fraction (specified by *quantile*) of the total
        reads in every cell.
    min_counts : `int`, optional (default: 1)
        Cells with counts less than `min_counts` are filtered out during
        normalization.
    key_n_counts : `str`, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    layers : `str` or list of `str`, optional (default: `[]`)
        List of layers to normalize. Set to `'all'` to normalize all layers.
    use_rep : `str` or `None`, optional (default: `None`)
        Specifies how to normalize layers.
        If `None`, after normalization, for each layer in *layers* each cell has a total count equal
        to the median of the *counts_per_cell* before normalization of the layer.
        If `'after'`, for each layer in *layers* each cell has a total count equal
        to counts_per_cell_after.
        If `'X'`, for each layer in *layers* each cell has a total count equal
        to the median of the *counts_per_cell* of data.X before normalization.

    Returns
    -------
    Returns or updates `adata` with normalized version of the original
    `adata.X`, depending on `copy`.

    Examples
    --------
    >>> adata = AnnData(np.array([[1, 0, 1], [3, 0, 1], [5, 6, 1]]))
    >>> sc.pp.normalize_quantile(adata, quantile=0.7)
    >>> print(adata.X)
    [[0.         1.        ]
     [0.         1.        ]
     [0.85714287 0.14285715]]
    """
    if quantile < 0 or quantile > 1:
        raise ValueError('Choose quantile between 0 and 1.')

    if key_n_counts is None: key_n_counts = 'n_counts'

    adata = data.copy() if copy else data

    if quantile < 1:
        logg.msg('normalizing by count per cell for \
                  genes that make up less than quantile * total count per cell', r=True)
        X = adata.X

        counts_per_cell = counts_per_cell if counts_per_cell is not None else X.sum(1)
        counts_per_cell = np.ravel(counts_per_cell)

        gene_subset = (X>counts_per_cell[:, None]*quantile).sum(0)
        gene_subset = (np.ravel(gene_subset) == 0)
        adata._inplace_subset_var(gene_subset)
        del gene_subset
    else:
        logg.msg('normalizing by total count per cell', r=True)

    X = adata.X

    counts_per_cell = counts_per_cell if (counts_per_cell is not None and quantile == 1) else X.sum(1)
    counts_per_cell = np.ravel(counts_per_cell)

    adata.obs[key_n_counts] = counts_per_cell
    cell_subset = counts_per_cell >= min_counts
    adata._inplace_subset_obs(cell_subset)
    del X
    counts_per_cell = counts_per_cell[cell_subset]

    if use_rep == 'after':
        after = counts_per_cell_after
    elif use_rep == 'X':
        after = np.median(counts_per_cell)
    elif use_rep is None:
        after = None
    else: raise ValueError('use_rep should be "after", "X" or None')

    _normalize_data(adata.X, counts_per_cell, counts_per_cell_after)

    layers = adata.layers.keys() if layers == 'all' else layers
    for layer in layers:
        L = adata.layers[layer]
        counts = np.ravel(L.sum(1))
        _normalize_data(L, counts, after)

    logg.msg('    finished', t=True, end=': ')
    logg.msg('normalized adata.X and added', no_indent=True)
    logg.msg('    \'{}\', counts per cell before normalization (adata.obs)'
        .format(key_n_counts))

    return adata if copy else None

def normalize_total(data, counts_per_cell_after=None, counts_per_cell=None,
                    key_n_counts=None, copy=False, layers=[], use_rep=None,
                    min_counts=1):
    """Normalize total counts per cell.

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    key_n_counts : `str`, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    layers : `str` or list of `str`, optional (default: `[]`)
        List of layers to normalize. Set to `'all'` to normalize all layers.
    use_rep : `str` or `None`, optional (default: `None`)
        Specifies how to normalize layers.
        If `None`, after normalization, for each layer in *layers* each cell has a total count equal
        to the median of the *counts_per_cell* before normalization of the layer.
        If `'after'`, for each layer in *layers* each cell has a total count equal
        to counts_per_cell_after.
        If `'X'`, for each layer in *layers* each cell has a total count equal
        to the median of the *counts_per_cell* of data.X before normalization.
    min_counts : `int`, optional (default: 1)
        Cells with counts less than `min_counts` are filtered out during
        normalization.

    Returns
    -------
    Returns or updates `adata` with normalized version of the original
    `adata.X`, depending on `copy`.

    Examples
    --------
    >>> adata = AnnData(data=np.array([[1, 0], [3, 0], [5, 6]]))
    >>> print(adata.X.sum(axis=1))
    [  1.   3.  11.]
    >>> sc.pp.normalize_per_cell(adata)
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts
    0       1.0
    1       3.0
    2      11.0
    [ 3.  3.  3.]
    >>> sc.pp.normalize_per_cell(adata, counts_per_cell_after=1,
    >>>                          key_n_counts='n_counts2')
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts  n_counts2
    0       1.0        3.0
    1       3.0        3.0
    2      11.0        3.0
    [ 1.  1.  1.]
    """
    return normalize_quantile(data=data, counts_per_cell_after=counts_per_cell_after,
                              counts_per_cell=counts_per_cell, key_n_counts=key_n_counts,
                              copy=copy, layers=layers, use_rep=use_rep, min_counts=min_counts,
                              quantile=1)
