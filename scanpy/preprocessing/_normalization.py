import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from .. import logging as logg
from ..utils import doc_params

doc_norm_descr = """\
Normalize counts per cell.

For `fraction=1`, this is standard total-count normalization, if choosing
`target_sum=1e6`, this is CPM normalization.

Normalize each cell by sum of counts over genes that make up less than fraction
(specified by *fraction*) of the total counts in every cell. These genes in each
cell will sum up to *target_sum*.

Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
[Zheng17]_ or SPRING [Weinreb17]_.\
"""

doc_params_bulk = """\
Parameters
----------
adata : :class:`~anndata.AnnData`
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
    to cells and columns to genes.
target_sum : `float` or `None`, optional (default: `None`)
    If `None`, after normalization, each observation (cell) has a total count
    equal to the median of total counts for observations (cells)
    before normalization.
fraction : `float`, optional (default: 1)
    Only use genes that make up less than fraction (specified by *fraction*)
    of the total count in every cell. So only these genes will sum up
    to *target_sum*.
key_added : `str`, optional (default: `None`)
    Name of the field in `adata.obs` where the total counts per cell are
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
"""

doc_norm_return = """\
Returns
-------
Returns dictionary with normalized copies of `adata.X` and `adata.layers`
or updates `adata` with normalized version of the original
`adata.X` and `adata.layers`, depending on `inplace`.
"""

doc_examples = """\
Example
--------
>>> adata = AnnData(np.array([[1, 0], [3, 0], [5, 6]]))
>>> print(adata.X.sum(axis=1))
[  1.   3.  11.]
>>> sc.pp.normalize_total(adata, key_added='n_counts')
>>> print(adata.obs)
>>> print(adata.X.sum(axis=1))
   n_counts
0       1.0
1       3.0
2      11.0
[ 3.  3.  3.]
>>> sc.pp.normalize_total(adata, target_sum=1,
>>>                       key_added='n_counts2')
>>> print(adata.obs)
>>> print(adata.X.sum(axis=1))
   n_counts  n_counts2
0       1.0        3.0
1       3.0        3.0
2      11.0        3.0
[ 1.  1.  1.]

An example using `fraction`.

>>> adata = AnnData(np.array([[1, 0, 1], [3, 0, 1], [5, 6, 1]]))
>>> sc.pp.normalize_total(adata, fraction=0.7)
>>> print(adata.X)
[[1.         0.         1.        ]
 [3.         0.         1.        ]
 [0.71428573 0.85714287 0.14285715]]

Genes 1 and 2 were normalized and now sum up to 1 in each cell.
"""


def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    after = np.median(counts[counts>0]) if after is None else after
    counts += (counts == 0)
    counts /= after
    if issparse(X):
        X = sparsefuncs.inplace_row_scale(X, 1/counts)
    else:
        X /= counts[:, None]
    return X if copy else None


@doc_params(norm_descr=doc_norm_descr, params_bulk=doc_params_bulk, norm_return=doc_norm_return,
            examples=doc_examples)
def normalize_total(adata, target_sum=None, fraction=1, key_added=None,
                    layers=None, layer_norm=None, inplace=True):
    """\
    {norm_descr}

    {params_bulk}

    {norm_return}

    {examples}
    """
    if fraction < 0 or fraction > 1:
        raise ValueError('Choose fraction between 0 and 1.')

    X = adata.X
    gene_subset = None
    if not inplace:
    # not recarray because need to support sparse
        dat = {}

    if fraction < 1:
        logg.msg('normalizing by count per cell for \
                  genes that make up less than fraction * total count per cell', r=True)
        X = adata.X

        counts_per_cell = X.sum(1)
        counts_per_cell = np.ravel(counts_per_cell)

        gene_subset = (X>counts_per_cell[:, None]*fraction).sum(0)
        gene_subset = (np.ravel(gene_subset) == 0)
    else:
        logg.msg('normalizing by total count per cell', r=True)

    X = X if gene_subset is None else adata[:, gene_subset].X
    counts_per_cell = X.sum(1)
    # get rid of adata view
    counts_per_cell = np.ravel(counts_per_cell).copy()
    del X
    del gene_subset

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
