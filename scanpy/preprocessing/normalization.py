import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs
from .. import logging as logg
from ..utils import doc_params
from .docs import doc_norm_bulk, doc_norm_quant, doc_norm_return, doc_ex_quant, doc_ex_total

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

@doc_params(norm_bulk=doc_norm_bulk, norm_quant=doc_norm_quant, norm_return=doc_norm_return, ex_quant=doc_ex_quant)
def normalize_quantile(adata, target_sum=None, quantile=1, key_added=None,
                       layers=[], layer_norm=None, inplace=True):
    """\
    {norm_bulk}
    {norm_quant}

    {norm_return}

    {ex_quant}
    """
    if quantile < 0 or quantile > 1:
        raise ValueError('Choose quantile between 0 and 1.')

    X = adata.X
    gene_subset = None
    if not inplace:
    # not recarray because need to support sparse
        dat = {}

    if quantile < 1:
        logg.msg('normalizing by count per cell for \
                  genes that make up less than quantile * total count per cell', r=True)
        X = adata.X

        counts_per_cell = X.sum(1)
        counts_per_cell = np.ravel(counts_per_cell)

        gene_subset = (X>counts_per_cell[:, None]*quantile).sum(0)
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
        logg.warning('Some cells have total count of genes equal to zero')

    if layer_norm == 'after':
        after = target_sum
    elif layer_norm == 'X':
        after = np.median(counts_per_cell[cell_subset])
    elif layer_norm is None:
        after = None
    else: raise ValueError('layer_norm should be "after", "X" or None')
    del cell_subset

    if inplace:
        _normalize_data(adata.X, counts_per_cell, target_sum)
    else:
        dat['X'] = _normalize_data(adata.X, counts_per_cell, target_sum, copy=True)

    layers = adata.layers.keys() if layers == 'all' else layers
    for layer in layers:
        L = adata.layers[layer]
        counts = np.ravel(L.sum(1))
        if inplace:
            _normalize_data(L, counts, after)
        else:
            dat[layer] = _normalize_data(L, counts, after, copy=True)

    logg.msg('    finished', t=True, end=': ')
    logg.msg('normalized adata.X')
    if key_added is not None:
        logg.msg('and added \'{}\', counts per cell before normalization (adata.obs)'
            .format(key_added))

    return dat if not inplace else None

@doc_params(norm_bulk=doc_norm_bulk, norm_return=doc_norm_return, ex_total=doc_ex_total)
def normalize_total(adata, target_sum=None, key_added=None, layers=[], layer_norm=None, inplace=True):
    """\
    {norm_bulk}

    {norm_return}

    {ex_total}
    """
    return normalize_quantile(adata=adata, target_sum=target_sum,
                              key_added=key_added, layers=layers,
                              layer_norm=layer_norm, quantile=1, inplace=inplace)
