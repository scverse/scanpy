import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs


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

def normalize_quantile(data, counts_per_cell_after=None, quantile=1, min_counts=1,
                       key_n_counts=None, copy=False, layers=[], use_rep=None):

    if quantile < 0 or quantile > 1:
        raise ValueError('Choose quantile between 0 and 1.')

    if key_n_counts is None: key_n_counts = 'n_counts'

    adata = data.copy() if copy else data

    if quantile < 1:
        X = adata.X
        counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
        gene_subset = (X>counts_per_cell[:, None]*quantile).sum(0)
        gene_subset = gene_subset.A1 == 0 if issparse(X) else gene_subset == 0
        adata._inplace_subset_var(gene_subset)

    X = adata.X
    counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
    adata.obs[key_n_counts] = counts_per_cell
    cell_subset = counts_per_cell >= min_counts
    adata._inplace_subset_obs(cell_subset)
    del X
    counts_per_cell = counts_per_cell[cell_subset]
    _normalize_data(adata.X, counts_per_cell, counts_per_cell_after)

    layers = adata.layers.keys() if layers == 'all' else layers
    if use_rep == 'after':
        after = counts_per_cell_after
    elif use_rep == 'X':
        after = np.median(counts_per_cell)
    elif use_rep is None:
        after = None

    for layer in layers:
        L = adata.layers[layer]
        counts = L.sum(1).A1 if issparse(L) else L.sum(1)
        _normalize_data(L, counts, after)

    return adata if copy else None

def normalize_total(data, counts_per_cell_after=None, counts_per_cell=None,
                    key_n_counts=None, copy=False, layers=[], use_rep=None,
                    min_counts=1):
    return normalize_quantile(data=data, counts_per_cell_after=counts_per_cell_after,
                              key_n_counts=key_n_counts, copy=copy, layers=layers,
                              use_rep=use_rep, min_counts=min_counts, quantile=1)
