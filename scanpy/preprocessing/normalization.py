import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs


def _normalize_data(X, after=None, counts):
    after = np.median(counts) if after is None else after
    counts += (counts == 0)
    counts /= after
    if issparse(X):
        X = sparsefuncs.inplace_row_scale(X, 1/counts)
    else:
        X = X/counts
    return X

def normalize_quantile(data, counts_per_cell_after=None, quantile=1, min_counts=1,
                       key_n_counts=None, copy=False, layers=[], use_rep=None):

    if quantile < 0 or quantile > 1:
        raise ValueError('Choose quantile between 0 and 1.')

    adata = data.copy() if copy else data

    X = adata.X

    counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
    gene_subset = (X>=counts_per_cell[:, None]*quantile).sum(0).A1 == X.shape[0]

    adata._inplace_subset_var(gene_subset)

    X = adata.X

    
