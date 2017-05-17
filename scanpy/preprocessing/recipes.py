# Author: F. Alex Wolf (http://falexwolf.de)
"""
Preprocessing recipes from the literature
"""

from .. import settings as sett
from . import simple as pp
from .. import plotting as pl

def recipe_weinreb16(adata, mean_threshold=0.01, cv_threshold=2,
              n_pcs=50, svd_solver='randomized', random_state=0, copy=False):
    """
    Normalization and filtering as of Weinreb et al. (2016).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    svd_solver : str, optional (default: 'randomized')
        SVD solver to use. Either 'arpack' for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or 'randomized' for the randomized algorithm
        due to Halko (2009).
    random_state : int, optional (default: 0)
        Change to use different intial states for the optimization.
    copy : bool (default: False)
        Return a copy if true.

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016).
    """
    sett.m(0, 'preprocess: weinreb16, X has shape n_samples x n_variables =',
           adata.X.shape[0], 'x', adata.X.shape[1])
    if copy: adata = adata.copy()
    adata.X = pp.normalize_per_cell_weinreb16(adata.X,
                                           max_fraction=0.05,
                                           mult_with_mean=True)
    # filter out genes with mean expression < 0.1 and coefficient of variance <
    # cv_threshold
    gene_subset = pp.filter_genes_cv_deprecated(adata.X, mean_threshold, cv_threshold)
    # three alternative ways of slicing
    adata.inplace_subset_var(gene_subset)  # this modifies the object itself
    # adata = adata[:, gene_subset]  # this does a copy
    # adata[:, ~gene_subset] = None  # this doesn't work yet
    # compute zscore of filtered matrix and compute PCA
    X_pca = pp.pca(zscore_deprecated(adata.X),
                   n_comps=n_pcs, svd_solver=svd_solver, random_state=random_state)
    # update adata
    adata.smp['X_pca'] = X_pca
    sett.m(0, 'X_pca (computed from z-scored X) has shape n_samples x n_comps =',
           X_pca.shape[0], 'x', X_pca.shape[1])
    return adata if copy else None

# name for backwards compat
weinreb16 = recipe_weinreb16


def recipe_zheng17(adata, n_top_genes=1000, copy=False):
    if copy: adata = adata.copy()
    pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
    pp.normalize_per_cell(adata)          # normalize with total UMI count per cell
    filter = pp.filter_genes_dispersion(adata.X,
                                        flavor='cell_ranger',
                                        n_top_genes=n_top_genes,
                                        log=False)
    pl.filter_genes_dispersion(filter, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter.gene_subset]
    adata.inplace_subset_var(filter.gene_subset)  # filter genes
    pp.normalize_per_cell(adata)  # need to redo normalization after filtering
    pp.log1p(adata)  # log transform: X = log(X + 1)
    pp.scale(adata, zero_center=True)
    return adata if copy else None
