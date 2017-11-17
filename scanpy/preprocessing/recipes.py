# Author: Alex Wolf (http://falexwolf.de)
"""Preprocessing recipes from the literature
"""

from . import simple as pp


def recipe_weinreb16(adata, mean_threshold=0.01, cv_threshold=2,
                     n_pcs=50, svd_solver='randomized', random_state=0, copy=False):
    """Normalization and filtering as of [Weinreb17]_.

    Expects logarithmized data.

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
    """
    from scipy.sparse import issparse
    if issparse(adata.X):
        raise ValueError('recipe_weinreb16() does not support sparse matrices.')
    if copy: adata = adata.copy()
    adata.X = pp.normalize_per_cell_weinreb16_deprecated(adata.X,
                                                         max_fraction=0.05,
                                                         mult_with_mean=True)
    # filter out genes with mean expression < 0.1 and coefficient of variance <
    # cv_threshold
    gene_subset = pp.filter_genes_cv_deprecated(adata.X, mean_threshold, cv_threshold)
    # three alternative ways of slicing
    adata._inplace_subset_var(gene_subset)  # this modifies the object itself
    # adata = adata[:, gene_subset]  # this does a copy
    # adata[:, ~gene_subset] = None  # this doesn't work yet
    # compute zscore of filtered matrix and compute PCA
    X_pca = pp.pca(pp.zscore_deprecated(adata.X),
                   n_comps=n_pcs, svd_solver=svd_solver, random_state=random_state)
    # update adata
    adata.smpm['X_pca'] = X_pca
    return adata if copy else None


def recipe_zheng17(adata, n_top_genes=1000, zero_center=True, plot=False, copy=False):
    """Normalization and filtering as of [Zheng17]_.

    Expects non-logarithmized data.

    This reproduces the preprocessing of the reference below, at the time, the
    Cell Ranger R Kit preprocessing of 10X Genomics.

    Parameters
    ----------
    n_top_genes : int, optional (default: 1000)
        Number of genes to keep.
    zero_center : bool, optional (default: True)
        Zero center the data matrix. Only switch this to False if you have
        serious memory problems.
    plot : bool, optional (default: True)
        Show a plot of the gene dispersion vs. mean relation.
    copy : bool, optional (default: False)
        Return a copy of adata instead of updating the passed object.
    """
    if copy: adata = adata.copy()
    pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
    pp.normalize_per_cell(adata,  # normalize with total UMI count per cell
                          key_n_counts='n_counts_all')
    filter_result = pp.filter_genes_dispersion(adata.X,
                                               flavor='cell_ranger',
                                               n_top_genes=n_top_genes,
                                               log=False)
    if plot:
        from .. import plotting as pl  # should not import at the top of the file
        pl.filter_genes_dispersion(filter_result, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter.gene_subset]
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    pp.normalize_per_cell(adata)  # need to redo normalization after filtering
    pp.log1p(adata)  # log transform: X = log(X + 1)
    pp.scale(adata, zero_center=zero_center)
    return adata if copy else None
