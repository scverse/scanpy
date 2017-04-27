# Author: F. Alex Wolf (http://falexwolf.de)
"""
Preprocessing recipes from the literature
"""

from .. import settings as sett
from .simple import normalize_per_cell_weinreb16, filter_genes_cv, pca, zscore

def weinreb16(adata, mean_threshold=0.01, cv_threshold=2,
              n_pcs=50, svd_solver='randomized', random_state=0):
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

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016).
    """
    sett.m(0, 'preprocess: weinreb16, X has shape n_samples x n_variables =',
           adata.X.shape[0], 'x', adata.X.shape[1])
    # row normalize
    adata.X = cell_norm_weinreb16(adata.X, max_fraction=0.05, mult_with_mean=True)
    # filter out genes with mean expression < 0.1 and
    # coefficient of variance < cv_threshold
    gene_filter = gene_filter_cv(adata.X, mean_threshold, cv_threshold)
    adata = adata[:, gene_filter]
    # compute zscore of filtered matrix and compute PCA
    X_pca = pca(zscore(adata.X),
                n_comps=n_pcs, svd_solver=svd_solver, random_state=random_state)
    # update adata
    adata['X_pca'] = X_pca
    sett.m(0, 'X_pca (computed from z-scored X) has shape n_samples x n_comps =',
           adata['X_pca'].shape[0], 'x', adata['X_pca'].shape[1])
    return adata
