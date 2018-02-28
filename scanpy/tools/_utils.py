from .. import logging as logg
from .pca import pca

def preprocess_with_pca(adata, n_pcs=None, random_state=0):
    """
    Parameters
    ----------
    n_pcs : `int` or `None`, optional (default: `None`)
        If `n_pcs=0`, do not preprocess with PCA.
        If `None` and there is a PCA version of the data, use this.
        If an integer, compute the PCA.
    """
    if n_pcs == 0:
        logg.info('    using data matrix X directly (no PCA)')
        return adata.X
    elif n_pcs is None and 'X_pca' in adata.obsm_keys():
        logg.info('    using \'X_pca\' with n_pcs = {} for tSNE'
                  .format(adata.obsm['X_pca'].shape[1]))
        return adata.obsm['X_pca']
    elif ('X_pca' in adata.obsm_keys()
          and adata.obsm['X_pca'].shape[1] >= n_pcs):
        logg.info('    using \'X_pca\' with n_pcs = {} for tSNE'
                  .format(n_pcs))
        return adata.obsm['X_pca'][:, :n_pcs]
    else:
        if adata.X.shape[1] > n_pcs:
            logg.info('    computing \'X_pca\' with n_pcs = {}'.format(n_pcs))
            logg.hint('avoid this by setting n_pcs = 0')
            X = pca(adata.X, n_comps=n_pcs, random_state=random_state)
            adata.obsm['X_pca'] = X
            return X
        else:
            logg.warn('Data matrix has lower dimension than passed number of PCs. '
                      'Using data matrix X directly (no PCA).')
            return adata.X
