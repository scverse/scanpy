from ..preprocessing import pca
from typing import Optional, Union
import numpy as np
from anndata import AnnData
from .. import logging as logg
from .._utils import AnyRandom

def multi_batch_pca(adata: AnnData,                     
                    n_comps: int = None,                    
                    weights: Optional[np.ndarray] = None,                     
                    batch_key: str = 'batch',                     
                    svd_solver: str = 'arpack',                    
                    random_state: AnyRandom = 0,                    
                    return_info: bool = False,                    
                    dtype: str = 'float32',                    
                    **kwargs):    
    """\        
    Multi-batch principal component analysis                     
    Parameters    
    ----------    
    adata        
        The (annotated) AnnData object with .X of shape `n_obs` × `n_vars`.    
    
    weights        
        The weight of each batch in computation of aggregated princial componenets        
        of the adata.    
    batch_key        
        The key value for division of the adata.obs into batches.    
    svd_solver        
        SVD solver to use:
            
        `'arpack'` (the default) 
          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)        
        `'randomized'`          
          for the randomized algorithm due to Halko (2009).        
        `'auto'`          
          chooses automatically depending on the size of the problem.        
        `'lobpcg'`          
          An alternative SciPy solver.  
        
        Efficient computation of the principal components of a sparse matrix        
        currently only works with the `'arpack`' or `'lobpcg'` solvers.
        
    random_state        
        Change to use different initial states for the optimization.    

        
    Returns
    -------
    adata : anndata.AnnData

        `.obsm['X_multi_batch_pca']`
             PCA representation of data.
        `.varm['multi_batch_PCs']`
             The principal components containing the loadings.
        `.uns['multi_batch_pca']['variance_ratio']`
             Ratio of explained variance.
        `.uns['multi_batch_pca']['variance']`
             Explained variance, equivalent to the eigenvalues of the
             covariance matrix.
    """     
        
    if not isinstance(data, AnnData):
        raise ValueError('`adata` must be an `AnnData` object.')    
        
    start = logg.info(f'computing multi-batch PCA')      
    
    batches = list(adata.obs[batch_key].unique())        
    
    if weights is None:        
        weights = np.ones(len(batches))            
    elif len(weights) != len(batches):
        raise ValueError('Length of `weights` must be equal to the number of unique values in `adata.obs[batch_key]`.')            

    n_comps = min(adata.n_vars, n_comps)        
    means = np.zeros(adata.n_vars)    
    scales = np.ones(adata.n_obs)
    

    batch_inds = []
    for i, batch in enumerate(batches):        
        batch_inds.appaned(adata.obs[batch_key] == batch)        
        adata_batch = adata[batch_inds[-1]]
        
        means += weights[i] * np.ravel(adata_batch.X.mean(0))        
        scales[batch_inds[-1]] = np.sqrt(adata_batch.n_obs / weights[i])      
        
    means /= sum(weights)         
    
    from sklearn.decomposition import PCA    
    if not issparse(X):

        
        X_proc = (adata.X - means)/scales[:, None]
      
        pca_ = PCA(
            n_components=n_comps, svd_solver=svd_solver, random_state=random_state
        )
        X_pca = pca_.fit_transform(X_proc)
                
        adata.obsm['X_multi_batch_pca'] = X_pca
        adata.varm['multi_batch_PCs'] = pca_.components_.T
        adata.uns['multi_batch_pca']['variance'] = pca_.explained_variance_
        adata.uns['multi_batch_pca']['variance_ratio'] = pca_.explained_variance_ratio_

    
    else:     
        if svd_solver == "randomized":
            # This is for backwards compat. Better behaviour would be to either error or use arpack.
            logg.warning(
                "svd_solver 'randomized' does not work with sparse input. Densifying the array. "
                "This may take a very large amount of memory."
            )

        adata.X = check_array(adata.X, accept_sparse=['csr', 'csc'])
            
        mdot = means.dot    
        mmat = mdot    
        mhdot = means.T.dot    
        mhmat = mhdot
        Xdot = X.dot    
        Xmat = Xdot    
        XHdot = X.T.conj().dot    
        XHmat = XHdot    
        ones = np.ones(X.shape[0])[None, :].dot
    
        def matvec(x):        
            return (Xdot(x) - mdot(x)) / scales[:, None]
        def matmat(x):        
            return (Xmat(x) - mmat(x)) / scales[:, None]
        def rmatvec(x):        
            return (XHdot(x) - mhdot(ones(x)) / scales[:, None]
        def rmatmat(x):        
            return (XHmat(x) - mhmat(ones(x)) / scales[:, None]
    
        XL = LinearOperator(
            matvec=matvec,
            dtype=X.dtype, 
            matmat=matmat, 
            shape=X.shape, 
            rmatvec=rmatvec, 
            rmatmat=rmatmat
            )
    
        u, s, v = svds(XL, solver=solver, k=npcs, v0=random_init)
        u, v = svd_flip(u, v)
        idx = np.argsort(-s)
        v = v[idx, :]

        X_pca = (u * s)[:, idx]
        ev = s[idx] ** 2 / (X.shape[0] - 1)

        total_var = _get_mean_var(X)[1].sum()
        ev_ratio = ev / total_var

        
        adata.obsm['X_multi_batch_pca'] = X_pca
        adata.varm['multi_batch_PCs'] = v.T 
        adata.uns['multi_batch_pca']['variance'] = ev 
        adata.uns['multi_batch_pca']['variance_ratio'] = ev_ratio 
        
    logg.info('finished', time=start)
    logg.debug(
        'and added\n'
        '    \'X_pca\', the PCA coordinates (adata.obs)\n'
        '    \'PC1\', \'PC2\', ..., the loadings (adata.var)\n'
        '    \'pca_variance\', the variance / eigenvalues (adata.uns)\n'
        '    \'pca_variance_ratio\', the variance ratio (adata.uns)'
    )         
    return adata