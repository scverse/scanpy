# Author: F. Alex Wolf (http://falexwolf.de)
"""
Preprocessing functions that take X as an argument
-> "simple" preprocessing functions
"""

from ..classes.ann_data import AnnData
import numpy as np

def filter_cells(X, min_reads):
    """
    Filter out cells with total UMI count < min_reads.

    Paramaters
    ----------
    X : np.ndarray
        Data matrix. Rows correspond to cells and columns to genes.
    min_reads : int
        Minimum number of reads required for a cell to survive filtering.

    Returns
    -------
    X : np.ndarray
        Filtered data matrix.
    cell_filter : np.ndarray
        Boolean mask that reports filtering. True means that the cell is
        kept. False means the cell is removed.
    """
    total_counts = np.sum(X, axis=1)
    cell_filter = total_counts >= min_reads
    return X[cell_filter], cell_filter

def gene_filter_cv(X, Ecutoff, cvFilter):
    """
    Filter genes by coefficient of variance and mean.
    """
    mean_filter = np.mean(X, axis=0) > Ecutoff
    var_filter = np.std(X, axis=0) / (np.mean(X,axis=0) + .0001) > cvFilter
    gene_filter = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_filter

def gene_filter_fano(X, Ecutoff, Vcutoff):
    """
    Filter genes by fano factor and mean.
    """
    mean_filter = np.mean(X,axis=0) > Ecutoff
    var_filter = np.var(X,axis=0) / (np.mean(X,axis=0)+.0001) > Vcutoff
    gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
    return gene_filter

def log(X):
    """ 
    Apply logarithm to count data.

    Shifted by one to map 0 to 0.
    """
    if isinstance(X, AnnData):
        return AnnData(log(X.X), X.smp, X.var, **X.add)
    X = np.log(X + 1)
    return X

def pca(X, n_comps=50, zero_center=True, svd_solver='randomized', random_state=0):
    """
    Return PCA representation of data.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    n_comps : int, optional (default: 50)
        Number of PCs to compute.
    zero_center : bool, optional (default: True)
        If True, compute standard PCA from Covariance matrix. If False, omit
        zero-centering variables, which allows to handle sparse input efficiently.
        For sparse intput, automatically defaults to False.
    svd_solver : str, optional (default: 'randomized')
        SVD solver to use. Either “arpack” for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or “randomized” for the randomized algorithm
        due to Halko (2009).

    Returns
    -------
    Y : np.ndarray
        Data projected on n_comps PCs.
    """
    from .. import settings as sett
    if X.shape[1] < n_comps:
        n_comps = X.shape[1]-1
        sett.m(0, 'reducing number of computed PCs to', 
               n_comps, 'as dim of data is only', X.shape[1])
    try:
        from sklearn.decomposition import PCA, TruncatedSVD
        sett.mt(0, 'compute PCA with n_comps =', n_comps)
        from scipy.sparse import issparse
        if zero_center and not issparse(X):
            Y = PCA(n_components=n_comps, svd_solver=svd_solver).fit_transform(X)
        else:
            sett.m(0, '... without zero-centering')
            Y = TruncatedSVD(n_components=n_comps).fit_transform(X)
        sett.mt(0, 'finished')
        sett.m(1, '--> to speed this up, set option exact=False')
    except ImportError:
        Y = _pca_fallback(X, n_comps=n_comps)
        sett.mt(0, 'preprocess: computed PCA using fallback code\n',
                '--> can be sped up by installing package scikit-learn\n',
                '    or by setting the option exact=False')
    return Y

def row_norm(X, max_fraction=1, mult_with_mean=False):
    """
    Normalize so that every cell has the same total read count. 

    Used, for example, by Haghverdi et al. (2016) or Weinreb et al. (2016).

    Using Euclidian distance after this normalisation will yield the same result
    as using cosine distance.
        eucl_dist(Y1,Y2) = (Y1-Y2)^2 
        = Y1^2 +Y2^2 - 2Y1*Y2 = 1 + 1 - 2 Y1*Y2 = 2*(1-(Y1*Y2))
        = 2*(1-(X1*X2)/(|X1|*|X2|)) = 2*cosine_dist(X1,X2)
    
    Parameters
    ----------
    X : np.ndarray
        Expression matrix. Rows correspond to cells and columns to genes.
    max_fraction : float, optional
        Only use genes that make up less than max_fraction of the total
        reads in every cell.
    mult_with_mean: bool, optional
        Multiply the result with the mean of total counts.
        
    Returns
    -------
    Xnormalized : np.ndarray
        Normalized version of the original expression matrix. 
    """
    from .. import settings as sett
    if max_fraction < 0 or max_fraction > 1:
        raise ValueError('choose max_fraction between 0 and 1')
    total_counts = np.sum(X,axis=1)
    if max_fraction == 1:
        X_norm = X / total_counts[:,np.newaxis]
        return X_norm
    # restrict computation of counts to genes that make up less than
    # constrain_theshold of the total reads
    tc_tiled = np.tile(total_counts[:, np.newaxis],(1, X.shape[1]))
    included = np.all(X <= tc_tiled * max_fraction, axis=0)    
    tc_include = np.sum(X[:,included], axis=1)
    tc_tiled = np.tile(tc_include[:, np.newaxis],(1, X.shape[1])) + 1e-6
    X_norm =  X / tc_tiled 
    if mult_with_mean:
        X_norm *= np.mean(total_counts)
    return X_norm

def subsample(adata, subsample, seed=0):
    """ 
    Subsample.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    subsample : int
        Subsample to a fraction of 1/subsample of the data.
    seed : int
        Root to change subsampling.
            
    Returns
    -------
    adata : dict containing modified entries
        'row_names', 'expindices', 'explabels', 'expcolors'
    """
    from .. import utils
    _, smp_indices = utils.subsample(adata.X, subsample, seed)
    adata = adata[smp_indices, ]
    for key in ['X_pca']:
        if key in adata:
            adata[key] = adata[key][smp_indices]
    for k in adata.smp_keys():
        if k + '_masks' in adata:
            adata[k + '_masks'] = adata[k + '_masks'][:, smp_indices]
    adata['subsample'] = True
    return adata

def zscore(X):
    """
    Z-score standardize each column of X.
    
    Parameters
    ----------
    X : np.ndarray
        Data matrix. Rows correspond to cells and columns to genes.
        
    Returns
    -------
    XZ : np.ndarray 
        Z-score standardized version of the data matrix.
    """
    means = np.tile(np.mean(X,axis=0)[None,:],(X.shape[0],1))
    stds = np.tile(np.std(X,axis=0)[None,:],(X.shape[0],1))
    return (X - means) / (stds + .0001)

#--------------------------------------------------------------------------------
# Helper Functions
#--------------------------------------------------------------------------------

def _pca_fallback(data, n_comps=2):
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    C = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since C is symmetric, 
    # the performance gain is substantial
    # evals, evecs = np.linalg.eigh(C)
    evals, evecs = sp.sparse.linalg.eigsh(C, k=n_comps)
    # sort eigenvalues in decreasing order
    idcs = np.argsort(evals)[::-1]
    evecs = evecs[:, idcs]
    evals = evals[idcs]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or n_comps)
    evecs = evecs[:, :n_comps]
    # project data points on eigenvectors
    return np.dot(evecs.T, data.T).T

