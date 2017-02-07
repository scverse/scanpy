# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocessing functions that take X as an argument
-> "simple" preprocessing functions
"""

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

def filter_genes_cv(X, Ecutoff, cvFilter):
    mean_filter = np.mean(X,axis=0)> Ecutoff
    var_filter = np.std(X,axis=0) / (np.mean(X,axis=0)+.0001) > cvFilter
    gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
    return X[:,gene_filter], gene_filter

def filter_genes_fano(X, Ecutoff, Vcutoff):
    mean_filter = np.mean(X,axis=0) > Ecutoff
    var_filter = np.var(X,axis=0) / (np.mean(X,axis=0)+.0001) > Vcutoff
    gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
    return X[:,gene_filter], gene_filter

def log(X):
    """ 
    Apply logarithm to count data.

    Shifted by one to avoid negative infinities.
    """
    X = np.log(X + 1)
    return X

def pca(X, nr_comps=50, exact=True):
    """
    Return PCA representation of data.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    nr_comps : int
        Number of PCs to compute.

    Returns
    -------
    Y : np.ndarray
        Data projected on nr_comps PCs.
    """
    from .. import settings as sett
    if X.shape[1] < nr_comps:
        nr_comps = X.shape[1]-1
        sett.m(0, 'reducing number of computed PCs to', 
               nr_comps, 'as dim of data is only', X.shape[1])
    # deal with multiple PCA implementations
    try:
        from sklearn.decomposition import PCA
        if exact:
            # run deterministic PCA
            svd_solver = 'arpack'
        else:
            # run randomized, more efficient version
            svd_solver = 'randomized'
        p = PCA(n_components=nr_comps, svd_solver=svd_solver)
        Y = p.fit_transform(X)
        sett.mt(0, 'computed PCA with nr_comps =', nr_comps)
        sett.m(1, '--> to speed this up, set option exact=False')
    except ImportError:
        Y = _pca_fallback(X, nr_comps=nr_comps, exact=exact)
        sett.mt(0,'preprocess: computed PCA using fallback code\n',
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
    tc_tiled = np.tile(total_counts[:,np.newaxis],(1,X.shape[1]))
    included = np.all(X <= tc_tiled * max_fraction, axis=0)    
    tc_include = np.sum(X[:,included],axis=1)
    tc_tiled = np.tile(tc_include[:,np.newaxis],(1,X.shape[1])) + 1e-6
    X_norm =  X / tc_tiled 
    if mult_with_mean:
        X_norm *= np.mean(total_counts)
    return X_norm

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

def _pca_fallback(data, nr_comps=2, exact=False):
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    C = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since C is symmetric, 
    # the performance gain is substantial
    if exact:
        evals, evecs = np.linalg.eigh(C)
    else:
        evals, evecs = sp.sparse.linalg.eigsh(C, k=nr_comps)
    # sort eigenvalues in decreasing order
    idcs = np.argsort(evals)[::-1]
    evecs = evecs[:, idcs]
    evals = evals[idcs]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or nr_comps)
    evecs = evecs[:, :nr_comps]
    # project data points on eigenvectors
    return np.dot(evecs.T, data.T).T

