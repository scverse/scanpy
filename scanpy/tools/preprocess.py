# coding: utf-8
"""
Preprocess Data
===============

From package Scanpy (https://github.com/falexwolf/scanpy).
Written in Python 3 (compatible with 2).
Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).

Normalization and filtering functions that can be called directly or via
>>> preprocess(ddata,preprocess_key)
Here, ddata is a data dictionary and preprocess_key is a string that
identifies the preprocessing function.

Beta Version, still lacks many important cases.
"""

import numpy as np
import scipy as sp
from .. import settings as sett
from .. import utils

def preprocess(ddata,preprocess_key,*args,**kwargs):
    """
    Preprocess data with a set of available functions in this module.

    A function is selected based on the 'preprocess key' and
    is passed the optional arguments args and kwargs.

    Parameters
    ----------
    ddata : dict
       Data dictionary containing at least a data matrix 'X'.
    preprocess_key : str
       String that identifies the normalization function in this module.

    Returns
    -------
    ddata : dict
       Data dictionary that stores the preprocessed data matrix.
    """

    # dictionary of all globally defined functions to preprocess
    # data for examples etc.
    preprocess_functions = globals()
    if preprocess_key in preprocess_functions:
        return preprocess_functions[preprocess_key](ddata,*args,**kwargs)
    else:
        raise ValueError('Do not know preprocess function' + preprocess_key
                         + 'try one of: \n' + str(preprocess_functions))

# ------------------------------------------------------------------------------
# Simple preprocessing functions
# ------------------------------------------------------------------------------

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

def log(ddata):
    """ 
    Apply logarithm to pseudocounts.

    Used, for example, by Haghverdi et al. (2016), doi:10.1038/nmeth.3971.
    """
    # apply logarithm to plain count data, shifted by one to avoid 
    # negative infinities, as done in Haghverdi et al. (2016)
    ddata['X'] = np.log(ddata['X']+1)

    # if root cell is defined as expression vector
    if 'xiroot' in ddata and type(ddata['xiroot']) == np.ndarray:
        ddata['xiroot'] = np.log(ddata['xiroot'] + 1)

    return ddata

def pca(X, n_components=2, exact=True):
    """
    Return PCA representation of data.

    Mind that the sklearn implementation is not completely deterministic,
    even if the 'arpack' solver is used. # Alex now thinks it's determinisitic.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    n_components : int
        Number of PCs to compute.

    Returns
    -------
    Y : np.ndarray
        Data projected on n_components PCs.
    """
    # deal with multiple PCA implementations
    try:
        from sklearn.decomposition import PCA
        if exact:
            # run deterministic PCA
            svd_solver = 'arpack'
        else:
            # run randomized, more efficient version
            svd_solver = 'randomized'
        p = PCA(n_components=n_components, svd_solver=svd_solver)
        sett.m(0,'compute PCA using sklearn')
        # sett.m(0,'--> mind that this is not 100% deterministic')
        # Alex now thinks it's determinisitic.
        sett.m(1,'--> to speed this up, set option exact=False')
        Y = p.fit_transform(X)
    except ImportError:
        sett.m(0,'compute PCA using fallback code\n',
                 '--> can be sped up by installing package scikit-learn\n',
                 '    or by setting the option exact=False')
        Y = _pca_fallback(X, n_components=n_components, exact=exact)
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
    sett.m(0,'preprocess: normalizing rows of data matrix')
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

def subsample(ddata, subsample, seed=0):
    """ 
    Subsample.

    Parameters
    ----------
    ddata : data dictionary
    subsample : int
        Inverse fraction to sample to.
    seed : int
        Root to change subsampling.
            
    Returns
    -------
    ddata : dict containing modified entries
        'rownames', 'expindices', 'explabels', 'expcolors'
    """
    X, row_indices = utils.subsample(ddata['X'],subsample,seed)
    for key in ['rownames']:
        if key in ddata and len(ddata[key]) == ddata['X'].shape[0]:
            ddata[key] = ddata[key][row_indices]
    if 'groupmasks' in ddata:
        ddata['groupmasks'] = ddata['groupmasks'][:,row_indices]
    ddata['X'] = X
    ddata['subsample'] = True
    return ddata

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
# Whole preprocessing workflows from the literature
#--------------------------------------------------------------------------------

def weinreb16(ddata):
    """
    Normalization and filtering as of Weinreb et al. (2016).

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """

    sett.m(0, 'preprocess: weinreb16')

    meanFilter = 0.01
    cvFilter = 2
    numPCs = 50

    X = ddata['X']

    # filter out cells with fewer than 1000 UMIs
    # X, cell_filter = filter_cells(X,1000)
    # ddata['rownames'] = ddata['rownames'][cell_filter]

    # row normalize
    X = row_norm(X, max_fraction=0.05, mult_with_mean=True)

    # filter out genes with mean expression < 0.1 and coefficient of variance <
    # cvFilter
    _, gene_filter = filter_genes_cv(X, meanFilter, cvFilter)

    # compute zscore of filtered matrix
    X_z = zscore(X[:, gene_filter])

    # PCA
    X_pca = pca(X_z, n_components=numPCs)

    ddata['X'] = X_pca
    sett.m(0, 'after PCA, X has shape', 
           ddata['X'].shape[0], 'x', ddata['X'].shape[1])

    return ddata

#--------------------------------------------------------------------------------
# Helper Functions
#--------------------------------------------------------------------------------

def _pca_fallback(data, n_components=2, exact=False):
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
        evals, evecs = sp.sparse.linalg.eigsh(C, k=n_components)
    # sort eigenvalues in decreasing order
    idcs = np.argsort(evals)[::-1]
    evecs = evecs[:, idcs]
    evals = evals[idcs]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or n_components)
    evecs = evecs[:, :n_components]
    # project data points on eigenvectors
    return np.dot(evecs.T, data.T).T


