# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocessing functions that take ddata as argument
-> subsampling
-> whole preprocessing workflows from the literature
"""

from .. import settings as sett
from .simple import *

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
        'row_names', 'expindices', 'explabels', 'expcolors'
    """
    from .. import utils
    X, row_indices = utils.subsample(ddata['X'],subsample,seed)
    for key in ['row_names', 'Xpca']:
        if key in ddata and len(ddata[key]) == ddata['X'].shape[0]:
            ddata[key] = ddata[key][row_indices]
    if 'rowcat' in ddata:
        for k in ddata['rowcat']:
            ddata['rowcat'][k] = ddata['rowcat'][k][row_indices]
            ddata[k + '_masks'] = ddata[k + '_masks'][:, row_indices]
    ddata['X'] = X
    ddata['subsample'] = True
    return ddata

def weinreb16(ddata):
    """
    Normalization and filtering as of Weinreb et al. (2016).

    Parameters
    ----------
    ddata : dict
        Data dictionary.

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """
    sett.m(0, 'preprocess: weinreb16')
    meanFilter = 0.01
    cvFilter = 2
    nr_pcs = 50
    X = ddata['X']
    # row normalize
    X = row_norm(X, max_fraction=0.05, mult_with_mean=True)
    # filter out genes with mean expression < 0.1 and coefficient of variance <
    # cvFilter
    X, gene_filter = filter_genes_cv(X, meanFilter, cvFilter)
    # compute zscore of filtered matrix
    Xz = zscore(X)
    # PCA
    Xpca = pca(Xz, nr_comps=nr_pcs)
    # update dictionary
    ddata['X'] = X
    ddata['Xpca'] = Xpca
    ddata['col_names'] = ddata['col_names'][gene_filter]
    sett.m(0, 'Xpca has shape', 
           ddata['Xpca'].shape[0], 'x', ddata['Xpca'].shape[1])
    return ddata

