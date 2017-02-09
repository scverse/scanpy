# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocessing functions that take adata as argument
-> subsampling
-> whole preprocessing workflows from the literature
"""

from .. import settings as sett
from .simple import *

def subsample(adata, subsample, seed=0):
    """ 
    Subsample.

    Parameters
    ----------
    adata : data dictionary
    subsample : int
        Inverse fraction to sample to.
    seed : int
        Root to change subsampling.
            
    Returns
    -------
    adata : dict containing modified entries
        'row_names', 'expindices', 'explabels', 'expcolors'
    """
    from .. import utils
    _, smp_indices = utils.subsample(adata.X,subsample,seed)
    adata = adata[smp_indices]
    for key in ['Xpca']:
        if key in adata:
            adata[key] = adata[key][smp_indices]
    for k in adata.smp_keys():
        adata[k + '_masks'] = adata[k + '_masks'][:, smp_indices]
    adata['subsample'] = True
    return adata

def weinreb16(adata):
    """
    Normalization and filtering as of Weinreb et al. (2016).

    Parameters
    ----------
    adata : dict
        Data dictionary.

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """
    sett.m(0, 'preprocess: weinreb16')
    meanFilter = 0.01
    cvFilter = 2
    nr_pcs = 50

    ddata = adata.to_dict()
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
    ddata['var_names'] = ddata['var_names'][gene_filter]
    sett.m(0, 'Xpca has shape', 
           ddata['Xpca'].shape[0], 'x', ddata['Xpca'].shape[1])
    from ..ann_data import AnnData
    adata = AnnData(ddata)
#     print(adata.X)
#     #return adata

#     # TODO: this should all be organized more nicely
#     # ones retrieve, do everythin on X, then feed back into adata
#     X = adata.X
#     # row normalize
#     X = row_norm(X, max_fraction=0.05, mult_with_mean=True)
#     # filter out genes with mean expression < 0.1 and coefficient of variance <
#     # cvFilter
#     X, gene_filter = filter_genes_cv(X, meanFilter, cvFilter)
#     # compute zscore of filtered matrix
#     Xz = zscore(X)
#     # PCA
#     Xpca = pca(Xz, nr_comps=nr_pcs)
#     # update adata
#     adata.X = X
#     adata = adata.var_names[gene_filter] # filter genes
#     adata['Xpca'] = Xpca
#     sett.m(0, 'Xpca has shape', 
#            adata['Xpca'].shape[0], 'x', adata['Xpca'].shape[1])
#     print(adata.X)

    return adata
