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
    adata : AnnData
        Annotated data matrix
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
    adata = adata[smp_indices, ]
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
    adata : AnnData
        Annotated data matrix.

    Reference
    ---------
    Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """
    sett.m(0, 'preprocess: weinreb16')
    mean_filter = 0.01
    cv_filter = 2
    nr_pcs = 50

    # row normalize
    adata.X = row_norm(adata.X, max_fraction=0.05, mult_with_mean=True)
    # filter out genes with mean expression < 0.1
    # and coefficient of variance < cv_filter
    gene_filter = gene_filter_cv(adata.X, mean_filter, cv_filter)
    adata = adata[:, gene_filter]  # filter genes

    # compute zscore of filtered matrix and create PCA
    Xz = zscore(adata.X)
    Xpca = pca(Xz, nr_comps=nr_pcs)

    # update adata
    adata['Xpca'] = Xpca
    sett.m(0, 'Xpca has shape', adata['Xpca'].shape[0], 'x', adata['Xpca'].shape[1])

    return adata
