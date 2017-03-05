# coding: utf-8
# Author: F. Alex Wolf (http://falexwolf.de)
"""
Cluster using DBSCAN

Using the scikit-learn implementation.
"""

import numpy as np
from .. import settings as sett
from .. import utils

def dbscan(adata, eps=None, min_samples=None):
    """
    Cluster using DBSCAN.

    Parameters
    ----------
    eps : float or None, optional
        The maximum distance between samples for being considered as in the same
        neighborhood. Clusters are "grown" from samples that have more than
        min_samples points in their neighborhood. Increasing eps therefore 
        allows clusters to spread over wider regions.
    min_samples : int or None, optional
        The number of samples (or total weight) in a neighborhood for a point
        to be considered as a core point. This includes the point itself.

    References
    ----------
    Ester et al. (1996), "A Density-Based Algorithm for Discovering Clusters in
    Large Spatial Databases with Noise".
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231.

    Pedregosa et al. (2011) ...
    """
    if 'X_tsne' in adata:
        X = adata['X_tsne']
    elif 'X_pca' in adata:
        X = adata['X_pca']
    else:
        raise ValueError('perform tSNE first')
    if eps is None:
        # average area per point
        avg_area_per_point = ((np.max(X[:, 0]) - np.min(X[:, 0]))
                              * (np.max(X[:, 1]) - np.min(X[:, 1])) / X.shape[0])
        eps = 3*np.sqrt(avg_area_per_point)
        # reduce a bit further
        sett.m(0, '... using eps', eps)
    if min_samples is None:
        min_samples = int(X.shape[0] / 120)
        sett.m(0, '... using min_samples', min_samples)
    from sklearn.cluster import DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    labels = db.labels_
    mask = labels == -1
    sett.m(0, 'found', len(np.unique(labels))-1, 'clusters')
    labels = labels.astype(str)
    labels[mask] = '?'
    # loop_over_labels = (label for label in np.unique(labels) if label >= 0)
    adata.smp['dbscan_groups'] = labels
    adata['dbscan_groups_names'] = np.unique(labels)
    return adata

def plot_dbscan(adata,
         basis='tsne',
         smp=None,
         names=None,
         comps=None,
         cont=None,
         layout='2d',
         legendloc='right margin',
         cmap=None,
         pal=None,
         right_margin=None,
         size=3,
         titles=None):
    """
    Plot results of DPT analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    smps = ['dbscan_groups']
    if smp is not None:
        smps += smp.split(',')
    from .. import plotting as plott
    smps = plott.scatter(adata,
                basis=basis,
                smp=smps,
                names=names,
                comps=comps,
                cont=cont,
                layout=layout,
                legendloc=legendloc,
                cmap=cmap,
                pal=pal,
                right_margin=right_margin,
                size=size,
                titles=titles)

    writekey = sett.basekey + '_dbscan_'+ basis
    writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
    plott.savefig(writekey)

    if not sett.savefigs and sett.autoshow:
        from ..compat.matplotlib import pyplot as pl
        pl.show()
