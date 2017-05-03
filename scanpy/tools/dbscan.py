# Author: F. Alex Wolf (http://falexwolf.de)
"""Cluster using DBSCAN

Using the scikit-learn implementation.
"""

import numpy as np
from .. import settings as sett


def dbscan(adata, basis='tsne', n_comps=10, eps=None, min_samples=None, copy=False):
    """Cluster using DBSCAN

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
    sett.mt(0, 'starting DBSCAN', start=True)
    adata = adata.copy() if copy else adata
    found = False
    if 'X_tsne' in adata.smp and basis == 'tsne':
        X = adata.smp['X_tsne']
        found = True
    if 'X_pca' in adata.smp and basis == 'pca':
        X = adata.smp['X_pca']
        found = True
    if basis not in {'tsne', 'pca'}:
        raise ValueError('`basis` needs to be "tsne" or "pca"')
    if not found: raise ValueError('Run {} first.'.format(basis))
    if eps is None:
        if n_comps == 2:
            # average area per point
            avg_area_per_point = ((np.max(X[:, 0]) - np.min(X[:, 0]))
                                  * (np.max(X[:, 1]) - np.min(X[:, 1])) / X.shape[0])
            eps = 3*np.sqrt(avg_area_per_point)  # reduce a bit further
        else:
            eps = 5
        sett.m(0, '... using eps', eps)
    if min_samples is None:
        min_samples = int(X.shape[0] / 120)
        sett.m(0, '... using min_samples', min_samples)
    from sklearn.cluster import DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    labels = db.labels_
    sett.m(0, 'found', len(np.unique(labels))-1, 'clusters')
    dont_know = labels == -1
    labels = labels.astype(str)
    labels[dont_know] = '?'
    # loop_over_labels = (label for label in np.unique(labels) if label >= 0)
    adata.smp['dbscan_groups'] = labels
    adata.add['dbscan_groups_names'] = np.unique(labels)
    sett.mt(0, 'finished, added\n'
            '    "dbscan_groups", the cluster labels (adata.smp)\n'
            '    "dbscan_groups_names", the unique cluster labels (adata.smp_names)')
    return adata if copy else None
