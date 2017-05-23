# Author: F. Alex Wolf (http://falexwolf.de)
"""Cluster using DBSCAN

Using the scikit-learn implementation.
"""

import numpy as np
from .. import settings as sett
from .. import logging as logg

def dbscan(adata, basis='tsne', n_comps=2, eps=None, min_samples=None, n_jobs=None, copy=False):
    """Cluster cells using DBSCAN

    This wraps sklearn.cluster.DBSCAN and shares most of the parameters.

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
    n_jobs : int (default: None)
        Number of threads to use. Defaults to sett.n_jobs.
    copy : bool (default: False)

    References
    ----------
    Ester et al. (1996), "A Density-Based Algorithm for Discovering Clusters in
    Large Spatial Databases with Noise".
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231.

    Pedregosa et al. (2011) ...
    """
    logg.m('starting DBSCAN', r=True)
    adata = adata.copy() if copy else adata
    if basis not in {'tsne', 'pca'}:
        raise ValueError('`basis` needs to be "tsne" or "pca"')
    if 'X_tsne' in adata.smp and basis == 'tsne':
        X = adata.smp['X_tsne'][:, :n_comps]
    elif 'X_pca' in adata.smp and basis == 'pca':
        X = adata.smp['X_pca'][:, :n_comps]
    else:
        raise ValueError('Run {} first.'.format(basis))
    n_jobs = sett.n_jobs if n_jobs is None else n_jobs
    range_1 = np.max(X[:, 0]) - np.min(X[:, 0])
    range_2 = np.max(X[:, 1]) - np.min(X[:, 1])

    if eps is None:
        if n_comps == 2:
            avg_area_per_point = (range_1 * range_2 / X.shape[0])
            logg.m('... the "drawing range" is', range_1, 'x', range_2,
                   'with the average area per point', avg_area_per_point)
            eps = 1.7 * np.sqrt(avg_area_per_point)
        else:
            eps = 5
    if min_samples is None: min_samples = 30
    logg.m('... using eps =', eps, end=', ')
    logg.m('min_samples =', min_samples, end=', ')
    logg.m('basis =', basis, end=', ')
    logg.m('n_comps =', basis, end=', ')
    logg.m('n_jobs =', n_jobs)  #, end=', ')
    logg.m('increase `min_samples` if you find too many clusters', v='hint')
    logg.m('reduce eps if "everything is connected"', v='hint')
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=min_samples, n_jobs=n_jobs)
    nn.fit(X)
    D = nn.radius_neighbors_graph(radius=eps, mode='distance')
    db = DBSCAN(eps=eps, min_samples=min_samples,
                n_jobs=n_jobs, metric='precomputed').fit(D)
    labels = db.labels_
    dont_know = labels == -1
    labels = labels.astype(str)
    labels[dont_know] = '?'
    # loop_over_labels = (label for label in np.unique(labels) if label >= 0)
    adata.smp['dbscan_groups'] = labels
    from natsort import natsorted
    adata.add['dbscan_groups_names'] = np.array(natsorted(np.unique(labels)))[:-1]
    logg.m('    finished', t=True, end=' ')
    logg.m('and found', len(np.unique(labels))-1, 'clusters, added\n'
           '    "dbscan_groups", the cluster labels (adata.smp)\n'
           '    "dbscan_groups_names", the unique cluster labels (adata.add)')
    return adata if copy else None
