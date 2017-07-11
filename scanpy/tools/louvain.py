# Author: F. Alex Wolf (http://falexwolf.de)
"""Cluster using Louvain community detection algorithm

Uses the pip packages "louvain" and igraph.
"""

import numpy as np
from .. import settings as sett
from .. import logging as logg

def louvain(adata, resolution=None, flavor='vtraag', directed=True, copy=False):
    """Cluster cells using Louvain Community detection.

    This has been suggested for single-cell transcriptomics by
    Levine et al., Cell 162, 184-197 (2015).

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
    copy : bool (default: False)

    References
    ----------
    Algorithm: Blondel et al., J. Stat. Mech., P10008 (2008)
    Package: Csardi et al., InterJournal Complex Systems, 1695 (2006)
    """
    logg.m('starting Louvain clustering', r=True)
    adata = adata.copy() if copy else adata
    import igraph as ig
    import louvain
    data = adata.add['distance']
    sources, targets = data.nonzero()
    weights = data[sources, targets]
    weights = np.array(weights)[0]  # need to convert sparse matrix into a form appropriate for igraph
    if flavor == 'igraph' and resolution is not None:
        logg.warn('`resolution` parameter has no effect for flavor "igraph"')
    if directed and flavor == 'igraph':
        directed = False
    if not directed: logg.info('    using the undirected graph')
    g = ig.Graph(list(zip(sources, targets)), directed=directed, edge_attrs={'weight': weights})
    if flavor == 'vtraag':
        if resolution is None: resolution = 1
        part = louvain.find_partition(g, method='RBConfiguration', resolution_parameter=resolution)
    elif flavor == 'igraph':
        part = g.community_multilevel()
    else:
        raise ValueError('flavor needs to be "vrtraag" or "igraph"')
    groups = np.array(part.membership, dtype='U')
    adata.smp['louvain_groups'] = groups
    from natsort import natsorted
    adata.add['louvain_groups_names'] = np.array(natsorted(np.unique(groups)))
    logg.m('    finished', t=True, end=' ')
    logg.m('and found', len(part), 'clusters, added\n'
           '    "louvain_groups", the cluster labels (adata.smp)\n'
           '    "louvain_groups_names", the unique cluster labels (adata.add)')
    return adata if copy else None
