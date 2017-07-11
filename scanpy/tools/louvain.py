# Author: F. Alex Wolf (http://falexwolf.de)
"""Cluster using Louvain community detection algorithm

Uses the pip packages "louvain" and igraph.
"""

import numpy as np
from .. import settings as sett
from .. import logging as logg
from .. import data_structs

def louvain(adata,
            n_neighbors=30,
            n_pcs=50,
            resolution=None,
            flavor='vtraag',
            directed=True,
            recompute_graph=False,
            n_jobs=None,
            copy=False):
    """Cluster cells using Louvain Community detection.

    This has been suggested for single-cell transcriptomics by
    Levine et al., Cell 162, 184-197 (2015).

    Parameters
    ----------
    n_neighbors : int, optional (default: 30)
        Number of neighbors to use for construction of knn graph.
    resolution : float or None, optional
        For the default flavor, you provide a resolution, which defaults to 1.0.
    copy : bool (default: False)

    References
    ----------
    Algorithm: Blondel et al., J. Stat. Mech., P10008 (2008)
    Package: Csardi et al., InterJournal Complex Systems, 1695 (2006)
    """
    logg.m('run Louvain clustering', r=True)
    adata = adata.copy() if copy else adata
    import igraph as ig
    import louvain
    if 'distance' not in adata.add or recompute_graph:
        graph = data_structs.DataGraph(adata,
                                       k=n_neighbors,
                                       n_pcs=n_pcs,
                                       n_jobs=n_jobs)
        graph.compute_distance_matrix()
        adata.add['distance'] = graph.Dsq
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
