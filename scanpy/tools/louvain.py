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
            flavor='igraph',
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
    # if 'distance' not in adata.add or recompute_graph:
    #     graph = data_structs.DataGraph(adata,
    #                                    k=n_neighbors,
    #                                    n_pcs=n_pcs,
    #                                    n_jobs=n_jobs)
    #     graph.compute_distance_matrix()
    #     adata.add['distance'] = graph.Dsq
    if 'Ktilde' not in adata.add or recompute_graph:
        graph = data_structs.DataGraph(adata,
                                       k=n_neighbors,
                                       n_pcs=n_pcs,
                                       n_jobs=n_jobs)
        graph.compute_transition_matrix()
        adata.add['Ktilde'] = graph.Ktilde
    adjacency = adata.add['Ktilde']
    if flavor in {'vtraag', 'igraph'}:
        sources, targets = adjacency.nonzero()
        weights = adjacency[sources, targets]
        weights = np.array(weights)[0]  # need to convert sparse matrix into a form appropriate for igraph
        if flavor == 'igraph' and resolution is not None:
            logg.warn('`resolution` parameter has no effect for flavor "igraph"')
        if directed and flavor == 'igraph':
            directed = False
        if not directed: logg.info('    using the undirected graph')
        g = ig.Graph(list(zip(sources, targets)), directed=directed, edge_attrs={'weight': weights})
        if flavor == 'vtraag':
            if resolution is None: resolution = 1
            # part = louvain.find_partition(g, method='RBConfiguration',
            #                               resolution_parameter=resolution)
            part = louvain.find_partition(g, louvain.ModularityVertexPartition)
                                          # resolution_parameter=resolution)
        elif flavor == 'igraph':
            part = g.community_multilevel()
        groups = np.array(part.membership, dtype='U')
    elif flavor == 'taynaud':
        import networkx as nx
        import community
        g = nx.Graph(adata.add['distance'])
        partition = community.best_partition(g)
        groups = np.zeros(len(partition), dtype=int)
        for k, v in partition.items(): groups[k] = v
        groups = groups.astype('U')
    else:
        raise ValueError('flavor needs to be "vtraag" or "igraph"')
    adata.smp['louvain_groups'] = groups
    from natsort import natsorted
    adata.add['louvain_groups_names'] = np.array(natsorted(np.unique(groups)))
    logg.m('    finished', t=True, end=' ')
    logg.m('and found', len(adata.add['louvain_groups_names']), 'clusters, added\n'
           '    "louvain_groups", the cluster labels (adata.smp)\n'
           '    "louvain_groups_names", the unique cluster labels (adata.add)')
    return adata if copy else None
