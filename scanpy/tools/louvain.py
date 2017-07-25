# Author: F. Alex Wolf (http://falexwolf.de)
"""Cluster using Louvain community detection algorithm

Uses the pip packages "louvain" and "igraph".
"""

import numpy as np
from .. import utils
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

    The basic method has been suggested for single-cell transcriptomics by
    Levine et al., Cell 162, 184-197 (2015).

    Wolf et al. (2017) suggested to compute, in addition, an "attachedness"
    matrix.

    Parameters
    ----------
    n_neighbors : int, optional (default: 30)
        Number of neighbors to use for construction of knn graph.
    resolution : float or None, optional
        For the default flavor, you provide a resolution, which defaults to 1.0.
    copy : bool (default: False)

    References
    ----------
    - implementation of Louvain algorithm: Traag, doi:10.5281/zenodo.35117 (2017)
    - Louvain algorithm: Blondel et al., J. Stat. Mech., P10008 (2008)
    - base graph package: Csardi et al., InterJournal Complex Systems, 1695 (2006)
    - basic suggestion for single-cell: Levine et al., Cell 162, 184-197 (2015)
    - combination with "attachedness" matrix: Wolf et al., bioRxiv (2017)
    """
    logg.m('run Louvain clustering', r=True)
    adata = adata.copy() if copy else adata
    # TODO: it makes not much of a difference, but we should keep this here
    # and make it an option at some point
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
        if flavor == 'igraph' and resolution is not None:
            logg.warn('`resolution` parameter has no effect for flavor "igraph"')
        if directed and flavor == 'igraph':
            directed = False
        if not directed: logg.info('    using the undirected graph')
        g = utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if flavor == 'vtraag':
            import louvain
            if resolution is None: resolution = 1
            try:
                louvain.set_rng_seed(0)
                part = louvain.find_partition(g, louvain.RBConfigurationVertexPartition,
                                              resolution_parameter=resolution)
                adata.add['louvain_quality'] = part.quality()
            except AttributeError:
                logg.warn('Did not find louvain package >= 0.6 on your system, '
                          'the result will therefore not be 100% reproducible, but '
                          'is influenced by randomness in the community detection '
                          'algorithm. Still you get very meaningful results!\n'
                          'If you want 100% reproducible results, but 0.6 is not yet '
                          'available via "pip install louvain", '
                          'either get the latest (development) version from '
                          'https://github.com/vtraag/louvain-igraph or use the option '
                          '`flavor=igraph` in sc.tl.louvain(). '
                          'The latter does not a provide `resolution` parameter, though.')
                part = louvain.find_partition(g, method='RBConfiguration',
                                              resolution_parameter=resolution)
        elif flavor == 'igraph':
            part = g.community_multilevel()
        groups = np.array(part.membership, dtype='U')
    elif flavor == 'taynaud':
        # this is deprecated
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
