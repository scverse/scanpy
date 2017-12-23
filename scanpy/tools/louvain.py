"""Cluster cells using Louvain community detection algorithm.

Uses the pip package "louvain" by V. Traag.
"""

import numpy as np
import pandas as pd
from natsort import natsorted
from .. import utils
from .. import settings
from .. import logging as logg
from ..data_structs.data_graph import add_or_update_graph_in_adata


def louvain(adata,
            n_neighbors=None,
            resolution=None,
            n_pcs=50,
            random_state=0,
            restrict_to=None,
            key_added=None,
            flavor='vtraag',
            directed=True,
            recompute_pca=False,
            recompute_distances=False,
            recompute_graph=False,
            n_dcs=None,
            n_jobs=None,
            copy=False):
    """Cluster cells into subgroups [Blondel08]_ [Levine15]_ [Traag17]_.

    Cluster cells using the Louvain algorithm [Blondel08]_ in the implementation
    of [Traag17]_. The Louvain algorithm has been proposed for single-cell
    analysis by [Levine15]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        The annotated data matrix.
    n_neighbors : `int`, optional (default: 30)
        Number of neighbors to use for construction of knn graph.
    resolution : `float` or `None`, optional (default: 1)
        For the default flavor ('vtraag'), you can provide a resolution (higher
        resolution means finding more and smaller clusters), which defaults to
        1.0.
    n_pcs : int, optional (default: 50)
        Number of PCs to use for computation of data point graph.
    random_state : int, optional (default: 0)
        Change the initialization of the optimization.
    key_added : str, optional (default: `None`)
        Key under which to add the cluster labels.
    restrict_to : tuple, optional (default: None)
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain (obs key, list of categories).
    flavor : {'vtraag', 'igraph'}
        Choose between to packages for computing the clustering. 'vtraag' is
        much more powerful.
    copy : `bool` (default: False)
        Copy adata or modify it inplace.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    louvain_groups : `pd.Series` (``adata.obs``, dtype `category`)
        Array of dim (number of samples) that stores the subgroup id ('0',
        '1', ...) for each cell.
    """
    logg.info('running Louvain clustering', r=True)
    adata = adata.copy() if copy else adata
    add_or_update_graph_in_adata(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        n_dcs=n_dcs,
        recompute_pca=recompute_pca,
        recompute_distances=recompute_distances,
        recompute_graph=recompute_graph,
        n_jobs=n_jobs)
    adjacency = adata.uns['data_graph_norm_weights']
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        if not isinstance(restrict_categories[0], str):
            raise ValueError('You need to use strings to label categories, '
                             'e.g. \'1\' instead of 1.')
        restrict_indices = adata.obs[restrict_key].isin(restrict_categories).values
        adjacency = adjacency[restrict_indices, :]
        adjacency = adjacency[:, restrict_indices]
    if flavor in {'vtraag', 'igraph'}:
        if flavor == 'igraph' and resolution is not None:
            logg.warn('`resolution` parameter has no effect for flavor "igraph"')
        if directed and flavor == 'igraph':
            directed = False
        if not directed: logg.m('    using the undirected graph', v=4)
        g = utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if flavor == 'vtraag':
            import louvain
            if resolution is None: resolution = 1
            try:
                logg.info('    using the "louvain" package of Traag (2017)')
                louvain.set_rng_seed(random_state)
                part = louvain.find_partition(g, louvain.RBConfigurationVertexPartition,
                                              resolution_parameter=resolution)
                # adata.uns['louvain_quality'] = part.quality()
            except AttributeError:
                logg.warn('Did not find package louvain>=0.6, '
                          'the clustering result will therefore not '
                          'be 100% reproducible, '
                          'but still meaningful. '
                          'If you want 100% reproducible results, '
                          'update via "pip install louvain --upgrade".')
                part = louvain.find_partition(g, method='RBConfiguration',
                                              resolution_parameter=resolution)
        elif flavor == 'igraph':
            part = g.community_multilevel()
        groups = np.array(part.membership)
    elif flavor == 'taynaud':
        # this is deprecated
        import networkx as nx
        import community
        g = nx.Graph(adata.uns['data_graph_distance_local'])
        partition = community.best_partition(g)
        groups = np.zeros(len(partition), dtype=int)
        for k, v in partition.items(): groups[k] = v
    else:
        raise ValueError('`flavor` needs to be "vtraag" or "igraph" or "taynaud".')
    unique_groups = np.unique(groups)
    n_clusters = len(unique_groups)
    if restrict_to is None:
        groups = groups.astype('U')
        adata.obs['louvain_groups'] = pd.Categorical(
            values=groups,
            categories=natsorted(unique_groups.astype('U')))
        key_added = 'louvain_groups' if key_added is None else key_added
    else:
        key_added = restrict_key + '_R' if key_added is None else key_added
        groups += 1
        adata.obs[key_added] = adata.obs[restrict_key].astype('U')
        adata.obs[key_added] += ','
        adata.obs[key_added].iloc[restrict_indices] += groups.astype('U')
        adata.obs[key_added].iloc[~restrict_indices] += '0'
        adata.obs[key_added] = adata.obs[key_added].astype(
            'category', categories=natsorted(adata.obs[key_added].unique()))
    adata.uns['louvain_params'] = np.array((resolution, random_state,),
                                           dtype=[('resolution', float), ('random_state', int)])
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('found {} clusters and added\n'
              '    \'{}\', the cluster labels (adata.obs, dtype=category)'
              .format(n_clusters, key_added))
    return adata if copy else None
