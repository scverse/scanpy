# Author: Alex Wolf (http://falexwolf.de)
"""Graph drawing for the single-cell graph.

References
----------
- General: https://en.wikipedia.org/wiki/Force-directed_graph_drawing
- Suggested for drawing knn-graphs in the context of single-cell
  transcriptomics: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
"""

import numpy as np
from .. import settings
from .. import utils
from .. import logging as logg
from ..data_structs.data_graph import add_or_update_graph_in_adata


def draw_graph(adata,
               layout='fr',
               root=None,
               n_neighbors=None,
               n_pcs=None,
               random_state=0,
               recompute_pca=False,
               recompute_distances=False,
               recompute_graph=False,
               n_jobs=None,
               copy=False,
               **kwargs):
    """Force-directed graph drawing [Fruchterman91]_ [Weinreb17]_ [Csardi06]_.

    Often a good alternative to tSNE, but runs considerably slower.

    `Force-directed graph drawing
    <https://en.wikipedia.org/wiki/Force-directed_graph_drawing>`__ describes a
    class of long-established algorithms for visualizing graphs. It has been
    suggested for visualizing single-cell data by [Weinreb17]_. Here, by
    default, the Fruchterman & Reingold [Fruchterman91]_ algorithm is used; many
    other layouts are available. Uses the igraph implementation [Csardi06]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    layout : `str`, optional (default: 'fr')
        Any valid `igraph layout
        <http://igraph.org/c/doc/igraph-Layout.html>`__. Of particular interest
        are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid Fruchterman Reingold,
        faster than 'fr'), 'kk' (Kamadi Kawai', slower than 'fr'), 'lgl' (Large
        Graph, very fast), 'drl' (Distributed Recursive Layout, pretty fast) and
        'rt' (Reingold Tilford tree layout).
    n_neighbors : `int` or `None` (default: `None`)
        Number of nearest neighbors in graph.
    n_pcs : `int` or `None` (default: `None`)
        Number of PCs used to compute distances.
    random_state : `int` or `None`, optional (default: 0)
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, no seed is set.
    **kwargs : further parameters
        Parameters of chosen igraph algorithm. See, e.g.,
        http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold.
    n_jobs : `int` or `None` (default: `sc.settings.n_jobs`)
        Number of jobs.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_draw_graph_`layout` : `np.ndarray` (`adata.smpm`, dtype `float`)
        Coordinates of graph layout.
    """
    logg.info('drawing single-cell graph using layout "{}"'.format(layout),
              r=True)
    avail_layouts = {'fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular'}
    if layout not in avail_layouts:
        raise ValueError('Provide a valid layout, one of {}.'.format(avail_layouts))
    adata = adata.copy() if copy else adata
    add_or_update_graph_in_adata(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        recompute_pca=recompute_pca,
        recompute_distances=recompute_distances,
        recompute_graph=recompute_graph,
        n_jobs=n_jobs)
    adjacency = adata.uns['data_graph_norm_weights']
    g = utils.get_igraph_from_adjacency(adjacency)
    if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
        np.random.seed(random_state)
        init_coords = np.random.random((adjacency.shape[0], 2)).tolist()
        ig_layout = g.layout(layout,  # weights='weight',
                             seed=init_coords, **kwargs)
    elif 'rt' in layout:
        if root is not None: root = [root]
        ig_layout = g.layout(layout, root=root, **kwargs)
    else:
        ig_layout = g.layout(layout, **kwargs)
    adata.uns['draw_graph_params'] = np.array(
        (layout, random_state,),
        dtype=[('layout', 'U20'), ('random_state', int)])
    smp_key = 'X_draw_graph_' + layout
    adata.smpm[smp_key] = np.array(ig_layout.coords)
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'{}\', graph_drawing coordinates (adata.smp)\n'
              '    \'draw_graph_params\', the parameters (adata.uns)'
              .format(smp_key))
    return adata if copy else None
