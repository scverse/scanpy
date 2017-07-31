# Author: F. Alex Wolf (http://falexwolf.de)
"""Graph drawing for the single-cell graph.

References
----------
- General: https://en.wikipedia.org/wiki/Force-directed_graph_drawing
- Suggested for drawing knn-graphs in the context of single-cell
  transcriptomics: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
"""

import numpy as np
from .. import utils
from ..data_structs.data_graph import add_or_update_graph_in_adata


def draw_graph(adata,
               layout='fr',
               root=None,
               n_neighbors=30,
               n_pcs=50,
               random_state=0,
               recompute_pca=False,
               recompute_distances=False,
               recompute_graph=False,
               adjacency=None,
               n_jobs=None,
               copy=False):
    """Visualize data using graph-drawing algorithms.

    Sometimes a good alternative to tSNE, but runs considerably slower.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layout : str, optional (default: 'fr')
        Any valid igraph layout: http://igraph.org/c/doc/igraph-Layout.html. Of
        particular interest are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid
        Fruchterman Reingold, faster than 'fr'), 'kk' (Kamadi Kawai', slower
        than 'fr'), 'lgl' (Large Graph, very fast), 'drl' (Distributed Recursive
        Layout, pretty fast) and 'rt' (Reingold Tilford tree layout).
    n_neighbors : int
        Number of nearest neighbors in graph.
    n_pcs : int
        Number of PCs used to compute distances.

    Returns
    -------
    Returns or updates adata depending on `copy` with
         `"X_draw_graph_" + layout`, the graph-drawing coordinates (adata.smp)

    References
    ----------
    - The package "igraph", which provides the drawing implementations used
      here: Csardi & Nepusz, InterJournal Complex Systems, 1695 (2006)
    - Suggestion to use the "spring" graph-drawing algorithm of the package D3js
      for single-cell data: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """
    from .. import logging as logg
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
    adjacency = adata.add['data_graph_norm_weights']
    g = utils.get_igraph_from_adjacency(adjacency)
    if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
        np.random.seed(random_state)
        init_coords = np.random.random((adjacency.shape[0], 2)).tolist()
        ig_layout = g.layout(layout,  # weights='weight',
                             seed=init_coords)
    elif 'rt' in layout:
        if root is not None: root = [root]
        ig_layout = g.layout(layout, root=root)
    else:
        ig_layout = g.layout(layout)
    if 'draw_graph_layout' in adata.add:
        adata.add['draw_graph_layout'] = list(adata.add['draw_graph_layout']) + [layout]
    else:
        adata.add['draw_graph_layout'] = [layout]
    smp_key = 'X_draw_graph_' + layout
    adata.smp[smp_key] = np.array(ig_layout.coords)
    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "{}", graph_drawing coordinates (adata.smp)\n'
           '    "draw_graph_layout", the chosen layout (adata.add)'
           .format(smp_key))
    return adata if copy else None
