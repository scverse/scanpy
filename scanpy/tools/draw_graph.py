# Author: F. Alex Wolf (http://falexwolf.de)
"""Standard graph drawing for the single-cell graph

References
----------
- General: https://en.wikipedia.org/wiki/Force-directed_graph_drawing
- Suggested for drawing knn-graphs in the context of single-cell
  transcriptomics: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
"""


def draw_graph(adata,
               layout='fr',
               n_neighbors=30,
               n_pcs=50,
               root=None,
               n_jobs=None,
               random_state=0,
               recompute_graph=False,
               adjacency=None,
               copy=False):
    """Visualize data using standard graph drawing algorithms.

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
    import numpy as np
    from .. import data_structs
    from .. import utils
    avail_layouts = {'fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular'}
    if layout not in avail_layouts:
        raise ValueError('Provide a valid layout, one of {}.'.format(avail_layouts))
    adata = adata.copy() if copy else adata
    if 'Ktilde' not in adata.add or recompute_graph:
        graph = data_structs.DataGraph(adata,
                                       k=n_neighbors,
                                       n_pcs=n_pcs,
                                       n_jobs=n_jobs)
        graph.compute_transition_matrix(recompute_distance=True)
        adata.add['Ktilde'] = graph.Ktilde
    elif n_neighbors is not None and not recompute_graph:
        logg.warn('`n_neighbors={}` has no effect (set `recompute_graph=True` to enable it)'
                  .format(n_neighbors))
    adjacency = adata.add['Ktilde']
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
