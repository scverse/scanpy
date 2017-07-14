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
               n_jobs=None,
               random_state=0,
               recompute_graph=False,
               copy=False):
    """Visualize data using standard graph drawing algorithms.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layout : str, optional (default: 'fruchterman_reingold')
        Any valid igraph layout: http://igraph.org/c/doc/igraph-Layout.html.  Of
        particular interest are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid
        Fruchterman Reingold, faster than 'fr'), 'kk' (Kamadi Kawai', slower
        than 'fr'), 'lgl' (Large Graph, very fast), 'drl' (Distributed Recursive
        Layout) and 'tree' (Reingold Tilford').
    n_neighbors : int
        Number of nearest neighbors in graph.

    Returns
    -------
    Returns or updates adata depending on `copy` with
         "X_draw_graph", graph-drawing coordinates (adata.smp)

    References
    ----------
    - The package "igraph", which provides the drawing implementations used
      here: Csardi & Nepusz, InterJournal Complex Systems, 1695 (2006)
    - Suggestion to use the "spring" graph-drawing algorithm of the package D3js
      for single-cell data: Weinreb et al., bioRxiv doi:10.1101/090332 (2016)
    """
    import numpy as np
    import igraph as ig
    from .. import logging as logg
    from .. import data_structs
    adata = adata.copy() if copy else adata
    logg.info('drawing single-cell graph using layout "{}"'.format(layout))
    if 'Ktilde' not in adata.add or recompute_graph:
        graph = data_structs.DataGraph(adata,
                                       k=n_neighbors,
                                       n_pcs=n_pcs,
                                       n_jobs=n_jobs)
        graph.compute_transition_matrix()
        adata.add['Ktilde'] = graph.Ktilde
    adjacency = adata.add['Ktilde']
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    weights = np.array(weights)[0]  # need to convert sparse matrix into a form appropriate for igraph
    g = ig.Graph(list(zip(sources, targets)),  #, edge_attrs={'weight': weights}
                 directed=True)
    if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
        np.random.seed(random_state)
        init_coords = np.random.random((adjacency.shape[0], 2)).tolist()
        ig_layout = g.layout(layout,  # weights='weight',
                             seed=init_coords)
    else:
        ig_layout = g.layout(layout)
    adata.smp['X_draw_graph'] = np.array(ig_layout.coords)
    adata.add['draw_graph_layout'] = layout
    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "X_draw_graph", graph_drawing coordinates (adata.smp)\n'
           '    "draw_graph_layout", the chosen layout (adata.add)')
    return adata if copy else None
