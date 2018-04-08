import numpy as np
from .. import settings
from .. import utils
from .. import logging as logg


def draw_graph(
        adata,
        layout='fr',
        root=None,
        random_state=0,
        n_jobs=None,
        adjacency=None,
        key_added_ext=None,
        proceed=False,
        use_paga=False,
        copy=False,
        **kwds):
    """Force-directed graph drawing [Fruchterman91]_ [Islam11]_ [Csardi06]_.

    An alternative to tSNE that often preserves the topology of the data
    better. However, it runs considerably slower. This requires to run
    :func:`~scanpy.api.pp.neighbors`, first.

    `Force-directed graph drawing
    <https://en.wikipedia.org/wiki/Force-directed_graph_drawing>`__ describes a
    class of long-established algorithms for visualizing graphs. It has been
    suggested for visualizing single-cell data by [Islam11]_. Here, by
    default, the Fruchterman & Reingold [Fruchterman91]_ algorithm is used; many
    other layouts are available. Uses the igraph implementation [Csardi06]_.

    Similar approaches have been used by [Zunder15]_ or [Weinreb17]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    layout : `str`, optional (default: 'fr')
        Any valid `igraph layout
        <http://igraph.org/c/doc/igraph-Layout.html>`_. Of particular interest
        are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid Fruchterman Reingold,
        faster than 'fr'), 'kk' (Kamadi Kawai', slower than 'fr'), 'lgl' (Large
        Graph, very fast), 'drl' (Distributed Recursive Layout, pretty fast) and
        'rt' (Reingold Tilford tree layout).
    root : `int` or `None`, optional (default: `None`)
        Root for tree layouts.
    random_state : `int` or `None`, optional (default: 0)
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, no seed is set.
    adjacency : sparse matrix or `None`, optional (default: `None`)
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']`.
    key_ext : `str`, optional (default: `None`)
        By default, append `layout`.
    proceed : `bool`, optional (default: `None`)
        Continue computation, starting off with 'X_draw_graph_`layout`'.
    use_paga : {'local', 'global', `False`}, optional (default: `False`)
        Use the PAGA coordinates.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.
    **kwds : further parameters
        Parameters of chosen igraph layout. See, e.g.,
        `fruchterman_reingold <http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold>`_. One of the most important ones is `maxiter`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following field.

    X_draw_graph_`layout` : `adata.obsm`
        Coordinates of graph layout.
    """
    logg.info('drawing single-cell graph using layout "{}"'.format(layout),
              r=True)
    avail_layouts = {'fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular'}
    if layout not in avail_layouts:
        raise ValueError('Provide a valid layout, one of {}.'.format(avail_layouts))
    adata = adata.copy() if copy else adata
    if adjacency is None and 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    if adjacency is None:
        adjacency = adata.uns['neighbors']['connectivities']
    key_added = 'X_draw_graph_' + (layout if key_added_ext is None else key_added_ext)
    if use_paga != 'local':
        g = utils.get_igraph_from_adjacency(adjacency)
    np.random.seed(random_state)
    all_coords = None
    if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
        if proceed:
            init_coords = adata.obsm[key_added]
            ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
        elif 'paga' in adata.uns and 'pos' in adata.uns['paga']:
            groups = adata.obs[adata.uns['paga']['groups']]
            all_pos = adata.uns['paga']['pos']
            if use_paga == 'global':
                init_coords = np.ones((adjacency.shape[0], 2))
                for i, pos in enumerate(all_pos):
                    subset = (groups == groups.cat.categories[i]).values
                    init_coords[subset] = pos
                ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
            elif use_paga == 'local':
                all_pos -= np.min(all_pos, axis=0)
                all_pos /= np.max(all_pos, axis=0)
                all_coords = np.zeros((adjacency.shape[0], 2))
                from scipy.spatial.distance import euclidean
                for i, pos in enumerate(all_pos):
                    subset = (groups == groups.cat.categories[i]).values
                    adjacency_sub = adjacency[subset][:, subset]
                    if adjacency_sub.shape[0] == 0: continue
                    g = utils.get_igraph_from_adjacency(adjacency_sub)
                    init_coords = np.random.random((adjacency_sub.shape[0], 2))
                    ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
                    coords = np.array(ig_layout.coords)
                    coords -= np.min(coords, axis=0)
                    coords /= np.max(coords, axis=0)
                    dists = [euclidean(pos, p) for j, p in enumerate(all_pos) if i != j]
                    min_dist = np.min(dists)
                    coords -= np.array([0.5, 0.5])
                    coords *= min_dist
                    coords += pos
                    all_coords[subset] = coords
        else:
            init_coords = np.random.random((adjacency.shape[0], 2))
            ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
    elif 'rt' in layout:
        if root is not None: root = [root]
        ig_layout = g.layout(layout, root=root, **kwds)
    else:
        ig_layout = g.layout(layout, **kwds)
    adata.uns['draw_graph'] = {}
    adata.uns['draw_graph']['params'] = {'layout': layout, 'random_state': random_state}
    key_added = 'X_draw_graph_' + (layout if key_added_ext is None else key_added_ext)
    adata.obsm[key_added] = (np.array(ig_layout.coords)
                           if all_coords is None else all_coords)
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'{}\', graph_drawing coordinates (adata.obs)'
              .format(key_added))
    return adata if copy else None
