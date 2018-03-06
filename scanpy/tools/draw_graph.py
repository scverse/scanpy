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
        key='neighbors_distances',
        key_ext='LAYOUT',
        copy=False,
        **kwargs):
    """Force-directed graph drawing [Fruchterman91]_ [Islam11]_ [Csardi06]_.

    Often a good alternative to tSNE, but runs considerably slower.

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
    key : `str`, optional (default: 'neighbors_similarities')
        Key for accessing the sparse adjacency matrix of the graph in
        `adata.uns`.
    key_ext : `str`, optional (default: 'LAYOUT')
        By default, append `layout`.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.
    **kwargs : further parameters
        Parameters of chosen igraph layout. See, e.g.,
        `fruchterman_reingold <http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold>`_.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_draw_graph_`layout` : `adata.obsm`
        Coordinates of graph layout.
    """
    logg.info('drawing single-cell graph using layout "{}"'.format(layout),
              r=True)
    avail_layouts = {'fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular'}
    if layout not in avail_layouts:
        raise ValueError('Provide a valid layout, one of {}.'.format(avail_layouts))
    adata = adata.copy() if copy else adata
    if key not in adata.uns:
        raise ValueError(
            '\'{}\' is not present in `adata.uns`. '
            'You need to run `pp.neighbors` first to compute a neighborhood graph.'
            .format(key))
    adjacency = adata.uns[key]
    g = utils.get_igraph_from_adjacency(adjacency)
    if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
        np.random.seed(random_state)
        init_coords = np.random.random((adjacency.shape[0], 2)).tolist()
        ig_layout = g.layout(layout, seed=init_coords, weights='weight', **kwargs)
    elif 'rt' in layout:
        if root is not None: root = [root]
        ig_layout = g.layout(layout, root=root, **kwargs)
    else:
        ig_layout = g.layout(layout, **kwargs)
    adata.uns['draw_graph_params'] = np.array(
        (layout, random_state,),
        dtype=[('layout', 'U20'), ('random_state', int)])
    obs_key = 'X_draw_graph_' + (layout if key_ext == 'LAYOUT' else key_ext)
    adata.obsm[obs_key] = np.array(ig_layout.coords)
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'{}\', graph_drawing coordinates (adata.obs)\n'
              '    \'draw_graph_params\', the parameters (adata.uns)'
              .format(obs_key))
    return adata if copy else None
