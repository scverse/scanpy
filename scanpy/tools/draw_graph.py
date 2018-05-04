import numpy as np
from .. import settings
from .. import utils
from .. import logging as logg


def draw_graph(
        adata,
        layout='fa',
        init_pos=None,
        root=None,
        random_state=0,
        n_jobs=None,
        adjacency=None,
        key_added_ext=None,
        copy=False,
        **kwds):
    """Force-directed graph drawing [Fruchterman91]_ [Islam11]_ [Chippada18]_.

    The default layout ('fa', `ForceAtlas2`) uses the package `fa2
    <https://github.com/bhargavchippada/forceatlas2>`_, which can be installed
    via `pip install fa2`.

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
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    layout : `str`, optional (default: 'fa')
        'fa' (`ForceAtlas2`) or any valid `igraph layout
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
    init_pos : {'paga', any valid 2d-`.obsm` key, `False`}, optional (default: `False`)
        Use precomputed coordinates for initialization.
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
    avail_layouts = {'fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa'}
    if layout not in avail_layouts:
        raise ValueError('Provide a valid layout, one of {}.'.format(avail_layouts))
    adata = adata.copy() if copy else adata
    if adjacency is None and 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    if adjacency is None:
        adjacency = adata.uns['neighbors']['connectivities']
    key_added = 'X_draw_graph_' + (layout if key_added_ext is None else key_added_ext)
    # init coordinates
    np.random.seed(random_state)
    if init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif init_pos == 'paga':
        if 'paga' in adata.uns and 'pos' in adata.uns['paga']:
            groups = adata.obs[adata.uns['paga']['groups']]
            pos = adata.uns['paga']['pos']
            confidence_coarse = adata.uns['paga']['confidence']
            init_coords = np.ones((adjacency.shape[0], 2))
            for i, group_pos in enumerate(pos):
                subset = (groups == groups.cat.categories[i]).values
                neighbors = confidence_coarse[i].nonzero()
                confidence = confidence_coarse[i][neighbors]
                nearest_neighbor = neighbors[1][np.argmax(confidence)]
                noise = np.random.random((len(subset[subset]), 2))
                dist = pos[i] - pos[nearest_neighbor]
                noise = noise * dist
                init_coords[subset] = group_pos - 0.5*dist + noise
        else:
            raise ValueError('Compute and plot PAGA first.')
    else:
        init_coords = np.random.random((adjacency.shape[0], 2))
    # actual drawing
    if layout == 'fa':
        try:
            from fa2 import ForceAtlas2
        except:
            raise ImportError('Please, install package \'fa2\'.')
        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=False,  # Dissuade hubs
            linLogMode=False,  # NOT IMPLEMENTED
            adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
            edgeWeightInfluence=1.0,
            # Performance
            jitterTolerance=1.0,  # Tolerance
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            multiThreaded=False,  # NOT IMPLEMENTED
            # Tuning
            scalingRatio=2.0,
            strongGravityMode=False,
            gravity=1.0,
            # Log
            verbose=False)
        if 'maxiter' in kwds:
            iterations = kwds['maxiter']
        elif 'iterations' in kwds:
            iterations = kwds['iterations']
        else:
            iterations = 500
        positions = forceatlas2.forceatlas2(
            adjacency, pos=init_coords, iterations=iterations)
        positions = np.array(positions)
    else:
        g = utils.get_igraph_from_adjacency(adjacency)
        if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
            ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
        elif 'rt' in layout:
            if root is not None: root = [root]
            ig_layout = g.layout(layout, root=root, **kwds)
        else:
            ig_layout = g.layout(layout, **kwds)
        positions = np.array(ig_layout.coords)
    adata.uns['draw_graph'] = {}
    adata.uns['draw_graph']['params'] = {'layout': layout, 'random_state': random_state}
    key_added = 'X_draw_graph_' + (layout if key_added_ext is None else key_added_ext)
    adata.obsm[key_added] = positions
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'{}\', graph_drawing coordinates (adata.obs)'
              .format(key_added))
    return adata if copy else None
