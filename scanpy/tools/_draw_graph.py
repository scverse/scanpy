from typing import Union, Optional

import numpy as np
import random
from anndata import AnnData
from scipy.sparse import spmatrix

from .. import _utils
from .. import logging as logg
from ._utils import get_init_pos_from_paga
from .._compat import Literal
from .._utils import AnyRandom, _choose_graph


_LAYOUTS = ('fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa')
_Layout = Literal[_LAYOUTS]


def draw_graph(
    adata: AnnData,
    layout: _Layout = 'fa',
    init_pos: Union[str, bool, None] = None,
    root: Optional[int] = None,
    random_state: AnyRandom = 0,
    n_jobs: Optional[int] = None,
    adjacency: Optional[spmatrix] = None,
    key_added_ext: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    obsp: Optional[str] = None,
    copy: bool = False,
    **kwds,
):
    """\
    Force-directed graph drawing [Islam11]_ [Jacomy14]_ [Chippada18]_.

    An alternative to tSNE that often preserves the topology of the data
    better. This requires to run :func:`~scanpy.pp.neighbors`, first.

    The default layout ('fa', `ForceAtlas2`) [Jacomy14]_ uses the package |fa2|_
    [Chippada18]_, which can be installed via `pip install fa2`.

    `Force-directed graph drawing`_ describes a class of long-established
    algorithms for visualizing graphs.
    It has been suggested for visualizing single-cell data by [Islam11]_.
    Many other layouts as implemented in igraph [Csardi06]_ are available.
    Similar approaches have been used by [Zunder15]_ or [Weinreb17]_.

    .. |fa2| replace:: `fa2`
    .. _fa2: https://github.com/bhargavchippada/forceatlas2
    .. _Force-directed graph drawing: https://en.wikipedia.org/wiki/Force-directed_graph_drawing

    Parameters
    ----------
    adata
        Annotated data matrix.
    layout
        'fa' (`ForceAtlas2`) or any valid `igraph layout
        <http://igraph.org/c/doc/igraph-Layout.html>`__. Of particular interest
        are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid Fruchterman Reingold,
        faster than 'fr'), 'kk' (Kamadi Kawai', slower than 'fr'), 'lgl' (Large
        Graph, very fast), 'drl' (Distributed Recursive Layout, pretty fast) and
        'rt' (Reingold Tilford tree layout).
    root
        Root for tree layouts.
    random_state
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, no seed is set.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    key_added_ext
        By default, append `layout`.
    proceed
        Continue computation, starting off with 'X_draw_graph_`layout`'.
    init_pos
        `'paga'`/`True`, `None`/`False`, or any valid 2d-`.obsm` key.
        Use precomputed coordinates for initialization.
        If `False`/`None` (the default), initialize randomly.
    neighbors_key
        If not specified, draw_graph looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, draw_graph looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Return a copy instead of writing to adata.
    **kwds
        Parameters of chosen igraph layout. See e.g. `fruchterman-reingold`_
        [Fruchterman91]_. One of the most important ones is `maxiter`.

        .. _fruchterman-reingold: http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following field.

    **X_draw_graph_layout** : `adata.obsm`
        Coordinates of graph layout. E.g. for layout='fa' (the default),
        the field is called 'X_draw_graph_fa'
    """
    start = logg.info(f'drawing single-cell graph using layout {layout!r}')
    if layout not in _LAYOUTS:
        raise ValueError(f'Provide a valid layout, one of {_LAYOUTS}.')
    adata = adata.copy() if copy else adata
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    # init coordinates
    if init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif init_pos == 'paga' or init_pos:
        init_coords = get_init_pos_from_paga(
            adata,
            adjacency,
            random_state=random_state,
            neighbors_key=neighbors_key,
            obsp=obsp,
        )
    else:
        np.random.seed(random_state)
        init_coords = np.random.random((adjacency.shape[0], 2))
    # see whether fa2 is installed
    if layout == 'fa':
        try:
            from fa2 import ForceAtlas2
        except ImportError:
            logg.warning(
                "Package 'fa2' is not installed, falling back to layout 'fr'."
                "To use the faster and better ForceAtlas2 layout, "
                "install package 'fa2' (`pip install fa2`)."
            )
            layout = 'fr'
    # actual drawing
    if layout == 'fa':
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
            verbose=False,
        )
        if 'maxiter' in kwds:
            iterations = kwds['maxiter']
        elif 'iterations' in kwds:
            iterations = kwds['iterations']
        else:
            iterations = 500
        positions = forceatlas2.forceatlas2(
            adjacency, pos=init_coords, iterations=iterations
        )
        positions = np.array(positions)
    else:
        # igraph doesn't use numpy seed
        random.seed(random_state)

        g = _utils.get_igraph_from_adjacency(adjacency)
        if layout in {'fr', 'drl', 'kk', 'grid_fr'}:
            ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
        elif 'rt' in layout:
            if root is not None:
                root = [root]
            ig_layout = g.layout(layout, root=root, **kwds)
        else:
            ig_layout = g.layout(layout, **kwds)
        positions = np.array(ig_layout.coords)
    adata.uns['draw_graph'] = {}
    adata.uns['draw_graph']['params'] = dict(layout=layout, random_state=random_state)
    key_added = f'X_draw_graph_{key_added_ext or layout}'
    adata.obsm[key_added] = positions
    logg.info(
        '    finished',
        time=start,
        deep=f'added\n    {key_added!r}, graph_drawing coordinates (adata.obsm)',
    )
    return adata if copy else None
