from __future__ import annotations

import random
from importlib.util import find_spec
from typing import TYPE_CHECKING, Literal

import numpy as np

from .. import _utils
from .. import logging as logg
from .._compat import old_positionals
from .._utils import _choose_graph, get_literal_vals
from ._utils import get_init_pos_from_paga

if TYPE_CHECKING:
    from typing import LiteralString, TypeVar

    from anndata import AnnData

    from .._compat import SpBase
    from .._utils.random import _LegacyRandom

    S = TypeVar("S", bound=LiteralString)


_Layout = Literal["fr", "drl", "kk", "grid_fr", "lgl", "rt", "rt_circular", "fa"]


@old_positionals(
    "init_pos",
    "root",
    "random_state",
    "n_jobs",
    "adjacency",
    "key_added_ext",
    "neighbors_key",
    "obsp",
    "copy",
)
def draw_graph(  # noqa: PLR0913
    adata: AnnData,
    layout: _Layout = "fa",
    *,
    init_pos: str | bool | None = None,
    root: int | None = None,
    random_state: _LegacyRandom = 0,
    n_jobs: int | None = None,
    adjacency: SpBase | None = None,
    key_added_ext: str | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    copy: bool = False,
    **kwds,
) -> AnnData | None:
    """Force-directed graph drawing :cite:p:`Islam2011,Jacomy2014,Chippada2018`.

    An alternative to tSNE that often preserves the topology of the data
    better. This requires running :func:`~scanpy.pp.neighbors`, first.

    The default layout ('fa', `ForceAtlas2`, :cite:t:`Jacomy2014`) uses the package |fa2-modified|_
    :cite:p:`Chippada2018`, which can be installed via `pip install fa2-modified`.

    `Force-directed graph drawing`_ describes a class of long-established
    algorithms for visualizing graphs.
    It was suggested for visualizing single-cell data by :cite:t:`Islam2011`.
    Many other layouts as implemented in igraph :cite:p:`Csardi2006` are available.
    Similar approaches have been used by :cite:t:`Zunder2015` or :cite:t:`Weinreb2017`.

    .. |fa2-modified| replace:: `fa2-modified`
    .. _fa2-modified: https://github.com/AminAlam/fa2_modified
    .. _Force-directed graph drawing: https://en.wikipedia.org/wiki/Force-directed_graph_drawing

    Parameters
    ----------
    adata
        Annotated data matrix.
    layout
        'fa' (`ForceAtlas2`) or any valid `igraph layout
        <https://igraph.org/c/doc/igraph-Layout.html>`__. Of particular interest
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
        If not specified, draw_graph looks at .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, draw_graph looks at
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Return a copy instead of writing to adata.
    **kwds
        Parameters of chosen igraph layout. See e.g.
        :meth:`~igraph.GraphBase.layout_fruchterman_reingold` :cite:p:`Fruchterman1991`.
        One of the most important ones is `maxiter`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsm['X_draw_graph_[layout | key_added_ext]']` : :class:`numpy.ndarray` (dtype `float`)
        Coordinates of graph layout. E.g. for `layout='fa'` (the default),
        the field is called `'X_draw_graph_fa'`. `key_added_ext` overwrites `layout`.
    `adata.uns['draw_graph']`: :class:`dict`
        `draw_graph` parameters.

    """
    start = logg.info(f"drawing single-cell graph using layout {layout!r}")
    if layout not in (layouts := get_literal_vals(_Layout)):
        msg = f"Provide a valid layout, one of {layouts}."
        raise ValueError(msg)
    adata = adata.copy() if copy else adata
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    # init coordinates
    if init_pos in adata.obsm:
        init_coords = adata.obsm[init_pos]
    elif init_pos == "paga" or init_pos:
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
    layout = coerce_fa2_layout(layout)
    # actual drawing
    if layout == "fa":
        positions = np.array(fa2_positions(adjacency, init_coords, **kwds))
    else:
        # igraph doesn't use numpy seed
        random.seed(random_state)

        g = _utils.get_igraph_from_adjacency(adjacency)
        if layout in {"fr", "drl", "kk", "grid_fr"}:
            ig_layout = g.layout(layout, seed=init_coords.tolist(), **kwds)
        elif "rt" in layout:
            if root is not None:
                root = [root]
            ig_layout = g.layout(layout, root=root, **kwds)
        else:
            ig_layout = g.layout(layout, **kwds)
        positions = np.array(ig_layout.coords)
    adata.uns["draw_graph"] = {}
    adata.uns["draw_graph"]["params"] = dict(layout=layout, random_state=random_state)
    key_added = f"X_draw_graph_{key_added_ext or layout}"
    adata.obsm[key_added] = positions
    logg.info(
        "    finished",
        time=start,
        deep=f"added\n    {key_added!r}, graph_drawing coordinates (adata.obsm)",
    )
    return adata if copy else None


def fa2_positions(
    adjacency: SpBase | np.ndarray, init_coords: np.ndarray, **kwds
) -> list[tuple[float, float]]:
    from fa2_modified import ForceAtlas2

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
    if "maxiter" in kwds:
        iterations = kwds["maxiter"]
    elif "iterations" in kwds:
        iterations = kwds["iterations"]
    else:
        iterations = 500
    return forceatlas2.forceatlas2(adjacency, pos=init_coords, iterations=iterations)


def coerce_fa2_layout(layout: S) -> S | Literal["fa", "fr"]:
    # see whether fa2 is installed
    if layout != "fa":
        return layout

    if find_spec("fa2_modified") is None:
        logg.warning(
            "Package 'fa2-modified' is not installed, falling back to layout 'fr'."
            "To use the faster and better ForceAtlas2 layout, "
            "install package 'fa2-modified' (`pip install fa2-modified`)."
        )
        return "fr"

    return "fa"
