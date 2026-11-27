from __future__ import annotations

import sys
import warnings
from collections.abc import Collection, Mapping, Sequence
from contextlib import nullcontext
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

if sys.version_info < (3, 15):
    from types import MappingProxyType as frozendict  # noqa: N813

import numpy as np
import pandas as pd
import scipy
from matplotlib import patheffects, rcParams, ticker
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from pandas.api.types import CategoricalDtype
from scverse_misc import Deprecation, deprecated

from scanpy.tools._draw_graph import coerce_fa2_layout, fa2_positions

from .... import _utils as _sc_utils
from .... import logging as logg
from ...._compat import CSBase
from ...._docs import doc_rng
from ...._settings import settings
from ...._utils import _doc_params
from ...._utils.random import (
    _accepts_legacy_random_state,
    _LegacyRng,
    _set_igraph_rng,
)
from .. import _utils, mpl_settings
from .._utils import matrix

if TYPE_CHECKING:
    from typing import Any, Literal

    import networkx as nx
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.collections import PathCollection
    from matplotlib.colors import Colormap
    from numpy.typing import NDArray

    from ...._utils.random import RNGLike, SeedLike
    from ....tools._draw_graph import _Layout as _LayoutWithoutEqTree
    from .._utils import _FontSize, _FontWeight, _LegendLoc

type _Layout = _LayoutWithoutEqTree | Literal["eq_tree"]


def paga_compare(  # noqa: PLR0913
    adata: AnnData,
    basis=None,
    *,
    edges=False,
    color=None,
    alpha=None,
    groups=None,
    components=None,
    projection: Literal["2d", "3d"] = "2d",
    legend_loc: _LegendLoc | None = "on data",
    legend_fontsize: float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight = "bold",
    legend_fontoutline=None,
    color_map=None,
    palette=None,
    frameon=False,
    size=None,
    title=None,
    right_margin=None,
    left_margin=0.05,
    show=None,
    save=None,
    title_graph=None,
    groups_graph=None,
    pos=None,
    **paga_graph_params,
):
    """Scatter and PAGA graph side-by-side.

    Consists in a scatter plot and the abstracted graph. See
    :func:`~scanpy.pl.paga` for all related parameters.

    See :func:`~scanpy.pl.paga_path` for visualizing gene changes along paths
    through the abstracted graph.

    Additional parameters are as follows.

    Parameters
    ----------
    adata
        Annotated data matrix.
    kwds_scatter
        Keywords for :func:`~scanpy.pl.scatter`.
    kwds_paga
        Keywords for :func:`~scanpy.pl.paga`.

    Returns
    -------
    A list of :class:`~matplotlib.axes.Axes` if `show` is `False`.

    Examples
    --------
    Compute a PAGA graph on the bundled PBMC dataset and show it next to the UMAP embedding.

    ..  exec-jupyter::

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.paga(adata, groups="bulk_labels")
        sc.pl.paga_compare(adata, basis="umap")

    """
    axs, _, _, _ = _utils.setup_axes(panels=[0, 1], right_margin=right_margin)
    if color is None:
        color = adata.uns["paga"]["groups"]
    suptitle = None  # common title for entire figure
    if title_graph is None:
        suptitle = color if title is None else title
        title, title_graph = "", ""
    if basis is None:
        basis = _default_compare_basis(adata)

    from .scatterplots import embedding

    embedding(
        adata,
        ax=axs[0],
        basis=basis,
        color=color,
        edges=edges,
        alpha=alpha,
        groups=groups,
        components=components,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        legend_fontoutline=legend_fontoutline,
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        size=size,
        title=title,
        show=False,
        save=False,
    )

    if pos is None:
        pos = _compare_pos(adata, basis=basis, color=color, components=components)
    xlim, ylim = axs[0].get_xlim(), axs[0].get_ylim()
    axs[1].set_xlim(xlim)
    axs[1].set_ylim(ylim)
    labels = paga_graph_params.pop("labels", groups_graph)
    font_params = dict(
        fontsize=legend_fontsize,
        fontweight=legend_fontweight,
        fontoutline=legend_fontoutline,
    )
    paga_graph_params.update({k: v for k, v in font_params.items() if v is not None})
    paga(
        adata,
        ax=axs[1],
        show=False,
        save=False,
        title=title_graph,
        labels=labels,
        colors=color,
        frameon=frameon,
        pos=pos,
        **paga_graph_params,
    )
    if suptitle is not None:
        plt.suptitle(suptitle)
    _utils.savefig_or_show("paga_compare", show=show, save=save)
    if show:
        return None
    return axs


def _default_compare_basis(adata: AnnData) -> str:
    for basis in ("draw_graph_fa", "umap", "tsne", "draw_graph_fr"):
        if f"X_{basis}" in adata.obsm:
            return basis
    return "umap"


def _compare_pos(
    adata: AnnData, *, basis: str, color, components
) -> NDArray[np.float64]:
    if color != adata.uns["paga"]["groups"]:
        return adata.uns["paga"]["pos"]

    from .scatterplots import _components_to_dimensions, _get_basis_arr

    # TODO: Use dimensions here
    _basis = _get_basis_arr(adata, basis)
    dims = _components_to_dimensions(
        components=components, dimensions=None, total_dims=_basis.shape[1]
    )[0]
    coords = _basis[:, dims]
    return (
        pd
        .DataFrame(coords, columns=["x", "y"], index=adata.obs_names)
        .groupby(adata.obs[color], observed=True)
        .median()
        .sort_index()
    ).to_numpy()


def _compute_pos(
    adjacency_solid: CSBase | np.ndarray,
    *,
    layout: _Layout | None,
    rng: np.random.Generator,
    init_pos: np.ndarray | None,
    adj_tree,
    root: int,
    layout_kwds: Mapping[str, Any] = frozendict({}),
) -> NDArray[np.float64]:
    import networkx as nx

    nx_g_solid = nx.Graph(adjacency_solid)
    if layout is None:
        layout = "fr"
    layout = coerce_fa2_layout(layout)
    if layout == "fa":
        if init_pos is None:
            init_coords = rng.random((adjacency_solid.shape[0], 2))
        else:
            init_coords = init_pos.copy()
        pos_list = fa2_positions(adjacency_solid, init_coords, **layout_kwds)
        pos = {n: (x, -y) for n, (x, y) in enumerate(pos_list)}
    elif layout == "eq_tree":
        nx_g_tree = nx.Graph(adj_tree)
        pos = hierarchy_pos(nx_g_tree, root)
        if len(pos) < adjacency_solid.shape[0]:
            msg = (
                "This is a forest and not a single tree. "
                "Try another `layout`, e.g., {'fr'}."
            )
            raise ValueError(msg)
    else:  # igraph layouts
        pos_list = _igraph_layout_coords(
            adjacency_solid,
            layout=layout,
            rng=rng,
            init_pos=init_pos,
            adj_tree=adj_tree,
            root=root,
            layout_kwds=layout_kwds,
        )
        pos = {n: (x, -y) for n, (x, y) in enumerate(pos_list)}
    if len(pos) == 1:
        pos[0] = (0.5, 0.5)
    pos_array = np.array([pos[n] for n in nx_g_solid])
    return pos_array


def _igraph_layout_coords(
    adjacency_solid: CSBase | np.ndarray,
    *,
    layout: str,
    rng: np.random.Generator,
    init_pos: np.ndarray | None,
    adj_tree,
    root: int,
    layout_kwds: Mapping[str, Any],
) -> list[list[float]]:
    import random

    if isinstance(rng, _LegacyRng):  # backwards compat
        random.seed(rng.bytes(8))
        ctx = nullcontext()
    else:
        ctx = _set_igraph_rng(rng)
    g = _sc_utils.get_igraph_from_adjacency(adjacency_solid, directed=False)
    with ctx:
        if "rt" in layout:
            g_tree = _sc_utils.get_igraph_from_adjacency(adj_tree)
            return g_tree.layout(
                layout, root=root if isinstance(root, list) else [root]
            ).coords
        if layout == "circle":
            return g.layout(layout).coords
        # I don't know why this is necessary
        if init_pos is None:
            init_coords = rng.random((adjacency_solid.shape[0], 2)).tolist()
        else:
            init_pos = init_pos.copy()
            # this is a super-weird hack that is necessary as igraph’s
            # layout function seems to do some strange stuff here
            init_pos[:, 1] *= -1
            init_coords = init_pos.tolist()
        try:
            return g.layout(
                layout, seed=init_coords, weights="weight", **layout_kwds
            ).coords
        except AttributeError:  # hack for empty graphs...
            return g.layout(layout, seed=init_coords, **layout_kwds).coords


class _Level(TypedDict):
    total: int
    current: int


def hierarchy_pos(
    g, /, root: int, levels_: Mapping[int, int] | None = None, width=1.0, height=1.0
) -> dict[int, tuple[float, float]]:
    """Tree layout for networkx graph.

    See https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
    answer by burubum.

    If there is a cycle that is reachable from root, then this will see
    infinite recursion.

    Parameters
    ----------
    G: the graph
    root: the root node
    levels: a dictionary
            key: level number (starting from 0)
            value: number of nodes in this level
    width: horizontal space allocated for drawing
    height: vertical space allocated for drawing

    """

    def make_levels(
        levels: dict[int, _Level],
        node: int = root,
        current_level: int = 0,
        parent: int | None = None,
    ) -> dict[int, _Level]:
        """Compute the number of nodes for each level."""
        if current_level not in levels:
            levels[current_level] = _Level(total=0, current=0)
        levels[current_level]["total"] += 1
        neighbors: list[int] = list(g.neighbors(node))
        if parent is not None:
            neighbors.remove(parent)
        for neighbor in neighbors:
            levels = make_levels(levels, neighbor, current_level + 1, node)
        return levels

    if levels_ is None:
        levels = make_levels({})
    else:
        levels = {k: _Level(total=0, current=0) for k, v in levels_.items()}

    def make_pos(
        pos: dict[int, tuple[float, float]],
        node: int = root,
        current_level: int = 0,
        parent: int | None = None,
        vert_loc: float = 0.0,
    ):
        dx = 1 / levels[current_level]["total"]
        left = dx / 2
        pos[node] = ((left + dx * levels[current_level]["current"]) * width, vert_loc)
        levels[current_level]["current"] += 1
        neighbors: list[int] = list(g.neighbors(node))
        if parent is not None:
            neighbors.remove(parent)
        for neighbor in neighbors:
            pos = make_pos(pos, neighbor, current_level + 1, node, vert_loc - vert_gap)
        return pos

    vert_gap = height / (max(levels.keys()) + 1)
    return make_pos({})


@_doc_params(rng=doc_rng)
@_accepts_legacy_random_state(0)
def paga(  # noqa: PLR0913
    adata: AnnData,
    *,
    threshold: float | None = None,
    color: str | Mapping[str | int, Mapping[Any, float]] | None = None,
    layout: _Layout | None = None,
    layout_kwds: Mapping[str, Any] = frozendict({}),
    init_pos: np.ndarray | None = None,
    root: int | str | Sequence[int] | None = 0,
    labels: str | Sequence[str] | Mapping[str, str] | None = None,
    single_component: bool = False,
    solid_edges: str = "connectivities",
    dashed_edges: str | None = None,
    transitions: str | None = None,
    fontsize: int | None = None,
    fontweight: str = "bold",
    fontoutline: int | None = None,
    text_kwds: Mapping[str, Any] = frozendict({}),
    node_size_scale: float = 1.0,
    node_size_power: float = 0.5,
    edge_width_scale: float = 1.0,
    min_edge_width: float | None = None,
    max_edge_width: float | None = None,
    arrowsize: int = 30,
    title: str | None = None,
    left_margin: float = 0.01,
    rng: SeedLike | RNGLike | None = None,
    pos: np.ndarray | Path | str | None = None,
    normalize_to_color: bool = False,
    cmap: str | Colormap | None = None,
    cax: Axes | None = None,
    colorbar=None,  # TODO: this seems to be unused
    cb_kwds: Mapping[str, Any] = frozendict({}),
    ncols: int | None = None,
    frameon: bool | None = None,
    add_pos: bool = True,
    export_to_gexf: bool = False,
    use_raw: bool = True,
    colors=None,  # backwards compat
    groups=None,  # backwards compat
    plot: bool = True,
    show: bool | None = None,
    ax: Axes | None = None,
    # deprecated
    save: bool | str | None = None,
) -> Axes | list[Axes] | None:
    r"""Plot the PAGA graph through thresholding low-connectivity edges.

    Compute a coarse-grained layout of the data. Reuse this by passing
    `init_pos='paga'` to :func:`~scanpy.tl.umap` or
    :func:`~scanpy.tl.draw_graph` and obtain embeddings with more meaningful
    global topology :cite:p:`Wolf2019`.

    This uses ForceAtlas2 or igraph's layout algorithms for most layouts :cite:p:`Csardi2006`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    threshold
        Do not draw edges for weights below this threshold. Set to 0 if you want
        all edges. Discarding low-connectivity edges helps in getting a much
        clearer picture of the graph.
    color
        Gene name or `obs` annotation defining the node colors.
        Also plots the degree of the abstracted graph when
        passing {{`'degree_dashed'`, `'degree_solid'`}}.

        Can be also used to visualize pie chart at each node in the following form:
        `{{<group name or index>: {{<color>: <fraction>, ...}}, ...}}`. If the fractions
        do not sum to 1, a new category called `'rest'` colored grey will be created.
    labels
        The node labels. If `None`, this defaults to the group labels stored in
        the categorical for which :func:`~scanpy.tl.paga` has been computed.
    pos
        Two-column array-like storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    layout
        Plotting layout that computes positions.
        `'fa'` stands for “ForceAtlas2”,
        `'fr'` stands for “Fruchterman-Reingold”,
        `'rt'` stands for “Reingold-Tilford”,
        `'eq_tree'` stands for “eqally spaced tree”.
        All but `'fa'` and `'eq_tree'` are igraph layouts.
        All other igraph layouts are also permitted.
        See also parameter `pos` and :func:`~scanpy.tl.draw_graph`.
    layout_kwds
        Keywords for the layout.
    init_pos
        Two-column array storing the x and y coordinates for initializing the
        layout.
    {rng}

        Applies to layouts with random initialization like `'fr'`.
    root
        If choosing a tree layout, this is the index of the root node or a list
        of root node indices. If this is a non-empty vector then the supplied
        node IDs are used as the roots of the trees (or a single tree if the
        graph is connected). If this is `None` or an empty list, the root
        vertices are automatically calculated based on topological sorting.
    transitions
        Key for `.uns['paga']` that specifies the matrix that stores the
        arrows, for instance `'transitions_confidence'`.
    solid_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    single_component
        Restrict to largest connected component.
    fontsize
        Font size for node labels.
    fontoutline
        Width of the white outline around fonts.
    text_kwds
        Keywords for :meth:`~matplotlib.axes.Axes.text`.
    node_size_scale
        Increase or decrease the size of the nodes.
    node_size_power
        The power with which groups sizes influence the radius of the nodes.
    edge_width_scale
        Edge with scale in units of `rcParams['lines.linewidth']`.
    min_edge_width
        Min width of solid edges.
    max_edge_width
        Max width of solid and dashed edges.
    arrowsize
       For directed graphs, choose the size of the arrow head head's length and
       width. See :py:class: `matplotlib.patches.FancyArrowPatch` for attribute
       `mutation_scale` for more info.
    export_to_gexf
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    normalize_to_color
        Whether to normalize categorical plots to `color` or the underlying
        grouping.
    cmap
        The color map.
    cax
        A matplotlib axes object for a potential colorbar.
    cb_kwds
        Keyword arguments for :class:`~matplotlib.colorbar.Colorbar`,
        for instance, `ticks`.
    ncols
        Number of panels shown per row when plotting multiple colors. If
        `None`, all panels are shown in a single row.
    add_pos
        Add the positions to `adata.uns['paga']`.
    title
        Provide a title.
    frameon
        Draw a frame around the PAGA graph.
    plot
        If `False`, do not create the figure, simply compute the layout.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.
    ax
        A matplotlib axes object.

    Returns
    -------
    If `show==False`, one or more :class:`~matplotlib.axes.Axes` objects.
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.

    Examples
    --------

    ..  exec-jupyter::

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata)

    You can increase node and edge sizes by specifying additional arguments.

    ..  exec-jupyter::

        sc.pl.paga(adata, node_size_scale=10, edge_width_scale=2)

    Notes
    -----
    When initializing the positions, note that – for some reason – igraph
    mirrors coordinates along the x axis... that is, you should increase the
    `maxiter` parameter by 1 if the layout is flipped.

    .. currentmodule:: scanpy

    See Also
    --------
    tl.paga
    pl.paga_compare
    pl.paga_path

    """
    rng = np.random.default_rng(rng)
    _validate_ncols(ncols, ax)
    if groups is not None:  # backwards compat
        labels = groups
        logg.warning("`groups` is deprecated in `pl.paga`: use `labels` instead")
    if colors is None:
        colors = color

    colors, labels, title = _normalize_paga_panels(
        adata, colors=colors, labels=labels, title=title
    )
    if frameon is None:
        frameon = mpl_settings.FRAMEON

    colorbars = _paga_colorbars(adata, colors, colorbar=colorbar)
    root = _resolve_root(root, labels)

    # define the adjacency matrices
    threshold = 0.01 if threshold is None else threshold  # default threshold
    adjacency_solid = _thresholded_adjacency(adata, solid_edges, threshold)
    adjacency_dashed = (
        None
        if dashed_edges is None
        else _thresholded_adjacency(adata, dashed_edges, threshold)
    )

    # compute positions
    if pos is None:
        adj_tree = (
            adata.uns["paga"]["connectivities_tree"]
            if layout in {"rt", "rt_circular", "eq_tree"}
            else None
        )
        pos = _compute_pos(
            adjacency_solid,
            layout=layout,
            rng=rng,
            init_pos=init_pos,
            layout_kwds=layout_kwds,
            adj_tree=adj_tree,
            root=root,
        )

    if plot:
        panels = _PagaPanels.create(ax, colors=colors, colorbars=colorbars, ncols=ncols)
        panel_graph = partial(
            _PagaGraph,
            adata=adata,
            adjacency_solid=adjacency_solid,
            adjacency_dashed=adjacency_dashed,
            pos=pos,
            transitions=transitions,
            threshold=threshold,
            fontsize=fontsize,
            fontweight=fontweight,
            fontoutline=fontoutline,
            text_kwds=text_kwds,
            node_size_scale=node_size_scale,
            node_size_power=node_size_power,
            edge_width_scale=edge_width_scale,
            min_edge_width=min_edge_width,
            max_edge_width=max_edge_width,
            normalize_to_color=normalize_to_color,
            frameon=frameon,
            cmap=cmap,
            use_raw=use_raw,
            export_to_gexf=export_to_gexf,
            single_component=single_component,
            arrowsize=arrowsize,
        )
        for icolor, c in enumerate(colors):
            if title[icolor] is not None:
                panels.axs[icolor].set_title(title[icolor])
            graph = panel_graph(
                ax=panels.axs[icolor],
                colors=colors if isinstance(colors, Mapping) else c,
                labels=labels[icolor],
            )
            sct = graph.draw()
            if colorbars[icolor]:
                panels.add_colorbar(sct, icolor, cax=cax, cb_kwds=cb_kwds)
        axs = panels.axs
    if add_pos:
        adata.uns["paga"]["pos"] = pos
        logg.hint("added 'pos', the PAGA positions (adata.uns['paga'])")

    if not plot:
        return None
    _utils.savefig_or_show("paga", show=show, save=save)
    if len(colors) == 1 and isinstance(axs, list):
        axs = axs[0]
    show = settings.autoshow if show is None else show
    if show:
        return None
    return axs


def _validate_ncols(ncols: int | None, ax: Axes | None) -> None:
    if ncols is None:
        return
    if isinstance(ncols, bool) or not isinstance(ncols, int) or ncols < 1:
        msg = "`ncols` must be a positive integer or `None`."
        raise ValueError(msg)
    if ax is not None:
        msg = "Cannot specify `ncols` when passing a pre-supplied `ax`."
        raise ValueError(msg)


def _normalize_paga_panels(adata: AnnData, *, colors, labels, title) -> tuple:
    """Turn `colors`, `labels`, and `title` into one entry per panel."""
    groups_key = adata.uns["paga"]["groups"]

    def is_flat(x):
        has_one_per_category = isinstance(x, Collection) and len(x) == len(
            adata.obs[groups_key].cat.categories
        )
        return has_one_per_category or x is None or isinstance(x, str)

    if isinstance(colors, Mapping) and isinstance(colors[next(iter(colors))], Mapping):
        # handle paga pie, remap string keys to integers
        names_to_ixs = {
            n: i for i, n in enumerate(adata.obs[groups_key].cat.categories)
        }
        colors = {names_to_ixs.get(n, n): v for n, v in colors.items()}
    if is_flat(colors):
        colors = [colors]

    # labels is a list that contains no lists
    if is_flat(labels):
        labels = [labels for _ in range(len(colors))]

    if title is None and len(colors) > 1:
        title = list(colors)
    elif isinstance(title, str):
        title = [title] * len(colors)
    elif title is None:
        title = [None] * len(colors)
    return colors, labels, title


def _paga_colorbars(adata: AnnData, colors, *, colorbar) -> list[bool]:
    if colorbar is not None:
        return [False for _ in colors]
    var_names = adata.var_names if adata.raw is None else adata.raw.var_names
    return [
        ((c in adata.obs and adata.obs[c].dtype.name != "category") or (c in var_names))
        for c in colors
    ]


def _resolve_root(root: int | str | Sequence[int] | None, labels):
    if isinstance(root, str):
        if root not in labels:
            msg = f"If `root` is a string, it needs to be one of {labels} not {root!r}."
            raise ValueError(msg)
        root = list(labels).index(root)
    if isinstance(root, Sequence) and root[0] in labels:
        root = [list(labels).index(r) for r in root]
    return root


def _thresholded_adjacency(adata: AnnData, key: str, threshold: float):
    adjacency = adata.uns["paga"][key].copy()
    if threshold > 0:
        adjacency.data[adjacency.data < threshold] = 0
        adjacency.eliminate_zeros()
    return adjacency


@dataclass(kw_only=True)
class _PagaPanels:
    """Axes for the panels of :func:`paga`, and their colorbar placement."""

    axs: list[Axes]
    n_panels: int
    ncols: int | None
    panel_pos: Any = None
    draw_region_width: float | None = None

    @classmethod
    def create(
        cls, ax: Axes | None, *, colors, colorbars: list[bool], ncols: int | None
    ) -> _PagaPanels:
        if ncols is not None:
            from .scatterplots import _panel_grid

            fig, grid = _panel_grid(
                hspace=0.25, wspace=0.1, ncols=ncols, num_panels=len(colors)
            )
            axs = [fig.add_subplot(grid[i]) for i in range(len(colors))]
            return cls(axs=axs, n_panels=len(colors), ncols=ncols)

        axs, panel_pos, draw_region_width, _figure_width = _utils.setup_axes(
            ax, panels=colors, colorbars=colorbars
        )
        if len(colors) == 1 and not isinstance(axs, list):
            axs = [axs]
        return cls(
            axs=axs,
            n_panels=len(colors),
            ncols=ncols,
            panel_pos=panel_pos,
            draw_region_width=draw_region_width,
        )

    def add_colorbar(
        self,
        sct: PathCollection,
        icolor: int,
        *,
        cax: Axes | Sequence[Axes] | None,
        cb_kwds: Mapping[str, Any],
    ) -> None:
        if self.ncols is not None and cax is None:
            colorbar_kwds = {"format": ticker.FuncFormatter(_utils.ticks_formatter)}
            colorbar_kwds.update(cb_kwds)
            ax = self.axs[icolor]
            _ = ax.figure.colorbar(sct, ax=ax, use_gridspec=False, **colorbar_kwds)
            return

        if self.ncols is not None:
            ax_cb = cax[icolor] if isinstance(cax, Sequence) else cax
        elif cax is None:
            ax_cb = self._add_colorbar_axes(icolor)
        else:
            ax_cb = cax[icolor]
        _ = plt.colorbar(
            sct, format=ticker.FuncFormatter(_utils.ticks_formatter), cax=ax_cb
        )

    def _add_colorbar_axes(self, icolor: int) -> Axes:
        bottom = self.panel_pos[0][0]
        height = self.panel_pos[1][0] - bottom
        width = 0.006 * self.draw_region_width / self.n_panels
        left = self.panel_pos[2][2 * icolor + 1] + 0.2 * width
        return plt.gcf().add_axes([left, bottom, width, height])


@dataclass(kw_only=True)
class _PagaGraph:
    """The PAGA graph in a single panel, colored by one color specification."""

    adata: AnnData
    ax: Axes
    colors: Any
    labels: Any
    adjacency_solid: CSBase
    adjacency_dashed: CSBase | None
    pos: np.ndarray | Path | str
    transitions: str | None
    threshold: float | None
    fontsize: int | None
    fontweight: str | None
    fontoutline: int | None
    text_kwds: Mapping[str, Any]
    node_size_scale: float
    node_size_power: float
    edge_width_scale: float
    min_edge_width: float | None
    max_edge_width: float | None
    normalize_to_color: bool
    frameon: bool
    cmap: str | Colormap | None
    use_raw: bool
    export_to_gexf: bool
    single_component: bool
    arrowsize: int

    groups_key: str = field(init=False)
    node_labels: Any = field(init=False)
    nx_g_solid: nx.Graph = field(init=False)
    nx_g_dashed: nx.Graph | None = field(init=False)
    pos_array: np.ndarray = field(init=False)
    pos_dict: dict[int, list[float]] = field(init=False)

    def __post_init__(self) -> None:
        import networkx as nx

        self.groups_key = self.adata.uns["paga"]["groups"]
        self.node_labels = self._resolve_node_labels()
        if (
            self.colors is None or self.colors == self.groups_key
        ) and self.groups_key is not None:
            self.colors = _paga_group_colors(self.adata, self.groups_key)
        self.nx_g_solid = nx.Graph(self.adjacency_solid)
        self.nx_g_dashed = (
            None if self.adjacency_dashed is None else nx.Graph(self.adjacency_dashed)
        )
        self.pos_array = (
            _read_gdf_pos(self.pos) if isinstance(self.pos, Path | str) else self.pos
        )
        self.pos_dict = {n: [p[0], p[1]] for n, p in enumerate(self.pos_array)}
        self.colors = self._resolve_colors()
        self._restrict_to_largest_component()

    def draw(self) -> PathCollection:
        self._draw_edges()
        if self.export_to_gexf:
            self._export_gexf()
        self.ax.set_frame_on(self.frameon)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        if isinstance(self.colors[0], Mapping):
            return self._draw_pie_nodes()
        return self._draw_scatter_nodes()

    def _resolve_node_labels(self):
        if (
            self.labels is not None
            and isinstance(self.labels, str)
            and self.labels != self.groups_key
        ):
            msg = (
                "Provide a list of group labels for the PAGA groups "
                f"{self.groups_key}, not {self.labels}."
            )
            raise ValueError(msg)
        if self.labels is None:
            return self.adata.obs[self.groups_key].cat.categories
        return self.labels

    def _resolve_colors(self):
        """Convert the color specification into one color or value per node."""
        colors = self.colors
        n_nodes = len(self.node_labels)
        # uniform color
        if isinstance(colors, str) and is_color_like(colors):
            colors = [colors for c in range(n_nodes)]
        # color degree of the graph
        if isinstance(colors, str) and colors.startswith("degree"):
            colors = self._degree_colors(colors)
        if isinstance(colors, str):
            colors = self._colors_from_key(colors)
        if len(colors) != n_nodes:
            msg = (
                f"Expected `colors` to be of length `{n_nodes}`, found `{len(colors)}`."
            )
            raise ValueError(msg)
        return colors

    def _degree_colors(self, kind: str) -> np.ndarray:
        # see also tools.paga.paga_degrees
        match kind:
            case "degree_dashed":
                if self.nx_g_dashed is None:
                    msg = (
                        "Coloring by `'degree_dashed'` requires passing `dashed_edges`."
                    )
                    raise NotImplementedError(msg)
                graph = self.nx_g_dashed
            case "degree_solid":
                graph = self.nx_g_solid
            case _:
                msg = '`degree` either "degree_dashed" or "degree_solid".'
                raise ValueError(msg)
        degrees = [d for _, d in graph.degree(weight="weight")]
        return (np.array(degrees) - np.min(degrees)) / (
            np.max(degrees) - np.min(degrees)
        )

    def _colors_from_key(self, key: str):
        """Aggregate a gene or `obs` annotation per PAGA group."""
        adata, groups_key = self.adata, self.groups_key
        cats = adata.obs[groups_key].cat.categories
        # plot gene expression
        var_names = adata.var_names if adata.raw is None else adata.raw.var_names
        if key in var_names:
            x_color = []
            for cat in cats:
                subset = (cat == adata.obs[groups_key]).to_numpy()
                if adata.raw is not None and self.use_raw:
                    adata_gene = adata.raw[:, key]
                else:
                    adata_gene = adata[:, key]
                x_color.append(np.mean(adata_gene.X[subset]))
            return x_color
        if key not in adata.obs:
            return key
        # plot continuous annotation
        if not isinstance(adata.obs[key].dtype, CategoricalDtype):
            x_color = []
            for cat in cats:
                subset = (cat == adata.obs[groups_key]).to_numpy()
                x_color.append(adata.obs.loc[subset, key].mean())
            return x_color
        # plot categorical annotation
        _asso_names, asso_matrix = _sc_utils.compute_association_matrix_of_groups(
            adata,
            prediction=groups_key,
            reference=key,
            normalization="reference" if self.normalize_to_color else "prediction",
        )
        _utils.add_colors_for_categorical_sample_annotation(adata, key)
        return _sc_utils.get_associated_colors_of_groups(
            adata.uns[f"{key}_colors"], asso_matrix
        )

    def _restrict_to_largest_component(self) -> None:
        import networkx as nx

        n_components, labels = scipy.sparse.csgraph.connected_components(
            self.adjacency_solid
        )
        if n_components <= 1:
            return
        if not self.single_component:
            logg.debug(
                "Graph has more than a single connected component. "
                "To restrict to this component, pass `single_component=True`."
            )
            return
        component_sizes = np.bincount(labels)
        largest_component = np.where(component_sizes == component_sizes.max())[0][0]
        mask = labels == largest_component
        self.adjacency_solid = self.adjacency_solid.tocsr()[mask, :]
        self.adjacency_solid = self.adjacency_solid.tocsc()[:, mask]
        self.colors = np.array(self.colors)[mask]
        self.node_labels = np.array(self.node_labels)[mask]
        cats_dropped = self.adata.obs[self.groups_key].cat.categories[~mask].tolist()
        logg.info(
            "Restricting graph to largest connected component by dropping categories\n"
            f"{cats_dropped}"
        )
        self.nx_g_solid = nx.Graph(self.adjacency_solid)
        if self.adjacency_dashed is not None:
            msg = "`single_component` only if `dashed_edges` is `None`."
            raise ValueError(msg)

    def _edge_widths(self, graph: nx.Graph, *, clip_min: bool = True) -> np.ndarray:
        base_edge_width = self.edge_width_scale * 5 * rcParams["lines.linewidth"]
        widths = [x[-1]["weight"] for x in graph.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        min_edge_width = self.min_edge_width if clip_min else None
        if min_edge_width is not None or self.max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, self.max_edge_width)
        return widths

    def _draw_edges(self) -> None:
        import networkx as nx

        draw = partial(nx.draw_networkx_edges, pos=self.pos_dict, ax=self.ax)
        if self.nx_g_dashed is not None:
            widths = self._edge_widths(self.nx_g_dashed, clip_min=False)
            draw(
                self.nx_g_dashed,
                width=widths,
                edge_color="grey",
                style="dashed",
                alpha=0.5,
            )
        # draw solid edges
        if self.transitions is None:
            widths = self._edge_widths(self.nx_g_solid)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                draw(self.nx_g_solid, width=widths, edge_color="black")
            return
        # draw directed edges
        adjacency_transitions = self.adata.uns["paga"][self.transitions].copy()
        threshold = 0.01 if self.threshold is None else self.threshold
        adjacency_transitions.data[adjacency_transitions.data < threshold] = 0
        adjacency_transitions.eliminate_zeros()
        g_dir = nx.DiGraph(adjacency_transitions.T)
        widths = self._edge_widths(g_dir)
        draw(g_dir, width=widths, edge_color="black", arrowsize=self.arrowsize)

    def _export_gexf(self) -> None:
        import networkx as nx

        if isinstance(self.colors[0], tuple):
            from matplotlib.colors import rgb2hex

            self.colors = [rgb2hex(c) for c in self.colors]
        for count, _n in enumerate(self.nx_g_solid.nodes()):
            node = self.nx_g_solid.node[count]
            node["label"] = str(self.node_labels[count])
            node["color"] = str(self.colors[count])
            node["viz"] = dict(
                position=dict(
                    x=1000 * self.pos_dict[count][0],
                    y=1000 * self.pos_dict[count][1],
                    z=0,
                )
            )
        filename = settings.writedir / "paga_graph.gexf"
        logg.warning(f"exporting to {filename}")
        settings.writedir.mkdir(parents=True, exist_ok=True)
        nx.write_gexf(self.nx_g_solid, settings.writedir / "paga_graph.gexf")

    def _node_sizes(self) -> np.ndarray:
        if self.groups_key is not None and f"{self.groups_key}_sizes" in self.adata.uns:
            groups_sizes = self.adata.uns[f"{self.groups_key}_sizes"]
        else:
            groups_sizes = np.ones(len(self.node_labels))
        base_scale_scatter = 2000
        n_nodes = self.adjacency_solid.shape[0]
        base_pie_size = (
            base_scale_scatter / (np.sqrt(n_nodes) + 10) * self.node_size_scale
        )
        median_group_size = np.median(groups_sizes)
        return base_pie_size * np.power(
            groups_sizes / median_group_size, self.node_size_power
        )

    def _text_kwds(self) -> dict[str, Any]:
        text_kwds = self.text_kwds
        if self.fontoutline is not None:
            text_kwds = dict(text_kwds)
            text_kwds["path_effects"] = [
                patheffects.withStroke(linewidth=self.fontoutline, foreground="w")
            ]
        return dict(
            verticalalignment="center",
            horizontalalignment="center",
            size=rcParams["legend.fontsize"]
            if self.fontsize is None
            else self.fontsize,
            fontweight=self.fontweight,
            **text_kwds,
        )

    def _draw_scatter_nodes(self) -> PathCollection:
        pos_array = self.pos_array
        sct = self.ax.scatter(
            pos_array[:, 0],
            pos_array[:, 1],
            c=self.colors[: len(pos_array)],
            edgecolors="face",
            s=self._node_sizes(),
            cmap=self.cmap,
        )
        text_kwds = self._text_kwds()
        for count, group in enumerate(self.node_labels):
            self.ax.text(pos_array[count, 0], pos_array[count, 1], group, **text_kwds)
        return sct

    def _draw_pie_nodes(self) -> PathCollection:
        sizes = self._node_sizes()
        text_kwds = self._text_kwds()
        xys = zip(self.pos_array[:, 0], self.pos_array[:, 1], strict=True)
        for ix, (xx, yy) in enumerate(xys):
            sct = _draw_pie_node(
                self.ax, xx, yy, self.colors[ix], ix=ix, size=sizes[ix]
            )
            if self.node_labels is not None:
                self.ax.text(xx, yy, self.node_labels[ix], **text_kwds)
        return sct


def _paga_group_colors(adata: AnnData, groups_key: str) -> list:
    if f"{groups_key}_colors" not in adata.uns or len(
        adata.obs[groups_key].cat.categories
    ) != len(adata.uns[f"{groups_key}_colors"]):
        _utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
    colors = adata.uns[f"{groups_key}_colors"]
    for iname, name in enumerate(adata.obs[groups_key].cat.categories):
        if name in settings.categories_to_ignore:
            colors[iname] = "grey"
    return colors


def _read_gdf_pos(pos: Path | str) -> np.ndarray:
    pos = Path(pos)
    if pos.suffix != ".gdf":
        msg = (
            "Currently only supporting reading positions from .gdf files. "
            "Consider generating them using, for instance, Gephi."
        )
        raise ValueError(msg)
    s = ""  # read the node definition from the file
    with pos.open() as f:
        f.readline()
        for line in f:
            if line.startswith("edgedef>"):
                break
            s += line
    from io import StringIO

    df = pd.read_csv(StringIO(s), header=-1)
    return df[[4, 5]].to_numpy()


def _draw_pie_node(
    ax: Axes, x: float, y: float, fracs_by_color, *, ix: int, size: float
):
    if not isinstance(fracs_by_color, Mapping):
        msg = (
            f"{fracs_by_color} is neither a dict of valid "
            "matplotlib colors nor a valid matplotlib color."
        )
        raise ValueError(msg)
    color_single = fracs_by_color.keys()
    fracs = [fracs_by_color[c] for c in color_single]
    total = sum(fracs)

    if total < 1:
        color_single = list(color_single)
        color_single.append("grey")
        fracs.append(1 - sum(fracs))
    elif not np.isclose(total, 1):
        msg = f"Expected fractions for node `{ix}` to be close to 1, found `{total}`."
        raise ValueError(msg)

    cumsum = np.cumsum(fracs)
    cumsum = cumsum / cumsum[-1]
    cumsum = [0, *cumsum.tolist()]

    for r1, r2, color in zip(cumsum[:-1], cumsum[1:], color_single, strict=True):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2, 20)
        xs = [0, *np.cos(angles).tolist()]
        ys = [0, *np.sin(angles).tolist()]

        xy = np.column_stack([xs, ys])
        s = np.abs(xy).max()

        sct = ax.scatter([x], [y], marker=xy, s=s**2 * size, color=color)
    return sct


def paga_path(  # noqa: PLR0913
    adata: AnnData,
    nodes: Sequence[str | int],
    keys: Sequence[str],
    *,
    use_raw: bool = True,
    annotations: Sequence[str] = ("dpt_pseudotime",),
    color_map: str | Colormap | None = None,
    color_maps_annotations: Mapping[str, str | Colormap] = frozendict(
        dict(dpt_pseudotime="Greys")
    ),
    palette_groups: Sequence[str] | None = None,
    n_avg: int = 1,
    groups_key: str | None = None,
    xlim: tuple[int | None, int | None] = (None, None),
    title: str | None = None,
    left_margin=None,
    ytick_fontsize: int | None = None,
    title_fontsize: int | None = None,
    show_node_names: bool = True,
    show_yticks: bool = True,
    show_colorbar: bool = True,
    legend_fontsize: float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight | None = None,
    normalize_to_zero_one: bool = False,
    as_heatmap: bool = True,
    return_data: bool = False,
    show: bool | None = None,
    ax: Axes | None = None,
    # deprecated
    save: bool | str | None = None,
) -> tuple[Axes, pd.DataFrame] | Axes | pd.DataFrame | None:
    r"""Gene expression and annotation changes along paths in the abstracted graph.

    Parameters
    ----------
    adata
        An annotated data matrix.
    nodes
        A path through nodes of the abstracted graph, that is, names or indices
        (within `.categories`) of groups that have been used to run PAGA.
    keys
        Either variables in `adata.var_names` or annotations in
        `adata.obs`. They are plotted using `color_map`.
    use_raw
        Use `adata.raw` for retrieving gene expressions if it has been set.
    annotations
        Plot these keys with `color_maps_annotations`. Need to be keys for
        `adata.obs`.
    color_map
        Matplotlib colormap.
    color_maps_annotations
        Color maps for plotting the annotations. Keys of the dictionary must
        appear in `annotations`.
    palette_groups
        Ususally, use the same `sc.pl.palettes...` as used for coloring the
        abstracted graph.
    n_avg
        Number of data points to include in computation of running average.
    groups_key
        Key of the grouping used to run PAGA. If `None`, defaults to
        `adata.uns['paga']['groups']`.
    as_heatmap
        Plot the timeseries as heatmap. If not plotting as heatmap,
        `annotations` have no effect.
    show_node_names
        Plot the node names on the nodes bar.
    show_colorbar
        Show the colorbar.
    show_yticks
        Show the y ticks.
    normalize_to_zero_one
        Shift and scale the running average to [0, 1] per gene.
    return_data
        Return the timeseries data in addition to the axes if `True`.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on \{`'.pdf'`, `'.png'`, `'.svg'`\}.
    ax
         A matplotlib axes object.

    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object, if `ax` is `None`, else `None`.
    If `return_data`, return the timeseries data in addition to an axes.

    """
    ax_was_none = ax is None
    path = _PagaPath(
        adata=adata,
        nodes=nodes,
        keys=keys,
        groups_key=groups_key,
        palette_groups=palette_groups,
        ax=ax,
        use_raw=use_raw,
        annotations=annotations,
        n_avg=n_avg,
        normalize_to_zero_one=normalize_to_zero_one,
        as_heatmap=as_heatmap,
        color_map=color_map,
        color_maps_annotations=color_maps_annotations,
        xlim=xlim,
        left_margin=left_margin,
        ytick_fontsize=ytick_fontsize,
        legend_fontsize=legend_fontsize,
        show_node_names=show_node_names,
        show_yticks=show_yticks,
        show_colorbar=show_colorbar,
    )
    path.collect_data()
    path.draw()
    if title is not None:
        path.ax.set_title(title, fontsize=title_fontsize)
    if show is None and not ax_was_none:
        show = False
    else:
        show = settings.autoshow if show is None else show
    _utils.savefig_or_show("paga_path", show=show, save=save)
    df = path.to_df() if return_data else None
    if not ax_was_none or show:
        return df
    return (path.ax, df) if return_data else path.ax


@dataclass(kw_only=True)
class _PagaPath:
    """Values of `keys` along a path through the PAGA graph, ordered by pseudotime."""

    adata: AnnData
    nodes: Sequence[str | int]
    keys: Sequence[str]
    groups_key: str | None
    palette_groups: Sequence[str] | None
    ax: Axes | None
    # data
    use_raw: bool
    annotations: Sequence[str]
    n_avg: int
    normalize_to_zero_one: bool
    # plotting
    as_heatmap: bool
    color_map: str | Colormap | None
    color_maps_annotations: Mapping[str, str | Colormap]
    xlim: tuple[int | None, int | None]
    left_margin: float | None
    ytick_fontsize: int | None
    legend_fontsize: float | _FontSize | None
    show_node_names: bool
    show_yticks: bool
    show_colorbar: bool

    groups_names: pd.Index = field(init=False)
    xs: list = field(init=False)
    """One array of values per key."""
    groups: list[int] = field(init=False)
    """Group index per data point."""
    x_tick_locs: list[int] = field(init=False)
    x_tick_labels: list[str | int] = field(init=False)
    anno_dict: dict[str, list] = field(init=False)
    """Values per data point for each annotation."""

    def __post_init__(self) -> None:
        if self.groups_key is None:
            if "groups" not in self.adata.uns["paga"]:
                msg = (
                    "Pass the key of the grouping with which you ran PAGA, "
                    "using the parameter `groups_key`."
                )
                raise KeyError(msg)
            self.groups_key = self.adata.uns["paga"]["groups"]
        self.groups_names = self.adata.obs[self.groups_key].cat.categories

        if "dpt_pseudotime" not in self.adata.obs.columns:
            msg = (
                "`pl.paga_path` requires computation of a pseudotime `tl.dpt` "
                "for ordering at single-cell resolution"
            )
            raise ValueError(msg)

        if self.palette_groups is None:
            _utils.add_colors_for_categorical_sample_annotation(
                self.adata, self.groups_key
            )
            self.palette_groups = self.adata.uns[f"{self.groups_key}_colors"]

        if self.ax is None:
            self.ax = plt.gca()

    def _moving_average(self, a):
        return _sc_utils.moving_average(a, self.n_avg)

    def collect_data(self) -> None:
        self.xs = []
        self.x_tick_locs = [0]
        self.x_tick_labels = []
        self.groups = []
        self.anno_dict = {anno: [] for anno in self.annotations}
        nodes_ints, nodes_strs = self._resolve_nodes()
        adata_x = (
            self.adata.raw
            if self.use_raw and self.adata.raw is not None
            else self.adata
        )

        for ikey, key in enumerate(self.keys):
            x = []
            for group, node_str in zip(nodes_ints, nodes_strs, strict=False):
                idcs = self._group_idcs(node_str=node_str, group=group)
                values = (
                    self.adata.obs[key].to_numpy()
                    if key in self.adata.obs
                    else adata_x[:, key].X
                )[idcs]
                x += (
                    values.toarray() if isinstance(values, CSBase) else values
                ).tolist()
                if ikey == 0:
                    self.groups += [group] * len(idcs)
                    self.x_tick_locs.append(len(x))
                    self._collect_annotations(idcs)
            if self.n_avg > 1:
                x = self._moving_average(x)
                if ikey == 0:
                    self._smooth_annotations()
            if self.normalize_to_zero_one:
                x -= np.min(x)
                x /= np.max(x)
            self.xs.append(x)
            if ikey == 0:
                self.x_tick_labels = [
                    self.groups_names[group]
                    if len(self.groups_names) > 0 and group not in self.groups_names
                    else group
                    for group in self.nodes
                ]

    def _resolve_nodes(self) -> tuple[Sequence[int], Sequence[str]]:
        if not isinstance(self.nodes[0], str):
            return self.nodes, [self.groups_names[node] for node in self.nodes]
        groups_names_set = set(self.groups_names)
        nodes_ints = []
        for node in self.nodes:
            if node not in groups_names_set:
                msg = (
                    f"Each node/group needs to be in {self.groups_names.tolist()} "
                    f"(groups_key={self.groups_key!r}) not {node!r}."
                )
                raise ValueError(msg)
            nodes_ints.append(self.groups_names.get_loc(node))
        return nodes_ints, self.nodes

    def _group_idcs(self, *, node_str: str, group: int) -> NDArray[np.intp]:
        """Get indices of cells in a group, sorted by pseudotime."""
        mask = (self.adata.obs[self.groups_key] == node_str).to_numpy()
        idcs = np.flatnonzero(mask)
        if len(idcs) == 0:
            msg = (
                "Did not find data points that match "
                f"`adata.obs[{self.groups_key!r}].values == {str(group)!r}`. "
                f"Check whether `adata.obs[{self.groups_key!r}]` "
                "actually contains what you expect."
            )
            raise ValueError(msg)
        idcs_group = np.argsort(self.adata.obs["dpt_pseudotime"].iloc[mask].to_numpy())
        return idcs[idcs_group]

    def _collect_annotations(self, idcs: NDArray[np.intp]) -> None:
        for anno in self.annotations:
            series = self.adata.obs[anno]
            if isinstance(series.dtype, CategoricalDtype):
                series = series.cat.codes
            self.anno_dict[anno] += series.iloc[idcs].to_list()

    def _smooth_annotations(self) -> None:
        for k in self.annotations:
            if not isinstance(self.anno_dict[k][0], str):
                self.anno_dict[k] = self._moving_average(self.anno_dict[k])

    @property
    def x(self) -> np.ndarray:
        return np.asarray(self.xs).squeeze()

    def draw(self) -> None:
        if not self.as_heatmap:
            self._draw_lines()
            return
        self._draw_heatmap()
        self._draw_groups_bar()
        self._draw_annotation_bars()

    def _draw_lines(self) -> None:
        for key, x in zip(self.keys, self.xs, strict=True):
            self.ax.plot(x[self.xlim[0] : self.xlim[1]], label=key)
        left_margin = 0.4 if self.left_margin is None else self.left_margin
        if len(self.keys) > 1:
            plt.legend(
                frameon=False,
                loc="center left",
                bbox_to_anchor=(-left_margin, 0.5),
                fontsize=self.legend_fontsize,
            )
        self.ax.set_xlabel(self.groups_key)
        plt.yticks([])
        if len(self.keys) == 1:
            plt.ylabel(f"{self.keys[0]} (a.u.)")

    def _draw_heatmap(self) -> None:
        ax, x = self.ax, self.x
        img = ax.imshow(x, aspect="auto", interpolation="nearest", cmap=self.color_map)
        if self.show_yticks:
            ax.set_yticks(range(len(x)))
            ax.set_yticklabels(self.keys, fontsize=self.ytick_fontsize)
        else:
            ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.tick_params(axis="both", which="both", length=0)
        ax.grid(visible=False)
        if self.show_colorbar:
            plt.colorbar(img, ax=ax)
        left_margin = 0.2 if self.left_margin is None else self.left_margin
        plt.subplots_adjust(left=left_margin)

    def _draw_bar(self, values, *, ianno: int, label: str, cmap) -> Axes:
        """Draw a bar below the heatmap (`ianno=-1` for the groups bar)."""
        ax_bounds = self.ax.get_position().bounds
        n_keys = len(self.keys)
        y_shift = ax_bounds[3] / n_keys if ianno <= 0 else ax_bounds[3] / n_keys / 2
        bar_axis = plt.axes((
            ax_bounds[0],
            ax_bounds[1] - (ianno + 2) * y_shift,
            ax_bounds[2],
            y_shift,
        ))
        bar_axis.imshow(
            np.array(values)[None, :],
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
        )
        if self.show_yticks:
            bar_axis.set_yticks(range(3), [label, "", ""], fontsize=self.ytick_fontsize)
        else:
            bar_axis.set_yticks([])
        bar_axis.set_frame_on(False)
        bar_axis.set_xticks([])
        bar_axis.grid(visible=False)
        return bar_axis

    def _draw_groups_bar(self) -> None:
        import matplotlib.colors

        groups = np.array(self.groups)
        groups_axis = self._draw_bar(
            groups,
            ianno=-1,
            label=self.groups_key,
            cmap=matplotlib.colors.ListedColormap(
                # the following line doesn't work because of normalization
                # adata.uns['paga_groups_colors'])
                self.palette_groups[int(np.min(groups)) : int(np.max(groups)) + 1]
            ),
        )
        if self.show_node_names:
            ypos = (groups_axis.get_ylim()[1] + groups_axis.get_ylim()[0]) / 2
            x_tick_locs = _sc_utils.moving_average(self.x_tick_locs, n=2)
            for loc, label in zip(x_tick_locs, self.x_tick_labels, strict=True):
                font = dict(horizontalalignment="center", verticalalignment="center")
                groups_axis.text(loc, ypos, label, fontdict=font)
        groups_axis.tick_params(axis="both", which="both", length=0)

    def _draw_annotation_bars(self) -> None:
        for ianno, anno in enumerate(self.annotations):
            if anno in self.color_maps_annotations:
                color_map_anno = self.color_maps_annotations[anno]
            elif isinstance(self.adata.obs[anno].dtype, CategoricalDtype):
                color_map_anno = "Vega10"
            else:
                color_map_anno = "Greys"
            anno_axis = self._draw_bar(
                self.anno_dict[anno], ianno=ianno, label=anno, cmap=color_map_anno
            )
            if self.show_yticks:
                anno_axis.tick_params(axis="both", which="both", length=0)

    def to_df(self) -> pd.DataFrame:
        df = pd.DataFrame(data=self.x.T, columns=self.keys)
        # groups is without moving average, yet
        df["groups"] = self._moving_average(self.groups)
        if "dpt_pseudotime" in self.anno_dict:
            df["distance"] = self.anno_dict["dpt_pseudotime"].T
        return df


@deprecated(Deprecation("1.13.0"))
def paga_adjacency(
    adata: AnnData,
    *,
    adjacency: str = "connectivities",
    adjacency_tree: str = "connectivities_tree",
    as_heatmap: bool = True,
    color_map: str | Colormap | None = None,
    show: bool | None = None,
    # deprecated
    save: bool | str | None = None,
) -> None:
    """Plot connectivity of paga groups."""
    connectivity = adata.uns[adjacency].toarray()
    connectivity_select = adata.uns[adjacency_tree]
    if as_heatmap:
        matrix(connectivity, color_map=color_map, show=False)
        for i in range(connectivity_select.shape[0]):
            neighbors = connectivity_select[i].nonzero()[1]
            plt.scatter([i for j in neighbors], neighbors, color="black", s=1)
    # as a stripplot
    else:
        plt.figure()
        for i, cs in enumerate(connectivity):
            x = [i for j, d in enumerate(cs) if i != j]
            y = [c for j, c in enumerate(cs) if i != j]
            plt.scatter(x, y, color="gray", s=1)
            neighbors = connectivity_select[i].nonzero()[1]
            plt.scatter([i for j in neighbors], cs[neighbors], color="black", s=1)
    _utils.savefig_or_show("paga_connectivity", show=show, save=save)
