from __future__ import annotations

import warnings
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scipy
from matplotlib import patheffects, rcParams, ticker
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from pandas.api.types import CategoricalDtype
from sklearn.utils import check_random_state

from scanpy.tools._draw_graph import coerce_fa2_layout, fa2_positions

from ... import _utils as _sc_utils
from ... import logging as logg
from ..._compat import CSBase, old_positionals
from ..._settings import settings
from .. import _utils
from .._utils import matrix

if TYPE_CHECKING:
    from typing import Any, Literal

    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.colors import Colormap

    from ..._compat import SpBase
    from ..._utils.random import _LegacyRandom
    from ...tools._draw_graph import _Layout as _LayoutWithoutEqTree
    from .._utils import _FontSize, _FontWeight, _LegendLoc

    _Layout = _LayoutWithoutEqTree | Literal["eq_tree"]


@old_positionals(
    "edges",
    "color",
    "alpha",
    "groups",
    "components",
    "projection",
    "legend_loc",
    "legend_fontsize",
    "legend_fontweight",
    "legend_fontoutline",
    "color_map",
    "palette",
    "frameon",
    "size",
    "title",
    "right_margin",
    "left_margin",
    "show",
    "save",
    "title_graph",
    "groups_graph",
)
def paga_compare(  # noqa: PLR0912, PLR0913
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

    """
    axs, _, _, _ = _utils.setup_axes(panels=[0, 1], right_margin=right_margin)
    if color is None:
        color = adata.uns["paga"]["groups"]
    suptitle = None  # common title for entire figure
    if title_graph is None:
        suptitle = color if title is None else title
        title, title_graph = "", ""
    if basis is None:
        if "X_draw_graph_fa" in adata.obsm:
            basis = "draw_graph_fa"
        elif "X_umap" in adata.obsm:
            basis = "umap"
        elif "X_tsne" in adata.obsm:
            basis = "tsne"
        elif "X_draw_graph_fr" in adata.obsm:
            basis = "draw_graph_fr"
        else:
            basis = "umap"

    from .scatterplots import _components_to_dimensions, _get_basis, embedding

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
        if color == adata.uns["paga"]["groups"]:
            # TODO: Use dimensions here
            _basis = _get_basis(adata, basis)
            dims = _components_to_dimensions(
                components=components, dimensions=None, total_dims=_basis.shape[1]
            )[0]
            coords = _basis[:, dims]
            pos = (
                pd.DataFrame(coords, columns=["x", "y"], index=adata.obs_names)
                .groupby(adata.obs[color], observed=True)
                .median()
                .sort_index()
            ).to_numpy()
        else:
            pos = adata.uns["paga"]["pos"]
    xlim, ylim = axs[0].get_xlim(), axs[0].get_ylim()
    axs[1].set_xlim(xlim)
    axs[1].set_ylim(ylim)
    if "labels" in paga_graph_params:
        labels = paga_graph_params.pop("labels")
    else:
        labels = groups_graph
    if legend_fontsize is not None:
        paga_graph_params["fontsize"] = legend_fontsize
    if legend_fontweight is not None:
        paga_graph_params["fontweight"] = legend_fontweight
    if legend_fontoutline is not None:
        paga_graph_params["fontoutline"] = legend_fontoutline
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


def _compute_pos(  # noqa: PLR0912
    adjacency_solid: SpBase | np.ndarray,
    *,
    layout: _Layout | None = None,
    random_state: _LegacyRandom = 0,
    init_pos: np.ndarray | None = None,
    adj_tree=None,
    root: int = 0,
    layout_kwds: Mapping[str, Any] = MappingProxyType({}),
):
    import random

    import networkx as nx

    random_state = check_random_state(random_state)

    nx_g_solid = nx.Graph(adjacency_solid)
    if layout is None:
        layout = "fr"
    layout = coerce_fa2_layout(layout)
    if layout == "fa":
        # np.random.seed(random_state)
        if init_pos is None:
            init_coords = random_state.random_sample((adjacency_solid.shape[0], 2))
        else:
            init_coords = init_pos.copy()
        pos_list = fa2_positions(adjacency_solid, init_coords, **layout_kwds)
        pos = {n: (x, -y) for n, (x, y) in enumerate(pos_list)}
    elif layout == "eq_tree":
        nx_g_tree = nx.Graph(adj_tree)
        pos = _utils.hierarchy_pos(nx_g_tree, root)
        if len(pos) < adjacency_solid.shape[0]:
            msg = (
                "This is a forest and not a single tree. "
                "Try another `layout`, e.g., {'fr'}."
            )
            raise ValueError(msg)
    else:
        # igraph layouts
        random.seed(random_state.bytes(8))
        g = _sc_utils.get_igraph_from_adjacency(adjacency_solid)
        if "rt" in layout:
            g_tree = _sc_utils.get_igraph_from_adjacency(adj_tree)
            pos_list = g_tree.layout(
                layout, root=root if isinstance(root, list) else [root]
            ).coords
        elif layout == "circle":
            pos_list = g.layout(layout).coords
        else:
            # I don't know why this is necessary
            # np.random.seed(random_state)
            if init_pos is None:
                init_coords = random_state.random_sample(
                    (adjacency_solid.shape[0], 2)
                ).tolist()
            else:
                init_pos = init_pos.copy()
                # this is a super-weird hack that is necessary as igraph’s
                # layout function seems to do some strange stuff here
                init_pos[:, 1] *= -1
                init_coords = init_pos.tolist()
            try:
                pos_list = g.layout(
                    layout, seed=init_coords, weights="weight", **layout_kwds
                ).coords
            except AttributeError:  # hack for empty graphs...
                pos_list = g.layout(layout, seed=init_coords, **layout_kwds).coords
        pos = {n: (x, -y) for n, (x, y) in enumerate(pos_list)}
    if len(pos) == 1:
        pos[0] = (0.5, 0.5)
    pos_array = np.array([pos[n] for count, n in enumerate(nx_g_solid)])
    return pos_array


@old_positionals(
    "threshold",
    "color",
    "layout",
    "layout_kwds",
    "init_pos",
    "root",
    "labels",
    "single_component",
    "solid_edges",
    "dashed_edges",
    "transitions",
    "fontsize",
    "fontweight",
    "fontoutline",
    "text_kwds",
    "node_size_scale",
    # 17 positionals are enough for backwards compat
)
def paga(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    *,
    threshold: float | None = None,
    color: str | Mapping[str | int, Mapping[Any, float]] | None = None,
    layout: _Layout | None = None,
    layout_kwds: Mapping[str, Any] = MappingProxyType({}),
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
    text_kwds: Mapping[str, Any] = MappingProxyType({}),
    node_size_scale: float = 1.0,
    node_size_power: float = 0.5,
    edge_width_scale: float = 1.0,
    min_edge_width: float | None = None,
    max_edge_width: float | None = None,
    arrowsize: int = 30,
    title: str | None = None,
    left_margin: float = 0.01,
    random_state: int | None = 0,
    pos: np.ndarray | Path | str | None = None,
    normalize_to_color: bool = False,
    cmap: str | Colormap | None = None,
    cax: Axes | None = None,
    colorbar=None,  # TODO: this seems to be unused
    cb_kwds: Mapping[str, Any] = MappingProxyType({}),
    frameon: bool | None = None,
    add_pos: bool = True,
    export_to_gexf: bool = False,
    use_raw: bool = True,
    colors=None,  # backwards compat
    groups=None,  # backwards compat
    plot: bool = True,
    show: bool | None = None,
    save: bool | str | None = None,
    ax: Axes | None = None,
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
        passing {`'degree_dashed'`, `'degree_solid'`}.

        Can be also used to visualize pie chart at each node in the following form:
        `{<group name or index>: {<color>: <fraction>, ...}, ...}`. If the fractions
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
    random_state
        For layouts with random initialization like `'fr'`, change this to use
        different intial states for the optimization. If `None`, the initial
        state is not reproducible.
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
        Infer the filetype if ending on \{`'.pdf'`, `'.png'`, `'.svg'`\}.
    ax
        A matplotlib axes object.

    Returns
    -------
    If `show==False`, one or more :class:`~matplotlib.axes.Axes` objects.
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata)

    You can increase node and edge sizes by specifying additional arguments.

    .. plot::
        :context: close-figs

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
    if groups is not None:  # backwards compat
        labels = groups
        logg.warning("`groups` is deprecated in `pl.paga`: use `labels` instead")
    if colors is None:
        colors = color

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

    if frameon is None:
        frameon = settings._frameon
    # labels is a list that contains no lists
    if is_flat(labels):
        labels = [labels for _ in range(len(colors))]

    if title is None and len(colors) > 1:
        title = list(colors)
    elif isinstance(title, str):
        title = [title] * len(colors)
    elif title is None:
        title = [None] * len(colors)

    if colorbar is None:
        var_names = adata.var_names if adata.raw is None else adata.raw.var_names
        colorbars = [
            (
                (c in adata.obs_keys() and adata.obs[c].dtype.name != "category")
                or (c in var_names)
            )
            for c in colors
        ]
    else:
        colorbars = [False for _ in colors]

    if isinstance(root, str):
        if root not in labels:
            msg = f"If `root` is a string, it needs to be one of {labels} not {root!r}."
            raise ValueError(msg)
        root = list(labels).index(root)
    if isinstance(root, Sequence) and root[0] in labels:
        root = [list(labels).index(r) for r in root]

    # define the adjacency matrices
    adjacency_solid = adata.uns["paga"][solid_edges].copy()
    adjacency_dashed = None
    if threshold is None:
        threshold = 0.01  # default threshold
    if threshold > 0:
        adjacency_solid.data[adjacency_solid.data < threshold] = 0
        adjacency_solid.eliminate_zeros()
    if dashed_edges is not None:
        adjacency_dashed = adata.uns["paga"][dashed_edges].copy()
        if threshold > 0:
            adjacency_dashed.data[adjacency_dashed.data < threshold] = 0
            adjacency_dashed.eliminate_zeros()

    # compute positions
    if pos is None:
        adj_tree = None
        if layout in {"rt", "rt_circular", "eq_tree"}:
            adj_tree = adata.uns["paga"]["connectivities_tree"]
        pos = _compute_pos(
            adjacency_solid,
            layout=layout,
            random_state=random_state,
            init_pos=init_pos,
            layout_kwds=layout_kwds,
            adj_tree=adj_tree,
            root=root,
        )

    if plot:
        axs, panel_pos, draw_region_width, figure_width = _utils.setup_axes(
            ax, panels=colors, colorbars=colorbars
        )

        if len(colors) == 1 and not isinstance(axs, list):
            axs = [axs]

        for icolor, c in enumerate(colors):
            if title[icolor] is not None:
                axs[icolor].set_title(title[icolor])
            sct = _paga_graph(
                adata,
                axs[icolor],
                colors=colors if isinstance(colors, Mapping) else c,
                solid_edges=solid_edges,
                dashed_edges=dashed_edges,
                transitions=transitions,
                threshold=threshold,
                adjacency_solid=adjacency_solid,
                adjacency_dashed=adjacency_dashed,
                root=root,
                labels=labels[icolor],
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
                colorbar=colorbars[icolor],
                cb_kwds=cb_kwds,
                use_raw=use_raw,
                title=title[icolor],
                export_to_gexf=export_to_gexf,
                single_component=single_component,
                arrowsize=arrowsize,
                pos=pos,
            )
            if colorbars[icolor]:
                if cax is None:
                    bottom = panel_pos[0][0]
                    height = panel_pos[1][0] - bottom
                    width = 0.006 * draw_region_width / len(colors)
                    left = panel_pos[2][2 * icolor + 1] + 0.2 * width
                    rectangle = [left, bottom, width, height]
                    fig = plt.gcf()
                    ax_cb = fig.add_axes(rectangle)
                else:
                    ax_cb = cax[icolor]

                _ = plt.colorbar(
                    sct,
                    format=ticker.FuncFormatter(_utils.ticks_formatter),
                    cax=ax_cb,
                )
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


def _paga_graph(  # noqa: PLR0912, PLR0913, PLR0915
    adata,
    ax,
    *,
    solid_edges=None,
    dashed_edges=None,
    adjacency_solid=None,
    adjacency_dashed=None,
    transitions=None,
    threshold=None,
    root=0,
    colors=None,
    labels=None,
    fontsize=None,
    fontweight=None,
    fontoutline=None,
    text_kwds: Mapping[str, Any] = MappingProxyType({}),
    node_size_scale=1.0,
    node_size_power=0.5,
    edge_width_scale=1.0,
    normalize_to_color="reference",
    title=None,
    pos=None,
    cmap=None,
    frameon=True,
    min_edge_width=None,
    max_edge_width=None,
    export_to_gexf=False,
    colorbar=None,
    use_raw=True,
    cb_kwds: Mapping[str, Any] = MappingProxyType({}),
    single_component=False,
    arrowsize=30,
):
    import networkx as nx

    node_labels = labels  # rename for clarity
    if (
        node_labels is not None
        and isinstance(node_labels, str)
        and node_labels != adata.uns["paga"]["groups"]
    ):
        msg = (
            "Provide a list of group labels for the PAGA groups "
            f"{adata.uns['paga']['groups']}, not {node_labels}."
        )
        raise ValueError(msg)
    groups_key = adata.uns["paga"]["groups"]
    if node_labels is None:
        node_labels = adata.obs[groups_key].cat.categories

    if (colors is None or colors == groups_key) and groups_key is not None:
        if groups_key + "_colors" not in adata.uns or len(
            adata.obs[groups_key].cat.categories
        ) != len(adata.uns[groups_key + "_colors"]):
            _utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        colors = adata.uns[groups_key + "_colors"]
        for iname, name in enumerate(adata.obs[groups_key].cat.categories):
            if name in settings.categories_to_ignore:
                colors[iname] = "grey"

    nx_g_solid = nx.Graph(adjacency_solid)
    if dashed_edges is not None:
        nx_g_dashed = nx.Graph(adjacency_dashed)

    # convert pos to array and dict
    if not isinstance(pos, Path | str):
        pos_array = pos
    else:
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
        pos_array = df[[4, 5]].values

    # convert to dictionary
    pos = {n: [p[0], p[1]] for n, p in enumerate(pos_array)}

    # uniform color
    if isinstance(colors, str) and is_color_like(colors):
        colors = [colors for c in range(len(node_labels))]

    # color degree of the graph
    if isinstance(colors, str) and colors.startswith("degree"):
        # see also tools.paga.paga_degrees
        if colors == "degree_dashed":
            colors = [d for _, d in nx_g_dashed.degree(weight="weight")]
        elif colors == "degree_solid":
            colors = [d for _, d in nx_g_solid.degree(weight="weight")]
        else:
            msg = '`degree` either "degree_dashed" or "degree_solid".'
            raise ValueError(msg)
        colors = (np.array(colors) - np.min(colors)) / (np.max(colors) - np.min(colors))

    # plot gene expression
    var_names = adata.var_names if adata.raw is None else adata.raw.var_names
    if isinstance(colors, str) and colors in var_names:
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for cat in cats:
            subset = (cat == adata.obs[groups_key]).values
            if adata.raw is not None and use_raw:
                adata_gene = adata.raw[:, colors]
            else:
                adata_gene = adata[:, colors]
            x_color.append(np.mean(adata_gene.X[subset]))
        colors = x_color

    # plot continuous annotation
    if (
        isinstance(colors, str)
        and colors in adata.obs
        and not isinstance(adata.obs[colors].dtype, CategoricalDtype)
    ):
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for cat in cats:
            subset = (cat == adata.obs[groups_key]).values
            x_color.append(adata.obs.loc[subset, colors].mean())
        colors = x_color

    # plot categorical annotation
    if (
        isinstance(colors, str)
        and colors in adata.obs
        and isinstance(adata.obs[colors].dtype, CategoricalDtype)
    ):
        asso_names, asso_matrix = _sc_utils.compute_association_matrix_of_groups(
            adata,
            prediction=groups_key,
            reference=colors,
            normalization="reference" if normalize_to_color else "prediction",
        )
        _utils.add_colors_for_categorical_sample_annotation(adata, colors)
        asso_colors = _sc_utils.get_associated_colors_of_groups(
            adata.uns[colors + "_colors"], asso_matrix
        )
        colors = asso_colors

    if len(colors) != len(node_labels):
        msg = (
            f"Expected `colors` to be of length `{len(node_labels)}`, "
            f"found `{len(colors)}`."
        )
        raise ValueError(msg)

    # count number of connected components
    n_components, labels = scipy.sparse.csgraph.connected_components(adjacency_solid)
    if n_components > 1 and not single_component:
        logg.debug(
            "Graph has more than a single connected component. "
            "To restrict to this component, pass `single_component=True`."
        )
    if n_components > 1 and single_component:
        component_sizes = np.bincount(labels)
        largest_component = np.where(component_sizes == component_sizes.max())[0][0]
        adjacency_solid = adjacency_solid.tocsr()[labels == largest_component, :]
        adjacency_solid = adjacency_solid.tocsc()[:, labels == largest_component]
        colors = np.array(colors)[labels == largest_component]
        node_labels = np.array(node_labels)[labels == largest_component]
        cats_dropped = (
            adata.obs[groups_key].cat.categories[labels != largest_component].tolist()
        )
        logg.info(
            "Restricting graph to largest connected component by dropping categories\n"
            f"{cats_dropped}"
        )
        nx_g_solid = nx.Graph(adjacency_solid)
        if dashed_edges is not None:
            msg = "`single_component` only if `dashed_edges` is `None`."
            raise ValueError(msg)

    # edge widths
    base_edge_width = edge_width_scale * 5 * rcParams["lines.linewidth"]

    # draw dashed edges
    if dashed_edges is not None:
        widths = [x[-1]["weight"] for x in nx_g_dashed.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if max_edge_width is not None:
            widths = np.clip(widths, None, max_edge_width)
        nx.draw_networkx_edges(
            nx_g_dashed,
            pos,
            ax=ax,
            width=widths,
            edge_color="grey",
            style="dashed",
            alpha=0.5,
        )

    # draw solid edges
    if transitions is None:
        widths = [x[-1]["weight"] for x in nx_g_solid.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nx.draw_networkx_edges(
                nx_g_solid, pos, ax=ax, width=widths, edge_color="black"
            )
    # draw directed edges
    else:
        adjacency_transitions = adata.uns["paga"][transitions].copy()
        if threshold is None:
            threshold = 0.01
        adjacency_transitions.data[adjacency_transitions.data < threshold] = 0
        adjacency_transitions.eliminate_zeros()
        g_dir = nx.DiGraph(adjacency_transitions.T)
        widths = [x[-1]["weight"] for x in g_dir.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
        nx.draw_networkx_edges(
            g_dir, pos, ax=ax, width=widths, edge_color="black", arrowsize=arrowsize
        )

    if export_to_gexf:
        if isinstance(colors[0], tuple):
            from matplotlib.colors import rgb2hex

            colors = [rgb2hex(c) for c in colors]
        for count, _n in enumerate(nx_g_solid.nodes()):
            nx_g_solid.node[count]["label"] = str(node_labels[count])
            nx_g_solid.node[count]["color"] = str(colors[count])
            nx_g_solid.node[count]["viz"] = dict(
                position=dict(
                    x=1000 * pos[count][0],
                    y=1000 * pos[count][1],
                    z=0,
                )
            )
        filename = settings.writedir / "paga_graph.gexf"
        logg.warning(f"exporting to {filename}")
        settings.writedir.mkdir(parents=True, exist_ok=True)
        nx.write_gexf(nx_g_solid, settings.writedir / "paga_graph.gexf")

    ax.set_frame_on(frameon)
    ax.set_xticks([])
    ax.set_yticks([])

    # groups sizes
    if groups_key is not None and groups_key + "_sizes" in adata.uns:
        groups_sizes = adata.uns[groups_key + "_sizes"]
    else:
        groups_sizes = np.ones(len(node_labels))
    base_scale_scatter = 2000
    base_pie_size = (
        base_scale_scatter / (np.sqrt(adjacency_solid.shape[0]) + 10) * node_size_scale
    )
    median_group_size = np.median(groups_sizes)
    groups_sizes = base_pie_size * np.power(
        groups_sizes / median_group_size, node_size_power
    )

    if fontsize is None:
        fontsize = rcParams["legend.fontsize"]
    if fontoutline is not None:
        text_kwds = dict(text_kwds)
        text_kwds["path_effects"] = [
            patheffects.withStroke(linewidth=fontoutline, foreground="w")
        ]
    # usual scatter plot
    if not isinstance(colors[0], Mapping):
        n_groups = len(pos_array)
        sct = ax.scatter(
            pos_array[:, 0],
            pos_array[:, 1],
            c=colors[:n_groups],
            edgecolors="face",
            s=groups_sizes,
            cmap=cmap,
        )
        for count, group in enumerate(node_labels):
            ax.text(
                pos_array[count, 0],
                pos_array[count, 1],
                group,
                verticalalignment="center",
                horizontalalignment="center",
                size=fontsize,
                fontweight=fontweight,
                **text_kwds,
            )
    # else pie chart plot
    else:
        for ix, (xx, yy) in enumerate(
            zip(pos_array[:, 0], pos_array[:, 1], strict=True)
        ):
            if not isinstance(colors[ix], Mapping):
                msg = (
                    f"{colors[ix]} is neither a dict of valid "
                    "matplotlib colors nor a valid matplotlib color."
                )
                raise ValueError(msg)
            color_single = colors[ix].keys()
            fracs = [colors[ix][c] for c in color_single]
            total = sum(fracs)

            if total < 1:
                color_single = list(color_single)
                color_single.append("grey")
                fracs.append(1 - sum(fracs))
            elif not np.isclose(total, 1):
                msg = (
                    f"Expected fractions for node `{ix}` to be "
                    f"close to 1, found `{total}`."
                )
                raise ValueError(msg)

            cumsum = np.cumsum(fracs)
            cumsum = cumsum / cumsum[-1]
            cumsum = [0, *cumsum.tolist()]

            for r1, r2, color in zip(
                cumsum[:-1], cumsum[1:], color_single, strict=True
            ):
                angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2, 20)
                x = [0, *np.cos(angles).tolist()]
                y = [0, *np.sin(angles).tolist()]

                xy = np.column_stack([x, y])
                s = np.abs(xy).max()

                sct = ax.scatter(
                    [xx], [yy], marker=xy, s=s**2 * groups_sizes[ix], color=color
                )

            if node_labels is not None:
                ax.text(
                    xx,
                    yy,
                    node_labels[ix],
                    verticalalignment="center",
                    horizontalalignment="center",
                    size=fontsize,
                    fontweight=fontweight,
                    **text_kwds,
                )

    return sct


@old_positionals(
    "use_raw",
    "annotations",
    "color_map",
    "color_maps_annotations",
    "palette_groups",
    "n_avg",
    "groups_key",
    "xlim",
    "title",
    "left_margin",
    "ytick_fontsize",
    "title_fontsize",
    "show_node_names",
    "show_yticks",
    "show_colorbar",
    "legend_fontsize",
    "legend_fontweight",
    "normalize_to_zero_one",
    "as_heatmap",
    "return_data",
    "show",
    "save",
    "ax",
)
def paga_path(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    nodes: Sequence[str | int],
    keys: Sequence[str],
    *,
    use_raw: bool = True,
    annotations: Sequence[str] = ("dpt_pseudotime",),
    color_map: str | Colormap | None = None,
    color_maps_annotations: Mapping[str, str | Colormap] = MappingProxyType(
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
    save: bool | str | None = None,
    ax: Axes | None = None,
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

    if groups_key is None:
        if "groups" not in adata.uns["paga"]:
            msg = (
                "Pass the key of the grouping with which you ran PAGA, "
                "using the parameter `groups_key`."
            )
            raise KeyError(msg)
        groups_key = adata.uns["paga"]["groups"]
    groups_names = adata.obs[groups_key].cat.categories

    if "dpt_pseudotime" not in adata.obs.columns:
        msg = (
            "`pl.paga_path` requires computation of a pseudotime `tl.dpt` "
            "for ordering at single-cell resolution"
        )
        raise ValueError(msg)

    if palette_groups is None:
        _utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        palette_groups = adata.uns[f"{groups_key}_colors"]

    def moving_average(a):
        return _sc_utils.moving_average(a, n_avg)

    ax = plt.gca() if ax is None else ax

    X = []
    x_tick_locs = [0]
    x_tick_labels = []
    groups = []
    anno_dict = {anno: [] for anno in annotations}
    if isinstance(nodes[0], str):
        nodes_ints = []
        groups_names_set = set(groups_names)
        for node in nodes:
            if node not in groups_names_set:
                msg = (
                    f"Each node/group needs to be in {groups_names.tolist()} "
                    f"({groups_key=!r}) not {node!r}."
                )
                raise ValueError(msg)
            nodes_ints.append(groups_names.get_loc(node))
        nodes_strs = nodes
    else:
        nodes_ints = nodes
        nodes_strs = [groups_names[node] for node in nodes]

    adata_X = adata
    if use_raw and adata.raw is not None:
        adata_X = adata.raw

    for ikey, key in enumerate(keys):
        x = []
        for igroup, group in enumerate(nodes_ints):
            idcs = np.arange(adata.n_obs)[
                adata.obs[groups_key].values == nodes_strs[igroup]
            ]
            if len(idcs) == 0:
                msg = (
                    "Did not find data points that match "
                    f"`adata.obs[{groups_key!r}].values == {str(group)!r}`. "
                    f"Check whether `adata.obs[{groups_key!r}]` "
                    "actually contains what you expect."
                )
                raise ValueError(msg)
            idcs_group = np.argsort(
                adata.obs["dpt_pseudotime"].values[
                    adata.obs[groups_key].values == nodes_strs[igroup]
                ]
            )
            idcs = idcs[idcs_group]
            values = (
                adata.obs[key].values if key in adata.obs_keys() else adata_X[:, key].X
            )[idcs]
            x += (values.toarray() if isinstance(values, CSBase) else values).tolist()
            if ikey == 0:
                groups += [group] * len(idcs)
                x_tick_locs.append(len(x))
                for anno in annotations:
                    series = adata.obs[anno]
                    if isinstance(series.dtype, CategoricalDtype):
                        series = series.cat.codes
                    anno_dict[anno] += list(series.values[idcs])
        if n_avg > 1:
            x = moving_average(x)
            if ikey == 0:
                for k in annotations:
                    if not isinstance(anno_dict[k][0], str):
                        anno_dict[k] = moving_average(anno_dict[k])
        if normalize_to_zero_one:
            x -= np.min(x)
            x /= np.max(x)
        X.append(x)
        if not as_heatmap:
            ax.plot(x[xlim[0] : xlim[1]], label=key)
        if ikey == 0:
            for group in nodes:
                if len(groups_names) > 0 and group not in groups_names:
                    label = groups_names[group]
                else:
                    label = group
                x_tick_labels.append(label)
    X = np.asarray(X).squeeze()
    if as_heatmap:
        img = ax.imshow(X, aspect="auto", interpolation="nearest", cmap=color_map)
        if show_yticks:
            ax.set_yticks(range(len(X)))
            ax.set_yticklabels(keys, fontsize=ytick_fontsize)
        else:
            ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.tick_params(axis="both", which="both", length=0)
        ax.grid(visible=False)
        if show_colorbar:
            plt.colorbar(img, ax=ax)
        left_margin = 0.2 if left_margin is None else left_margin
        plt.subplots_adjust(left=left_margin)
    else:
        left_margin = 0.4 if left_margin is None else left_margin
        if len(keys) > 1:
            plt.legend(
                frameon=False,
                loc="center left",
                bbox_to_anchor=(-left_margin, 0.5),
                fontsize=legend_fontsize,
            )
    xlabel = groups_key
    if not as_heatmap:
        ax.set_xlabel(xlabel)
        plt.yticks([])
        if len(keys) == 1:
            plt.ylabel(keys[0] + " (a.u.)")
    else:
        import matplotlib.colors

        # groups bar
        ax_bounds = ax.get_position().bounds
        groups_axis = plt.axes(
            (
                ax_bounds[0],
                ax_bounds[1] - ax_bounds[3] / len(keys),
                ax_bounds[2],
                ax_bounds[3] / len(keys),
            )
        )
        groups = np.array(groups)[None, :]
        groups_axis.imshow(
            groups,
            aspect="auto",
            interpolation="nearest",
            cmap=matplotlib.colors.ListedColormap(
                # the following line doesn't work because of normalization
                # adata.uns['paga_groups_colors'])
                palette_groups[np.min(groups).astype(int) :],
                N=int(np.max(groups) + 1 - np.min(groups)),
            ),
        )
        if show_yticks:
            groups_axis.set_yticklabels(["", xlabel, ""], fontsize=ytick_fontsize)
        else:
            groups_axis.set_yticks([])
        groups_axis.set_frame_on(False)
        if show_node_names:
            ypos = (groups_axis.get_ylim()[1] + groups_axis.get_ylim()[0]) / 2
            x_tick_locs = _sc_utils.moving_average(x_tick_locs, n=2)
            for loc, label in zip(x_tick_locs, x_tick_labels, strict=True):
                font = dict(horizontalalignment="center", verticalalignment="center")
                groups_axis.text(loc, ypos, label, fontdict=font)
        groups_axis.set_xticks([])
        groups_axis.grid(visible=False)
        groups_axis.tick_params(axis="both", which="both", length=0)
        # further annotations
        y_shift = ax_bounds[3] / len(keys)
        for ianno, anno in enumerate(annotations):
            if ianno > 0:
                y_shift = ax_bounds[3] / len(keys) / 2
            anno_axis = plt.axes(
                (
                    ax_bounds[0],
                    ax_bounds[1] - (ianno + 2) * y_shift,
                    ax_bounds[2],
                    y_shift,
                )
            )
            arr = np.array(anno_dict[anno])[None, :]
            if anno not in color_maps_annotations:
                color_map_anno = (
                    "Vega10"
                    if isinstance(adata.obs[anno].dtype, CategoricalDtype)
                    else "Greys"
                )
            else:
                color_map_anno = color_maps_annotations[anno]
            img = anno_axis.imshow(
                arr,
                aspect="auto",
                interpolation="nearest",
                cmap=color_map_anno,
            )
            if show_yticks:
                anno_axis.set_yticklabels(["", anno, ""], fontsize=ytick_fontsize)
                anno_axis.tick_params(axis="both", which="both", length=0)
            else:
                anno_axis.set_yticks([])
            anno_axis.set_frame_on(False)
            anno_axis.set_xticks([])
            anno_axis.grid(visible=False)
    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    if show is None and not ax_was_none:
        show = False
    else:
        show = settings.autoshow if show is None else show
    _utils.savefig_or_show("paga_path", show=show, save=save)
    if return_data:
        df = pd.DataFrame(data=X.T, columns=keys)
        df["groups"] = moving_average(groups)  # groups is without moving average, yet
        if "dpt_pseudotime" in anno_dict:
            df["distance"] = anno_dict["dpt_pseudotime"].T
    if not ax_was_none or show:
        return df if return_data else None
    return (ax, df) if return_data else ax


def paga_adjacency(
    adata: AnnData,
    *,
    adjacency: str = "connectivities",
    adjacency_tree: str = "connectivities_tree",
    as_heatmap: bool = True,
    color_map: str | Colormap | None = None,
    show: bool | None = None,
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
