import os
import numpy as np
import pandas as pd
import scipy
import warnings
from collections import Iterable
from typing import Optional, Union, List

import networkx as nx
from pandas.api.types import is_categorical_dtype
from matplotlib import pyplot as pl, rcParams, ticker
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like

from .. import _utils as utils
from .._utils import matrix
from ... import utils as sc_utils, logging as logg
from ..._settings import settings


def paga_compare(
        adata,
        basis=None,
        edges=False,
        color=None,
        alpha=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='on data',
        legend_fontsize=None,
        legend_fontweight='bold',
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
        **paga_graph_params):
    """Scatter and PAGA graph side-by-side.

    Consists in a scatter plot and the abstracted graph. See
    :func:`~scanpy.api.pl.paga` for all related parameters.

    See :func:`~scanpy.api.pl.paga_path` for visualizing gene changes along paths
    through the abstracted graph.

    Additional parameters are as follows.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    kwds_scatter : `dict`
        Keywords for :func:`~scanpy.api.pl.scatter`.
    kwds_paga : `dict`
        Keywords for :func:`~scanpy.api.pl.paga`.

    Returns
    -------
    A list of `matplotlib.axes.Axes` if `show` is `False`.
    """
    axs, _, _, _ = utils.setup_axes(panels=[0, 1],
                                    right_margin=right_margin)
    if color is None:
        color = adata.uns['paga']['groups']
    suptitle = None  # common title for entire figure
    if title_graph is None:
        suptitle = color if title is None else title
        title, title_graph = '', ''
    if basis is None:
        if 'X_draw_graph_fa' in adata.obsm.keys():
            basis = 'draw_graph_fa'
        elif 'X_umap' in adata.obsm.keys():
            basis = 'umap'
        elif 'X_tsne' in adata.obsm.keys():
            basis = 'tsne'
        elif 'X_draw_graph_fr' in adata.obsm.keys():
            basis = 'draw_graph_fr'
        else:
            basis = 'umap'
    from .scatterplots import plot_scatter
    plot_scatter(
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
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        size=size,
        title=title,
        show=False,
        save=False)
    if 'pos' not in paga_graph_params:
        if color == adata.uns['paga']['groups']:
            paga_graph_params['pos'] = utils._tmp_cluster_pos
        else:
            paga_graph_params['pos'] = adata.uns['paga']['pos']
    xlim, ylim = axs[0].get_xlim(), axs[0].get_ylim()
    axs[1].set_xlim(xlim)
    axs[1].set_ylim(ylim)
    if 'labels' in paga_graph_params:
        labels = paga_graph_params.pop('labels')
    else:
        labels = groups_graph
    paga(
        adata,
        ax=axs[1],
        show=False,
        save=False,
        title=title_graph,
        labels=labels,
        colors=color,
        frameon=frameon,
        **paga_graph_params)
    if suptitle is not None: pl.suptitle(suptitle)
    utils.savefig_or_show('paga_compare', show=show, save=save)
    if show == False: return axs


def _compute_pos(adjacency_solid, layout=None, random_state=0, init_pos=None, adj_tree=None, root=0, layout_kwds={}):
    nx_g_solid = nx.Graph(adjacency_solid)
    if layout is None:
        layout = 'fr'
    if layout == 'fa':
        try:
            from fa2 import ForceAtlas2
        except:
            logg.warn('Package \'fa2\' is not installed, falling back to layout \'fr\'.'
                      'To use the faster and better ForceAtlas2 layout, '
                      'install package \'fa2\' (`pip install fa2`).')
            layout = 'fr'
    if layout == 'fa':
        np.random.seed(random_state)
        if init_pos is None:
            init_coords = np.random.random((adjacency_solid.shape[0], 2))
        else:
            init_coords = init_pos.copy()
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
        if 'maxiter' in layout_kwds:
            iterations = layout_kwds['maxiter']
        elif 'iterations' in layout_kwds:
            iterations = layout_kwds['iterations']
        else:
            iterations = 500
        pos_list = forceatlas2.forceatlas2(
            adjacency_solid, pos=init_coords, iterations=iterations)
        pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
    elif layout == 'eq_tree':
        nx_g_tree = nx.Graph(adj_tree)
        pos = utils.hierarchy_pos(nx_g_tree, root)
        if len(pos) < adjacency_solid.shape[0]:
            raise ValueError('This is a forest and not a single tree. '
                             'Try another `layout`, e.g., {\'fr\'}.')
    else:
        # igraph layouts
        from ... import utils as sc_utils
        g = sc_utils.get_igraph_from_adjacency(adjacency_solid)
        if 'rt' in layout:
            g_tree = sc_utils.get_igraph_from_adjacency(adj_tree)
            pos_list = g_tree.layout(
                layout, root=root if isinstance(root, list) else [root]).coords
        elif layout == 'circle':
            pos_list = g.layout(layout).coords
        else:
            # I don't know why this is necessary
            np.random.seed(random_state)
            if init_pos is None:
                init_coords = np.random.random((adjacency_solid.shape[0], 2)).tolist()
            else:
                init_pos = init_pos.copy()
                # this is a super-weird hack that is necessary as igraphs layout function
                # seems to do some strange stuff, here
                init_pos[:, 1] *= -1
                init_coords = init_pos.tolist()
            try:
                pos_list = g.layout(
                    layout, seed=init_coords,
                    weights='weight', **layout_kwds).coords
            except:  # hack for excepting attribute error for empty graphs...
                pos_list = g.layout(
                    layout, seed=init_coords,
                    **layout_kwds).coords
        pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
    if len(pos) == 1: pos[0] = (0.5, 0.5)
    pos_array = np.array([pos[n] for count, n in enumerate(nx_g_solid)])
    return pos_array


def paga(
    adata,
    threshold=None,
    color=None,
    layout=None,
    layout_kwds={},
    init_pos=None,
    root=0,
    labels=None,
    single_component=False,
    solid_edges='connectivities',
    dashed_edges=None,
    transitions=None,
    fontsize=None,
    fontweight='bold',
    text_kwds={},
    node_size_scale=1,
    node_size_power=0.5,
    edge_width_scale=1,
    min_edge_width=None,
    max_edge_width=None,
    arrowsize=30,
    title=None,
    left_margin=0.01,
    random_state=0,
    pos=None,
    normalize_to_color=False,
    cmap=None,
    cax=None,
    colorbar=None,
    cb_kwds={},
    frameon=None,
    add_pos=True,
    export_to_gexf=False,
    use_raw=True,
    colors=None,   # backwards compat
    groups=None,  # backwards compat
    plot=True,
    show=None,
    save=None,
    ax=None,
) -> Union[Axes, List[Axes], None]:
    """Plot the PAGA graph through thresholding low-connectivity edges.

    Compute a coarse-grained layout of the data. Reuse this by passing
    `init_pos='paga'` to :func:`~scanpy.api.tl.umap` or
    :func:`~scanpy.api.tl.draw_graph` and obtain embeddings with more meaningful
    global topology [Wolf19]_.

    This uses ForceAtlas2 or igraph's layout algorithms for most layouts [Csardi06]_.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    threshold : `float` or `None`, optional (default: 0.01)
        Do not draw edges for weights below this threshold. Set to 0 if you want
        all edges. Discarding low-connectivity edges helps in getting a much
        clearer picture of the graph.
    color : gene name or obs. annotation, optional (default: `None`)
        The node colors. Also plots the degree of the abstracted graph when
        passing {'degree_dashed', 'degree_solid'}.
    labels : `None`, `str`, `list`, `dict`, optional (default: `None`)
        The node labels. If `None`, this defaults to the group labels stored in
        the categorical for which :func:`~scanpy.api.tl.paga` has been computed.
    pos : `np.ndarray`, filename of `.gdf` file,  optional (default: `None`)
        Two-column array-like storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    layout : {'fa', 'fr', 'rt', 'rt_circular', 'eq_tree', ...}, optional (default: 'fr')
        Plotting layout that computes positions. 'fa' stands for ForceAtlas2, 'fr' stands for
        Fruchterman-Reingold, 'rt' stands for Reingold Tilford. 'eq_tree' stands
        for 'eqally spaced tree'. All but 'fa' and 'eq_tree' are igraph
        layouts. All other igraph layouts are also permitted. See also parameter
        `pos` and :func:`~scanpy.api.tl.draw_graph`.
    init_pos : `np.ndarray`, optional (default: `None`)
        Two-column array storing the x and y coordinates for initializing the
        layout.
    random_state : `int` or `None`, optional (default: 0)
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, the initial
        state is not reproducible.
    root : `int`, `str` or list of `int`, optional (default: 0)
        If choosing a tree layout, this is the index of the root node or a list
        of root node indices. If this is a non-empty vector then the supplied
        node IDs are used as the roots of the trees (or a single tree if the
        graph is connected). If this is `None` or an empty list, the root
        vertices are automatically calculated based on topological sorting.
    transitions : `str` or `None`, optional (default: `None`)
        Key for `.uns['paga']` that specifies the matrix that - for instance
        `'transistions_confidence'` - that specifies the matrix that stores the
        arrows.
    solid_edges : `str`, optional (default: 'paga_connectivities')
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges : `str` or `None`, optional (default: `None`)
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    single_component : `bool`, optional (default: `False`)
        Restrict to largest connected component.
    fontsize : `int` (default: `None`)
        Font size for node labels.
    text_kwds : keywords for `matplotlib.text`
        See `here
        <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.text.html#matplotlib.axes.Axes.text>`_.
    node_size_scale : `float` (default: 1.0)
        Increase or decrease the size of the nodes.
    node_size_power : `float` (default: 0.5)
        The power with which groups sizes influence the radius of the nodes.
    edge_width_scale : `float`, optional (default: 5)
        Edge with scale in units of `rcParams['lines.linewidth']`.
    min_edge_width : `float`, optional (default: `None`)
        Min width of solid edges.
    max_edge_width : `float`, optional (default: `None`)
        Max width of solid and dashed edges.
    arrowsize : `int`, optional (default: 30)
       For directed graphs, choose the size of the arrow head head's length and
       width. See :py:class: `matplotlib.patches.FancyArrowPatch` for attribute
       `mutation_scale` for more info.
    export_to_gexf : `bool`, optional (default: `None`)
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    normalize_to_color : `bool`, optional (default: `False`)
        Whether to normalize categorical plots to `color` or the underlying
        grouping.
    cmap : color map
        The color map.
    cax : :class:`~matplotlib.axes.Axes`
        A matplotlib axes object for a potential colorbar.
    cb_kwds : colorbar keywords
        See `here
        <https://matplotlib.org/api/colorbar_api.html#matplotlib.colorbar.ColorbarBase>`__,
        for instance, `ticks`.
    add_pos : `bool`, optional (default: `True`)
        Add the positions to `adata.uns['paga']`.
    title : `str`, optional (default: `None`)
        Provide a title.
    frameon : `bool`, optional (default: `None`)
        Draw a frame around the PAGA graph.
    hide : `bool`, optional (default: `False`)
        Do not create a plot.
    plot : `bool`, optional (default: `True`)
        If `False`, do not create the figure, simply compute the layout.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \\{'.pdf', '.png', '.svg'\\}.
    ax : :class:`~matplotlib.axes.Axes`
        A matplotlib axes object.

    Returns
    -------
    If `show==False`, one or more :class:`~matplotlib.axes.Axes` objects.
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.

    Notes
    -----
    When initializing the positions, note that - for some reason - igraph
    mirrors coordinates along the x axis... that is, you should increase the
    `maxiter` parameter by 1 if the layout is flipped.

    See also
    --------
    tl.paga
    pl.paga_compare
    pl.paga_path
    """
    if groups is not None:  # backwards compat
        labels = groups
        logg.warn('`groups` is deprecated in `pl.paga`: use `labels` instead')
    if colors is None:
        colors = color
    # colors is a list that contains no lists
    groups_key = adata.uns['paga']['groups']
    if ((isinstance(colors, Iterable) and len(colors) == len(adata.obs[groups_key].cat.categories))
        or colors is None or isinstance(colors, str)):
        colors = [colors]

    if frameon is None:
        frameon = settings._frameon

    # labels is a list that contains no lists
    if ((isinstance(labels, Iterable) and len(labels) == len(adata.obs[groups_key].cat.categories))
        or labels is None or isinstance(labels, (str, dict))):
        labels = [labels for i in range(len(colors))]

    if title is None and len(colors) > 1:
        title = [c for c in colors]
    elif isinstance(title, str):
        title = [title for c in colors]
    elif title is None:
        title = [None for c in colors]

    if colorbar is None:
        var_names = adata.var_names if adata.raw is None else adata.raw.var_names
        colorbars = [True if c in var_names else False for c in colors]
    else:
        colorbars = [False for c in colors]

    if isinstance(root, str):
        if root in labels:
            root = list(labels).index(root)
        else:
            raise ValueError(
                'If `root` is a string, it needs to be one of {} not \'{}\'.'
                .format(labels, root))
    if isinstance(root, list) and root[0] in labels:
        root = [list(labels).index(r) for r in root]

    # define the adjacency matrices
    adjacency_solid = adata.uns['paga'][solid_edges].copy()
    adjacency_dashed = None
    if threshold is None:
        threshold = 0.01  # default threshold
    if threshold > 0:
        adjacency_solid.data[adjacency_solid.data < threshold] = 0
        adjacency_solid.eliminate_zeros()
    if dashed_edges is not None:
        adjacency_dashed = adata.uns['paga'][dashed_edges].copy()
        if threshold > 0:
            adjacency_dashed.data[adjacency_dashed.data < threshold] = 0
            adjacency_dashed.eliminate_zeros()

    # compute positions
    if pos is None:
        adj_tree = None
        if layout in {'rt', 'rt_circular', 'eq_tree'}:
            adj_tree = adata.uns['paga']['connectivities_tree']
        pos = _compute_pos(
            adjacency_solid, layout=layout, random_state=random_state, init_pos=init_pos, layout_kwds=layout_kwds, adj_tree=adj_tree, root=root)

    if plot:
        if ax is None:
            axs, panel_pos, draw_region_width, figure_width = utils.setup_axes(
                panels=colors, colorbars=colorbars)
        else:
            axs = ax

        if len(colors) == 1 and not isinstance(axs, list):
            axs = [axs]

        for icolor, c in enumerate(colors):
            if title[icolor] is not None:
                axs[icolor].set_title(title[icolor])
            sct = _paga_graph(
                adata,
                axs[icolor],
                colors=c,
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
                text_kwds=text_kwds,
                node_size_scale=node_size_scale,
                node_size_power=node_size_power,
                edge_width_scale=edge_width_scale,
                min_edge_width=min_edge_width,
                max_edge_width=max_edge_width,
                normalize_to_color=normalize_to_color,
                frameon=frameon,
                cmap=cmap,
                cax=cax,
                colorbar=colorbars[icolor],
                cb_kwds=cb_kwds,
                use_raw=use_raw,
                title=title[icolor],
                export_to_gexf=export_to_gexf,
                single_component=single_component,
                arrowsize=arrowsize,
                pos=pos)
            if colorbars[icolor]:
                bottom = panel_pos[0][0]
                height = panel_pos[1][0] - bottom
                width = 0.006 * draw_region_width / len(colors)
                left = panel_pos[2][2*icolor+1] + 0.2 * width
                rectangle = [left, bottom, width, height]
                fig = pl.gcf()
                ax_cb = fig.add_axes(rectangle)
                cb = pl.colorbar(sct, format=ticker.FuncFormatter(utils.ticks_formatter),
                                 cax=ax_cb)
    if add_pos:
        adata.uns['paga']['pos'] = pos
        logg.hint('added \'pos\', the PAGA positions (adata.uns[\'paga\'])')
    if plot:
        utils.savefig_or_show('paga', show=show, save=save)
        if len(colors) == 1 and isinstance(axs, list): axs = axs[0]
        return axs if show == False else None


def _paga_graph(
        adata,
        ax,
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
        text_kwds=None,
        node_size_scale=1,
        node_size_power=0.5,
        edge_width_scale=1,
        normalize_to_color='reference',
        title=None,
        pos=None,
        cmap=None,
        frameon=True,
        min_edge_width=None,
        max_edge_width=None,
        export_to_gexf=False,
        cax=None,
        colorbar=None,
        use_raw=True,
        cb_kwds={},
        single_component=False,
        arrowsize=30):
    node_labels = labels  # rename for clarity
    if (node_labels is not None
        and isinstance(node_labels, str)
        and node_labels != adata.uns['paga']['groups']):
        raise ValueError('Provide a list of group labels for the PAGA groups {}, not {}.'
                         .format(adata.uns['paga']['groups'], node_labels))
    groups_key = adata.uns['paga']['groups']
    if node_labels is None:
        node_labels = adata.obs[groups_key].cat.categories

    if (colors is None or colors == groups_key) and groups_key is not None:
        if (groups_key + '_colors' not in adata.uns
            or len(adata.obs[groups_key].cat.categories)
               != len(adata.uns[groups_key + '_colors'])):
            utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        colors = adata.uns[groups_key + '_colors']
        for iname, name in enumerate(adata.obs[groups_key].cat.categories):
            if name in settings.categories_to_ignore: colors[iname] = 'grey'

    nx_g_solid = nx.Graph(adjacency_solid)
    if dashed_edges is not None:
        nx_g_dashed = nx.Graph(adjacency_dashed)

    # convert pos to dict
    if isinstance(pos, str):
        if not pos.endswith('.gdf'):
            raise ValueError('Currently only supporting reading positions from .gdf files.'
                             'Consider generating them using, for instance, Gephi.')
        s = ''  # read the node definition from the file
        with open(pos) as f:
            f.readline()
            for line in f:
                if line.startswith('edgedef>'):
                    break
                s += line
        from io import StringIO
        df = pd.read_csv(StringIO(s), header=-1)
        pos = df[[4, 5]].values
    pos_array = pos
    # convert to dictionary
    pos = {n: [p[0], p[1]] for n, p in enumerate(pos)}

    # uniform color
    if isinstance(colors, str) and is_color_like(colors):
        colors = [colors for c in range(len(node_labels))]

    # color degree of the graph
    if isinstance(colors, str) and colors.startswith('degree'):
        # see also tools.paga.paga_degrees
        if colors == 'degree_dashed':
            colors = [d for _, d in nx_g_dashed.degree(weight='weight')]
        elif colors == 'degree_solid':
            colors = [d for _, d in nx_g_solid.degree(weight='weight')]
        else:
            raise ValueError('`degree` either "degree_dashed" or "degree_solid".')
        colors = (np.array(colors) - np.min(colors)) / (np.max(colors) - np.min(colors))

    # plot gene expression
    var_names = adata.var_names if adata.raw is None else adata.raw.var_names
    if isinstance(colors, str) and colors in var_names:
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for icat, cat in enumerate(cats):
            subset = (cat == adata.obs[groups_key]).values
            if adata.raw is not None and use_raw:
                adata_gene = adata.raw[:, colors]
            else:
                adata_gene = adata[:, colors]
            x_color.append(np.mean(adata_gene.X[subset]))
        colors = x_color

    # plot continuous annotation
    if (isinstance(colors, str) and colors in adata.obs
        and not is_categorical_dtype(adata.obs[colors])):
        x_color = []
        cats = adata.obs[groups_key].cat.categories
        for icat, cat in enumerate(cats):
            subset = (cat == adata.obs[groups_key]).values
            x_color.append(adata.obs.loc[subset, colors].mean())
        colors = x_color

    # plot categorical annotation
    if (isinstance(colors, str) and colors in adata.obs and
        is_categorical_dtype(adata.obs[colors])):
        from ... import utils as sc_utils
        asso_names, asso_matrix = sc_utils.compute_association_matrix_of_groups(
            adata, prediction=groups_key, reference=colors,
            normalization='reference' if normalize_to_color else 'prediction')
        utils.add_colors_for_categorical_sample_annotation(adata, colors)
        asso_colors = sc_utils.get_associated_colors_of_groups(
            adata.uns[colors + '_colors'], asso_matrix)
        colors = asso_colors

    if len(colors) < len(node_labels):
        print(node_labels, colors)
        raise ValueError(
            '`color` list need to be at least as long as `groups`/`node_labels` list.')

    # count number of connected components
    n_components, labels = scipy.sparse.csgraph.connected_components(adjacency_solid)
    if n_components > 1 and not single_component:
        logg.msg(
            'Graph has more than a single connected component. '
            'To restrict to this component, pass `single_component=True`.')
    if n_components > 1 and single_component:
        component_sizes = np.bincount(labels)
        largest_component = np.where(
            component_sizes == component_sizes.max())[0][0]
        adjacency_solid = adjacency_solid.tocsr()[labels == largest_component, :]
        adjacency_solid = adjacency_solid.tocsc()[:, labels == largest_component]
        colors = np.array(colors)[labels == largest_component]
        node_labels = np.array(node_labels)[labels == largest_component]
        logg.info(
            'Restricting graph to largest connected component by dropping categories\n'
            '{}'.format(
                adata.obs[groups_key].cat.categories[labels != largest_component].tolist()))
        nx_g_solid = nx.Graph(adjacency_solid)
        if dashed_edges is not None:
            raise ValueError('`single_component` only if `dashed_edges` is `None`.')

    # edge widths
    base_edge_width = edge_width_scale * 5 * rcParams['lines.linewidth']

    # draw dashed edges
    if dashed_edges is not None:
        widths = [x[-1]['weight'] for x in nx_g_dashed.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if max_edge_width is not None:
            widths = np.clip(widths, None, max_edge_width)
        nx.draw_networkx_edges(nx_g_dashed, pos, ax=ax, width=widths, edge_color='grey',
                               style='dashed', alpha=0.5)

    # draw solid edges
    if transitions is None:
        widths = [x[-1]['weight'] for x in nx_g_solid.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nx.draw_networkx_edges(nx_g_solid, pos, ax=ax, width=widths, edge_color='black')
    # draw directed edges
    else:
        adjacency_transitions = adata.uns['paga'][transitions].copy()
        if threshold is None: threshold = 0.01
        adjacency_transitions.data[adjacency_transitions.data < threshold] = 0
        adjacency_transitions.eliminate_zeros()
        g_dir = nx.DiGraph(adjacency_transitions.T)
        widths = [x[-1]['weight'] for x in g_dir.edges(data=True)]
        widths = base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
        nx.draw_networkx_edges(g_dir, pos, ax=ax, width=widths, edge_color='black', arrowsize=arrowsize)

    if export_to_gexf:
        if isinstance(colors[0], tuple):
            from matplotlib.colors import rgb2hex
            colors = [rgb2hex(c) for c in colors]
        for count, n in enumerate(nx_g_solid.nodes()):
            nx_g_solid.node[count]['label'] = str(node_labels[count])
            nx_g_solid.node[count]['color'] = str(colors[count])
            nx_g_solid.node[count]['viz'] = {
                'position': {'x': 1000*pos[count][0],
                             'y': 1000*pos[count][1],
                             'z': 0}}
        filename = settings.writedir + 'paga_graph.gexf'
        logg.msg('exporting to {}'.format(filename), v=1)
        if settings.writedir != '' and not os.path.exists(settings.writedir):
            os.makedirs(settings.writedir)
        nx.write_gexf(nx_g_solid, settings.writedir + 'paga_graph.gexf')

    ax.set_frame_on(frameon)
    ax.set_xticks([])
    ax.set_yticks([])

    # groups sizes
    if (groups_key is not None and groups_key + '_sizes' in adata.uns):
        groups_sizes = adata.uns[groups_key + '_sizes']
    else:
        groups_sizes = np.ones(len(node_labels))
    base_scale_scatter = 2000
    base_pie_size = (base_scale_scatter / (np.sqrt(adjacency_solid.shape[0]) + 10)
                     * node_size_scale)
    median_group_size = np.median(groups_sizes)
    groups_sizes = base_pie_size * np.power(
        groups_sizes / median_group_size, node_size_power)

    if fontsize is None:
        fontsize = rcParams['legend.fontsize']

    # usual scatter plot
    if not isinstance(colors[0], dict):
        n_groups = len(pos_array)
        sct = ax.scatter(
            pos_array[:, 0], pos_array[:, 1],
            c=colors[:n_groups], edgecolors='face', s=groups_sizes, cmap=cmap)
        for count, group in enumerate(node_labels):
            ax.text(pos_array[count, 0], pos_array[count, 1], group,
                    verticalalignment='center',
                    horizontalalignment='center',
                    size=fontsize, fontweight=fontweight, **text_kwds)
    # else pie chart plot
    else:
        # start with this dummy plot... otherwise strange behavior
        sct = ax.scatter(
            pos_array[:, 0], pos_array[:, 1],
            c='white', edgecolors='face', s=groups_sizes, cmap=cmap)
        trans = ax.transData.transform
        bbox = ax.get_position().get_points()
        ax_x_min = bbox[0, 0]
        ax_x_max = bbox[1, 0]
        ax_y_min = bbox[0, 1]
        ax_y_max = bbox[1, 1]
        ax_len_x = ax_x_max - ax_x_min
        ax_len_y = ax_y_max - ax_y_min
        trans2 = ax.transAxes.inverted().transform
        pie_axs = []
        for count, n in enumerate(nx_g_solid.nodes()):
            pie_size = groups_sizes[count] / base_scale_scatter
            x1, y1 = trans(pos[n])     # data coordinates
            xa, ya = trans2((x1, y1))  # axis coordinates
            xa = ax_x_min + (xa - pie_size/2) * ax_len_x
            ya = ax_y_min + (ya - pie_size/2) * ax_len_y
            # clip, the fruchterman layout sometimes places below figure
            if ya < 0: ya = 0
            if xa < 0: xa = 0
            pie_axs.append(pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y], frameon=False))
            pie_axs[count].set_xticks([])
            pie_axs[count].set_yticks([])
            if not isinstance(colors[count], dict):
                raise ValueError('{} is neither a dict of valid matplotlib colors '
                                 'nor a valid matplotlib color.'.format(colors[count]))
            color_single = colors[count].keys()
            fracs = [colors[count][c] for c in color_single]
            if sum(fracs) < 1:
                color_single = list(color_single)
                color_single.append('grey')
                fracs.append(1-sum(fracs))
            pie_axs[count].pie(fracs, colors=color_single)
        if node_labels is not None:
            for ia, a in enumerate(pie_axs):
                a.text(0.5, 0.5, node_labels[ia],
                       verticalalignment='center',
                       horizontalalignment='center',
                       transform=a.transAxes,
                       size=fontsize, fontweight=fontweight, **text_kwds)
    return sct


def paga_path(
    adata,
    nodes,
    keys,
    use_raw=True,
    annotations=['dpt_pseudotime'],
    color_map=None,
    color_maps_annotations={'dpt_pseudotime': 'Greys'},
    palette_groups=None,
    n_avg=1,
    groups_key=None,
    xlim=[None, None],
    title=None,
    left_margin=None,
    ytick_fontsize=None,
    title_fontsize=None,
    show_node_names=True,
    show_yticks=True,
    show_colorbar=True,
    legend_fontsize=None,
    legend_fontweight=None,
    normalize_to_zero_one=False,
    as_heatmap=True,
    return_data=False,
    show=None,
    save=None,
    ax=None,
) -> Optional[Axes]:
    """Gene expression and annotation changes along paths in the abstracted graph.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        An annotated data matrix.
    nodes : list of group names or their category indices
        A path through nodes of the abstracted graph, that is, names or indices
        (within `.categories`) of groups that have been used to run PAGA.
    keys : list of str
        Either variables in `adata.var_names` or annotations in
        `adata.obs`. They are plotted using `color_map`.
    use_raw : `bool`, optional (default: `True`)
        Use `adata.raw` for retrieving gene expressions if it has been set.
    annotations : list of annotations, optional (default: ['dpt_pseudotime'])
        Plot these keys with `color_maps_annotations`. Need to be keys for
        `adata.obs`.
    color_map : color map for plotting keys or `None`, optional (default: `None`)
        Matplotlib colormap.
    color_maps_annotations : dict storing color maps or `None`, optional (default: {'dpt_pseudotime': 'Greys'})
        Color maps for plotting the annotations. Keys of the dictionary must
        appear in `annotations`.
    palette_groups : list of colors or `None`, optional (default: `None`)
        Ususally, use the same `sc.pl.palettes...` as used for coloring the
        abstracted graph.
    n_avg : `int`, optional (default: 1)
        Number of data points to include in computation of running average.
    groups_key : `str`, optional (default: `None`)
        Key of the grouping used to run PAGA. If `None`, defaults to
        `adata.uns['paga']['groups']`.
    as_heatmap : `bool`, optional (default: `True`)
        Plot the timeseries as heatmap. If not plotting as heatmap,
        `annotations` have no effect.
    show_node_names : `bool`, optional (default: `True`)
        Plot the node names on the nodes bar.
    show_colorbar : `bool`, optional (default: `True`)
        Show the colorbar.
    show_yticks : `bool`, optional (default: `True`)
        Show the y ticks.
    normalize_to_zero_one : `bool`, optional (default: `True`)
        Shift and scale the running average to [0, 1] per gene.
    return_data : `bool`, optional (default: `False`)
        Return the timeseries data in addition to the axes if `True`.
    show : `bool`, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \\{'.pdf', '.png', '.svg'\\}.
    ax : :class:`~matplotlib.axes.Axes`
         A matplotlib axes object.

    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object, if `ax` is `None`, else `None`.
    If `return_data`, return the timeseries data in addition to an axes.
    """
    ax_was_none = ax is None

    if groups_key is None:
        if 'groups' not in adata.uns['paga']:
            raise KeyError(
                'Pass the key of the grouping with which you ran PAGA, '
                'using the parameter `groups_key`.')
        groups_key = adata.uns['paga']['groups']
    groups_names = adata.obs[groups_key].cat.categories

    if 'dpt_pseudotime' not in adata.obs.keys():
        raise ValueError(
            '`pl.paga_path` requires computation of a pseudotime `tl.dpt` '
            'for ordering at single-cell resolution')

    if palette_groups is None:
        utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        palette_groups = adata.uns[groups_key + '_colors']

    def moving_average(a):
        return sc_utils.moving_average(a, n_avg)

    ax = pl.gca() if ax is None else ax
    from matplotlib import transforms
    trans = transforms.blended_transform_factory(
        ax.transData, ax.transAxes)
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
                raise ValueError(
                    'Each node/group needs to be one of {} (`groups_key`=\'{}\') not \'{}\'.'
                    .format(groups_names.tolist(), groups_key, node))
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
                adata.obs[groups_key].values == nodes_strs[igroup]]
            if len(idcs) == 0:
                raise ValueError(
                    'Did not find data points that match '
                    '`adata.obs[{}].values == str({})`.'
                    'Check whether adata.obs[{}] actually contains what you expect.'
                    .format(groups_key, group, groups_key))
            idcs_group = np.argsort(adata.obs['dpt_pseudotime'].values[
                adata.obs[groups_key].values == nodes_strs[igroup]])
            idcs = idcs[idcs_group]
            if key in adata.obs_keys(): x += list(adata.obs[key].values[idcs])
            else: x += list(adata_X[:, key].X[idcs])
            if ikey == 0:
                groups += [group for i in range(len(idcs))]
                x_tick_locs.append(len(x))
                for anno in annotations:
                    series = adata.obs[anno]
                    if is_categorical_dtype(series): series = series.cat.codes
                    anno_dict[anno] += list(series.values[idcs])
        if n_avg > 1:
            old_len_x = len(x)
            x = moving_average(x)
            if ikey == 0:
                for key in annotations:
                    if not isinstance(anno_dict[key][0], str):
                        anno_dict[key] = moving_average(anno_dict[key])
        if normalize_to_zero_one:
            x -= np.min(x)
            x /= np.max(x)
        X.append(x)
        if not as_heatmap:
            ax.plot(x[xlim[0]:xlim[1]], label=key)
        if ikey == 0:
            for igroup, group in enumerate(nodes):
                if len(groups_names) > 0 and group not in groups_names:
                    label = groups_names[group]
                else:
                    label = group
                x_tick_labels.append(label)
    X = np.array(X)
    if as_heatmap:
        img = ax.imshow(X, aspect='auto', interpolation='nearest',
                        cmap=color_map)
        if show_yticks:
            ax.set_yticks(range(len(X)))
            ax.set_yticklabels(keys, fontsize=ytick_fontsize)
        else:
            ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.tick_params(axis='both', which='both', length=0)
        ax.grid(False)
        if show_colorbar:
            pl.colorbar(img, ax=ax)
        left_margin = 0.2 if left_margin is None else left_margin
        pl.subplots_adjust(left=left_margin)
    else:
        left_margin = 0.4 if left_margin is None else left_margin
        if len(keys) > 1:
            pl.legend(frameon=False, loc='center left',
                      bbox_to_anchor=(-left_margin, 0.5),
                      fontsize=legend_fontsize)
    xlabel = groups_key
    if not as_heatmap:
        ax.set_xlabel(xlabel)
        pl.yticks([])
        if len(keys) == 1: pl.ylabel(keys[0] + ' (a.u.)')
    else:
        import matplotlib.colors
        # groups bar
        ax_bounds = ax.get_position().bounds
        groups_axis = pl.axes([ax_bounds[0],
                               ax_bounds[1] - ax_bounds[3] / len(keys),
                               ax_bounds[2],
                               ax_bounds[3] / len(keys)])
        groups = np.array(groups)[None, :]
        groups_axis.imshow(groups, aspect='auto',
                           interpolation="nearest",
                           cmap=matplotlib.colors.ListedColormap(
                               # the following line doesn't work because of normalization
                               # adata.uns['paga_groups_colors'])
                               palette_groups[np.min(groups).astype(int):],
                               N=int(np.max(groups)+1-np.min(groups))))
        if show_yticks:
            groups_axis.set_yticklabels(['', xlabel, ''], fontsize=ytick_fontsize)
        else:
            groups_axis.set_yticks([])
        groups_axis.set_frame_on(False)
        if show_node_names:
            ypos = (groups_axis.get_ylim()[1] + groups_axis.get_ylim()[0])/2
            x_tick_locs = sc_utils.moving_average(x_tick_locs, n=2)
            for ilabel, label in enumerate(x_tick_labels):
                groups_axis.text(x_tick_locs[ilabel], ypos, x_tick_labels[ilabel],
                                 fontdict={'horizontalalignment': 'center',
                                           'verticalalignment': 'center'})
        groups_axis.set_xticks([])
        groups_axis.grid(False)
        groups_axis.tick_params(axis='both', which='both', length=0)
        # further annotations
        y_shift = ax_bounds[3] / len(keys)
        for ianno, anno in enumerate(annotations):
            if ianno > 0: y_shift = ax_bounds[3] / len(keys) / 2
            anno_axis = pl.axes([ax_bounds[0],
                                 ax_bounds[1] - (ianno+2) * y_shift,
                                 ax_bounds[2],
                                 y_shift])
            arr = np.array(anno_dict[anno])[None, :]
            if anno not in color_maps_annotations:
                color_map_anno = ('Vega10' if is_categorical_dtype(adata.obs[anno])
                                  else 'Greys')
            else:
                color_map_anno = color_maps_annotations[anno]
            img = anno_axis.imshow(arr, aspect='auto',
                                   interpolation='nearest',
                                   cmap=color_map_anno)
            if show_yticks:
                anno_axis.set_yticklabels(['', anno, ''],
                                          fontsize=ytick_fontsize)
                anno_axis.tick_params(axis='both', which='both', length=0)
            else:
                anno_axis.set_yticks([])
            anno_axis.set_frame_on(False)
            anno_axis.set_xticks([])
            anno_axis.grid(False)
    if title is not None: ax.set_title(title, fontsize=title_fontsize)
    if show is None and not ax_was_none: show = False
    else: show = settings.autoshow if show is None else show
    utils.savefig_or_show('paga_path', show=show, save=save)
    if return_data:
        df = pd.DataFrame(data=X.T, columns=keys)
        df['groups'] = moving_average(groups)  # groups is without moving average, yet
        if 'dpt_pseudotime' in anno_dict:
            df['distance'] = anno_dict['dpt_pseudotime'].T
        return ax, df if ax_was_none and show == False else df
    else:
        return ax if ax_was_none and show == False else None


def paga_adjacency(
        adata,
        adjacency='connectivities',
        adjacency_tree='connectivities_tree',
        as_heatmap=True,
        color_map=None,
        show=None,
        save=None):
    """Connectivity of paga groups.
    """
    connectivity = adata.uns[adjacency].toarray()
    connectivity_select = adata.uns[adjacency_tree]
    if as_heatmap:
        matrix(connectivity, color_map=color_map, show=False)
        for i in range(connectivity_select.shape[0]):
            neighbors = connectivity_select[i].nonzero()[1]
            pl.scatter([i for j in neighbors], neighbors, color='black', s=1)
    # as a stripplot
    else:
        pl.figure()
        for i, cs in enumerate(connectivity):
            x = [i for j, d in enumerate(cs) if i != j]
            y = [c for j, c in enumerate(cs) if i != j]
            pl.scatter(x, y, color='gray', s=1)
            neighbors = connectivity_select[i].nonzero()[1]
            pl.scatter([i for j in neighbors],
                       cs[neighbors], color='black', s=1)
    utils.savefig_or_show('paga_connectivity', show=show, save=save)
