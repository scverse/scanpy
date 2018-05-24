import os
import numpy as np
import pandas as pd
import scipy
from pandas.api.types import is_categorical_dtype
import networkx as nx
from matplotlib import pyplot as pl
from matplotlib.colors import is_color_like
from matplotlib import rcParams
from collections import Iterable

from .. import utils
from ... import utils as sc_utils
from ... import settings
from ... import logging as logg
from ..utils import matrix


def paga_compare(
        adata,
        basis='tsne',
        color=None,
        alpha=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='on data',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        size=None,
        title=None,
        right_margin=None,
        left_margin=0.05,
        show=None,
        save=None,
        title_graph=None,
        groups_graph=None,
        color_graph=None,
        **paga_graph_params):
    """Statisical graph abstraction.

    Consists in a scatter plot and the abstracted graph. See
    :func:`~sanpy.api.pl.paga_scatter` and :func:`~scanpy.api.pl.paga_graph` for
    most of the parameters.

    See :func:`~scanpy.api.pl.paga_path` for visualizing gene changes along paths
    through the abstracted graph.

    Additional parameters are as follows.

    Parameters
    ----------
    title : `str` or `None`, optional (default: `None`)
        Title for the scatter panel, or, if `title_graph is None`, title for the
        whole figure.
    title_graph : `str` or `None`, optional (default: `None`)
        Separate title for the abstracted graph.
    """
    axs, _, _, _ = utils.setup_axes(panels=[0, 1],
                                    right_margin=right_margin)  # dummy colors
    # set a common title for the figure
    suptitle = None
    if title is not None and title_graph is None:
        suptitle = title
        title = ''
        title_graph = ''
    elif title_graph is None:
        title_graph = 'abstracted graph'
    _paga_scatter(adata,
                basis=basis,
                color=color,
                alpha=alpha,
                groups=groups,
                components=components,
                projection=projection,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                legend_fontweight=legend_fontweight,
                color_map=color_map,
                palette=palette,
                right_margin=None,
                size=size,
                title=title,
                ax=axs[0],
                show=False,
                save=False)
    paga(adata, ax=axs[1], show=False, save=False, title=title_graph,
         labels=groups_graph, colors=color_graph, **paga_graph_params)
    if suptitle is not None: pl.suptitle(suptitle)
    utils.savefig_or_show('paga_compare', show=show, save=save)
    if show == False: return axs


def _paga_scatter(
        adata,
        basis='tsne',
        color=None,
        alpha=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        size=None,
        title=None,
        right_margin=None,
        show=None,
        save=None,
        ax=None):
    """Scatter plot of paga groups.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    show : bool, optional (default: None)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    If `show==False`, a list of `matplotlib.Axis` objects. Every second element
    corresponds to the 'right margin' drawing area for color bars and legends.
    """
    if color is None:
        color = [adata.uns['paga']['groups']]
    if not isinstance(color, list): color = [color]
    kwds = {}
    if 'draw_graph' in basis:
        from . import draw_graph
        scatter_func = draw_graph
        kwds['edges'] = True
    elif basis == 'umap':
        from . import umap
        scatter_func = umap
        kwds['edges'] = True
    else:
        from . import scatter
        scatter_func = scatter
        kwds['basis'] = basis
    axs = scatter_func(
        adata,
        color=color,
        alpha=alpha,
        groups=groups,
        components=components,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        ax=ax,
        show=False,
        **kwds)
    utils.savefig_or_show('paga_' + basis, show=show, save=save)
    if show == False: return axs


def paga(
        adata,
        threshold=None,
        layout=None,
        layout_kwds={},
        init_pos=None,
        root=0,
        labels=None,
        colors=None,
        single_component=False,
        solid_edges='confidence',
        dashed_edges=None,
        transitions=None,
        threshold_arrows=None,
        threshold_solid=None,
        threshold_dashed=None,
        fontsize=None,
        text_kwds={},
        node_size_scale=1,
        node_size_power=0.5,
        edge_width_scale=1,
        min_edge_width=None,
        max_edge_width=None,
        title=None,
        left_margin=0.01,
        random_state=0,
        pos=None,
        cmap=None,
        cax=None,
        colorbar=None,
        cb_kwds={},
        frameon=True,
        add_pos=True,
        export_to_gexf=False,
        color=None,   # backwards compat
        groups=None,  # backwards compat
        show=None,
        save=None,
        ax=None):
    """Plot the abstracted graph through thresholding low-confidence edges.

    This uses ForceAtlas2 or igraph's layout algorithms for most layouts [Csardi06]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    threshold : `float` or `None`, optional (default: 0.01)
        Do not draw edges for weights below this threshold. Set to 0 if you
        want all edges. Discarding low-confidence edges helps in getting a much
        clearer picture of the graph.
    labels : `None`, `str`, `list`, `dict`, optional (default: `None`)
        The node labels. If `None`, this defaults to the group labels stored in
        the categorical for which :func:`~scanpy.api.tl.paga` has been computed.
    colors : color string or iterable of `int` or `float` or color strings, {'degree_dashed', 'degree_solid'}, optional (default: `None`)
        The node colors. Besides lists, uniform colors this also automatically
        plots the degree of the abstracted graph when passing {'degree_dashed',
        'degree_solid'}.
    layout : {'fa', 'fr', 'rt', 'rt_circular', 'eq_tree', ...}, optional (default: 'fr')
        Plotting layout. 'fa' stands for ForceAtlas2, 'fr' stands for
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
    root : int, str or list of int, optional (default: 0)
        If choosing a tree layout, this is the index of the root node or a list
        of root node indices. If this is a non-empty vector then the supplied
        node IDs are used as the roots of the trees (or a single tree if the
        graph is connected). If this is `None` or an empty list, the root
        vertices are automatically calculated based on topological sorting.
    single_component : `bool`, optional (default: `False`)
        Restrict to largest connected component.
    solid_edges : `str`, optional (default: 'paga_confidence')
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges : `str` or `None`, optional (default: `None`)
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    threshold_solid : `float` or `None`, optional (default: `threshold`)
        Do not draw edges for weights below this threshold. Set to `None` if you
        want all edges.
    threshold_dashed : `float` or `None`, optional (default: `threshold`)
        Do not draw edges for weights below this threshold. Set to `None` if you
        want all edges.
    fontsize : int (default: None)
        Font size for node labels.
    text_kwds : keywords for text
        See `here
        <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.text.html#matplotlib.axes.Axes.text>`_.
    node_size_scale : float (default: 1.0)
        Increase or decrease the size of the nodes.
    node_size_power : float (default: 0.5)
        The power with which groups sizes influence the radius of the nodes.
    edge_width_scale : `float`, optional (default: 5)
        Edge with scale in units of `rcParams['lines.linewidth']`.
    min_edge_width : `float`, optional (default: `None`)
        Min width of solid edges.
    max_edge_width : `float`, optional (default: `None`)
        Max width of solid and dashed edges.
    pos : `np.ndarray`, filename of `.gdf` file,  optional (default: `None`)
        Two-column array/list storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    export_to_gexf : `bool`, optional (default: `None`)
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    cmap : color map
        The color map.
    cax : `matplotlib.Axes`
        A matplotlib axes object for a potential colorbar.
    cb_kwds : colorbar keywords
        See `here
        <https://matplotlib.org/api/colorbar_api.html#matplotlib.colorbar.ColorbarBase>`_,
        for instance, `ticks`.
    add_pos : `bool`, optional (default: `True`)
        Add the positions to `adata.uns['paga']`.
    title : `str`, optional (default: `None`)
         Provide a title.
    frameon : `bool`, optional (default: `True`)
         Draw a frame around the abstracted graph.
    show : `bool`, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : `matplotlib.Axes`
         A matplotlib axes object.

    Note
    ----

    When initializing the positions, note that - for some reason - igraph
    mirrors coordinates along the x axis... that is, you should increase the
    `maxiter` parameter by 1 if the layout is flipped.

    Returns
    -------
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.

    If `show==False`, one or more `matplotlib.Axis` objects. If at least one colorbar is
    drawn, a list with colorbar instances will be returned, too.
    """
    if groups is not None:  # backwards compat
        labels = groups
        logg.warn('`groups` is deprecated in `pl.paga`: use `labels` instead')
    if color is not None:
        colors = color
        logg.warn('`color` is deprecated in `pl.paga`: use `colors` instead')
    # colors is a list that contains no lists
    groups_key = adata.uns['paga']['groups']
    if ((isinstance(colors, Iterable) and len(colors) == len(adata.obs[groups_key].cat.categories))
        or colors is None or isinstance(colors, str)):
        colors = [colors]

    # labels is a list that contains no lists
    if ((isinstance(labels, Iterable) and len(labels) == len(adata.obs[groups_key].cat.categories))
        or labels is None or isinstance(labels, (str, dict))):
        labels = [labels]

    if title is None or isinstance(title, str):
        title = [title for name in labels]

    if ax is None:
        axs, _, _, _ = utils.setup_axes(panels=colors)
    else:
        axs = ax

    if len(colors) == 1 and not isinstance(axs, list):
        axs = [axs]

    for icolor, c in enumerate(colors):
        pos, _ = _paga_graph(
            adata,
            axs[icolor],
            layout=layout,
            layout_kwds=layout_kwds,
            init_pos=init_pos,
            solid_edges=solid_edges,
            dashed_edges=dashed_edges,
            transitions=transitions,
            threshold=threshold,
            threshold_arrows=threshold_arrows,
            threshold_solid=threshold_solid,
            threshold_dashed=threshold_dashed,
            root=root,
            colors=c,
            labels=labels[icolor],
            fontsize=fontsize,
            text_kwds=text_kwds,
            node_size_scale=node_size_scale,
            node_size_power=node_size_power,
            edge_width_scale=edge_width_scale,
            min_edge_width=min_edge_width,
            max_edge_width=max_edge_width,
            frameon=frameon,
            cmap=cmap,
            cax=cax,
            colorbar=colorbar,
            cb_kwds=cb_kwds,
            title=title[icolor],
            random_state=random_state,
            export_to_gexf=export_to_gexf,
            single_component=single_component,
            pos=pos)
    if add_pos:
        adata.uns['paga']['pos'] = pos
        logg.hint('added \'pos\', the PAGA positions (adata.uns[\'paga\'])')
    utils.savefig_or_show('paga', show=show, save=save)
    if len(colors) == 1 and isinstance(axs, list): axs = axs[0]
    return axs if show == False else None


def _paga_graph(
        adata,
        ax,
        layout=None,
        layout_kwds={},
        init_pos=None,
        solid_edges=None,
        dashed_edges=None,
        transitions=None,
        threshold=None,
        threshold_arrows=None,
        threshold_solid=None,
        threshold_dashed=None,
        root=0,
        colors=None,
        labels=None,
        fontsize=None,
        text_kwds=None,
        node_size_scale=1,
        node_size_power=0.5,
        edge_width_scale=1,
        title=None,
        pos=None,
        cmap=None,
        frameon=True,
        min_edge_width=None,
        max_edge_width=None,
        export_to_gexf=False,
        cax=None,
        colorbar=None,
        cb_kwds={},
        single_component=False,
        random_state=0):
    node_labels = labels  # rename for clarity
    if (node_labels is not None
        and isinstance(node_labels, str)
        and node_labels != adata.uns['paga']['groups']):
        raise ValueError('Provide a list of group labels for the PAGA groups {}, not {}.'
                         .format(adata.uns['paga']['groups'], node_labels))
    groups_key = adata.uns['paga']['groups']
    if node_labels is None:
        node_labels = adata.obs[groups_key].cat.categories

    if colors is None and groups_key is not None:
        if (groups_key + '_colors' not in adata.uns
            or len(adata.obs[groups_key].cat.categories)
               != len(adata.uns[groups_key + '_colors'])):
            utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        colors = adata.uns[groups_key + '_colors']
        for iname, name in enumerate(adata.obs[groups_key].cat.categories):
            if name in settings.categories_to_ignore: colors[iname] = 'grey'

    if isinstance(root, str):
        if root in node_labels:
            root = list(node_labels).index(root)
        else:
            raise ValueError(
                'If `root` is a string, it needs to be one of {} not \'{}\'.'
                .format(node_labels.tolist(), root))
    if isinstance(root, list) and root[0] in node_labels:
        root = [list(node_labels).index(r) for r in root]

    # define the objects
    adjacency_solid = adata.uns['paga'][solid_edges].copy()
    # set the the thresholds, either explicitly
    if threshold is not None:
        threshold_solid = threshold
        threshold_dashed = threshold
    # or to a default value
    else:
        if threshold_solid is None:
            threshold_solid = 0.01  # default threshold
        if threshold_dashed is None:
            threshold_dashed = 0.01  # default treshold
    if threshold_solid > 0:
        adjacency_solid.data[adjacency_solid.data < threshold_solid] = 0
        adjacency_solid.eliminate_zeros()
    nx_g_solid = nx.Graph(adjacency_solid)
    if dashed_edges is not None:
        adjacency_dashed = adata.uns['paga'][dashed_edges].copy()
        if threshold_dashed > 0:
            adjacency_dashed.data[adjacency_dashed.data < threshold_dashed] = 0
            adjacency_dashed.eliminate_zeros()
        nx_g_dashed = nx.Graph(adjacency_dashed)

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

    # plot numeric colors
    if isinstance(colors, Iterable) and not isinstance(colors[0], (str, dict)):
        import matplotlib
        norm = matplotlib.colors.Normalize()
        colors = norm(colors)
        if cmap is None: cmap = rcParams['image.cmap']
        cmap = matplotlib.cm.get_cmap(cmap)
        colors = [cmap(c) for c in colors]
        colorbar = True if colorbar is None else colorbar
    else:
        colorbar = False

    if len(colors) < len(node_labels):
        print(node_labels, colors)
        raise ValueError(
            '`color` list need to be at least as long as `goups`/`node_labels` list.')

    # count number of connected components
    n_components, labels = scipy.sparse.csgraph.connected_components(adjacency_solid)
    if n_components > 1 and not single_component:
        logg.info(
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

    # node positions from adjacency_solid
    if pos is None:
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
            nx_g_tree = nx_g_solid
            if solid_edges == 'connectivities':
                adj_tree = adata.uns['paga']['confidence_tree']
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
                g_tree = g
                if solid_edges == 'connectivities':
                    adj_tree = adata.uns['paga']['confidence_tree']
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
                pos_list = g.layout(
                    layout, seed=init_coords,
                    weights='weight', **layout_kwds).coords
            pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
        pos_array = np.array([pos[n] for count, n in enumerate(nx_g_solid)])
    else:
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

    if len(pos) == 1: pos[0] = (0.5, 0.5)

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
        nx.draw_networkx_edges(nx_g_solid, pos, ax=ax, width=widths, edge_color='black')
    # draw directed edges
    else:
        adjacency_transitions = adata.uns['paga'][transitions].copy()
        if threshold_arrows is None:
            threshold_arrows = 0.005
        adjacency_transitions.data[adjacency_transitions.data < threshold_arrows] = 0
        adjacency_transitions.eliminate_zeros()
        g_dir = nx.DiGraph(adjacency_transitions.T)
        widths = [x[-1]['weight'] for x in g_dir.edges(data=True)]
        widths = 100 * base_edge_width * np.array(widths)
        if min_edge_width is not None or max_edge_width is not None:
            widths = np.clip(widths, min_edge_width, max_edge_width)
        nx.draw_networkx_edges(g_dir, pos, ax=ax, width=widths, edge_color='black')

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

    # deal with empty graph
    # ax.plot(pos_array[:, 0], pos_array[:, 1], '.', c='white')

    # draw the nodes (pie charts)
    trans = ax.transData.transform
    bbox = ax.get_position().get_points()
    ax_x_min = bbox[0, 0]
    ax_x_max = bbox[1, 0]
    ax_y_min = bbox[0, 1]
    ax_y_max = bbox[1, 1]
    ax_len_x = ax_x_max - ax_x_min
    ax_len_y = ax_y_max - ax_y_min
    # print([ax_x_min, ax_x_max, ax_y_min, ax_y_max])
    # print([ax_len_x, ax_len_y])
    trans2 = ax.transAxes.inverted().transform
    ax.set_frame_on(frameon)
    ax.set_xticks([])
    ax.set_yticks([])
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
    # usual scatter plot
    if is_color_like(colors[0]):
        scatter = ax.scatter(pos_array[:, 0], pos_array[:, 1],
                             c=colors, edgecolors='face', s=groups_sizes)
        for count, group in enumerate(node_labels):
            ax.text(pos_array[count, 0], pos_array[count, 1], group,
                    verticalalignment='center',
                     horizontalalignment='center', size=fontsize, **text_kwds)
    # else pie chart plot
    else:
        force_labels_to_front = True  # TODO: solve this differently!
        for count, n in enumerate(nx_g_solid.nodes()):
            pie_size = groups_sizes[count] / base_scale_scatter
            xx, yy = trans(pos[n])     # data coordinates
            xa, ya = trans2((xx, yy))  # axis coordinates
            xa = ax_x_min + (xa - pie_size/2) * ax_len_x
            ya = ax_y_min + (ya - pie_size/2) * ax_len_y
            # clip, the fruchterman layout sometimes places below figure
            if ya < 0: ya = 0
            if xa < 0: xa = 0
            a = ax.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
            if not isinstance(colors[count], dict):
                raise ValueError('{} is neither a dict of valid matplotlib colors '
                                 'nor a valid matplotlib color.'.format(colors[count]))
            color_single = colors[count].keys()
            fracs = [colors[count][c] for c in color_single]
            if sum(fracs) < 1:
                color_single = list(color_single)
                color_single.append('grey')
                fracs.append(1-sum(fracs))
            a.pie(fracs, colors=color_single)
            if not force_labels_to_front and node_labels is not None:
                a.text(0.5, 0.5, node_labels[count],
                       verticalalignment='center',
                       horizontalalignment='center',
                       transform=a.transAxes,
                       size=fontsize)
        # TODO: this is a terrible hack, but if we use the solution above (`not
        # force_labels_to_front`), labels get hidden behind pies
        if force_labels_to_front and node_labels is not None:
            for count, n in enumerate(nx_g_solid.nodes()):
                pie_size = groups_sizes[count] / base_scale_scatter
                # all copy and paste from above
                xx, yy = trans(pos[n])     # data coordinates
                xa, ya = trans2((xx, yy))  # axis coordinates
                # make sure a new axis is created
                xa = ax_x_min + (xa - pie_size/2.0000001) * ax_len_x
                ya = ax_y_min + (ya - pie_size/2.0000001) * ax_len_y
                # clip, the fruchterman layout sometimes places below figure
                if ya < 0: ya = 0
                if xa < 0: xa = 0
                a = pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
                a.set_frame_on(False)
                a.set_xticks([])
                a.set_yticks([])
                a.text(0.5, 0.5, node_labels[count],
                       verticalalignment='center',
                       horizontalalignment='center',
                       transform=a.transAxes, size=fontsize)
    if title is not None: ax.set_title(title)
    cb = None
    if colorbar:
        cax = pl.axes([0.95, 0.1, 0.03, 0.8]) if cax is None else cax
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                              norm=norm, **cb_kwds)
    return pos_array, cb


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
        ax=None):
    """Gene expression and annotation changes along paths in the abstracted graph.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
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
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : `matplotlib.Axes`
         A matplotlib axes object.

    Returns
    -------
    A `matplotlib.Axes`, if `ax` is `None`, else `None`. If `return_data`,
    return the timeseries data in addition to an axes.
    """
    ax_was_none = ax is None

    if groups_key is None:
        if 'groups' not in adata.uns['paga']:
            raise KeyError(
                'Pass the key of the grouping with which you ran PAGA, '
                'using the parameter `groups_key`.')
        groups_key = adata.uns['paga']['groups']
    groups_names = adata.obs[groups_key].cat.categories

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
        adjacency='paga_confidence',
        adjacency_tree='paga_confidence_tree',
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
