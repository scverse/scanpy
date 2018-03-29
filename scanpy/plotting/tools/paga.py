import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
import networkx as nx
from matplotlib import pyplot as pl
from matplotlib.colors import is_color_like
from matplotlib import rcParams

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
    paga_scatter(adata,
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
         groups=groups_graph, color=color_graph, **paga_graph_params)
    if suptitle is not None: pl.suptitle(suptitle)
    utils.savefig_or_show('paga_compare', show=show, save=save)
    if show == False: return axs


def paga_scatter(
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
        solid_edges='confidence',
        dashed_edges=None,
        layout=None,
        root=0,
        groups=None,
        color=None,
        threshold_solid=None,
        threshold_dashed=1e-6,
        fontsize=None,
        node_size_scale=1,
        node_size_power=0.5,
        edge_width_scale=1,
        min_edge_width=None,
        max_edge_width=None,
        title='abstracted graph',
        left_margin=0.01,
        random_state=0,
        pos=None,
        cmap=None,
        frameon=True,
        rootlevel=None,
        return_pos=False,
        export_to_gexf=False,
        show=None,
        save=None,
        ax=None):
    """Plot the abstracted graph.

    This uses igraph's layout algorithms for most layouts [Csardi06]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    solid_edges : `str`, optional (default: 'paga_confidence')
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges : `str` or `None`, optional (default: `None`)
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    layout : {'fr', 'rt', 'rt_circular', 'eq_tree', ...}, optional (default: 'fr')
        Plotting layout. 'fr' stands for Fruchterman-Reingold, 'rt' stands for
        Reingold Tilford. 'eq_tree' stands for 'eqally spaced tree'. All but
        'eq_tree' use the igraph layout function. All other igraph layouts are
        also permitted. See also parameter `pos`.
    random_state : `int` or `None`, optional (default: 0)
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, the initial
        state is not reproducible.
    root : int, str or list of int, optional (default: 0)
        If choosing a tree layout, this is the index of the root node or root
        nodes. If this is a non-empty vector then the supplied node IDs are used
        as the roots of the trees (or a single tree if the graph is
        connected. If this is `None` or an empty list, the root vertices are
        automatically calculated based on topological sorting.
    groups : `str`, `list`, `dict`
        The node (groups) labels.
    color : color string or iterable, {'degree_dashed', 'degree_solid'}, optional (default: None)
        The node colors.  Besides cluster colors, lists and uniform colors this
        also acceppts {'degree_dashed', 'degree_solid'} which are plotted using
        continuous color map.
    threshold_solid : `float` or `None`, optional (default: `None`)
        Do not draw edges for weights below this threshold. Set to `None` if you
        want all edges.
    threshold_dashed : `float` or `None`, optional (default: 1e-6)
        Do not draw edges for weights below this threshold. Set to `None` if you
        want all edges.
    fontsize : int (default: None)
        Font size for node labels.
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
    pos : array-like, filename of `.gdf` file,  optional (default: `None`)
        Two-column array/list storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    export_to_gexf : `bool`, optional (default: `None`)
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    return_pos : `bool`, optional (default: `False`)
        Return the positions.
    title : `str`, optional (default: `None`)
         Provide title for panels either as `['title1', 'title2', ...]` or
         `'title1,title2,...'`.
    frameon : `bool`, optional (default: `True`)
         Draw a frame around the abstracted graph.
    show : `bool`, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : `matplotlib.Axes`
         A matplotlib axes object.

    Returns
    -------
    Adds `'paga_pos'` to `adata.uns`.

    If `show==False`, a list of `matplotlib.Axis` objects. Every second element
    corresponds to the 'right margin' drawing area for color bars and legends.

    If `return_pos` is `True`, in addition, the positions of the nodes.
    """

    # colors is a list that contains no lists
    if isinstance(color, list) and True not in [isinstance(c, list) for c in color]: color = [color]
    if color is None or isinstance(color, str): color = [color]

    # groups is a list that contains no lists
    if isinstance(groups, list) and True not in [isinstance(g, list) for g in groups]: groups = [groups]
    if groups is None or isinstance(groups, dict) or isinstance(groups, str): groups = [groups]

    if title is None or isinstance(title, str): title = [title for name in groups]

    if ax is None:
        axs, _, _, _ = utils.setup_axes(panels=color)
    else:
        axs = ax
    if len(color) == 1 and not isinstance(axs, list): axs = [axs]

    for icolor, c in enumerate(color):
        pos = _paga_graph(
            adata,
            axs[icolor],
            solid_edges=solid_edges,
            dashed_edges=dashed_edges,
            threshold_solid=threshold_solid,
            threshold_dashed=threshold_dashed,
            layout=layout,
            root=root,
            rootlevel=rootlevel,
            color=c,
            groups=groups[icolor],
            fontsize=fontsize,
            node_size_scale=node_size_scale,
            node_size_power=node_size_power,
            edge_width_scale=edge_width_scale,
            min_edge_width=min_edge_width,
            max_edge_width=max_edge_width,
            frameon=frameon,
            cmap=cmap,
            title=title[icolor],
            random_state=random_state,
            export_to_gexf=export_to_gexf,
            pos=pos)
    adata.uns['paga']['pos'] = pos
    utils.savefig_or_show('paga_graph', show=show, save=save)
    if len(color) == 1 and isinstance(axs, list): axs = axs[0]
    if return_pos:
        return (axs, pos) if show == False else pos
    else:
        return axs if show == False else None


def _paga_graph(
        adata,
        ax,
        solid_edges=None,
        dashed_edges=None,
        threshold_solid=None,
        threshold_dashed=1e-6,
        root=0,
        rootlevel=None,
        color=None,
        groups=None,
        fontsize=None,
        node_size_scale=1,
        node_size_power=0.5,
        edge_width_scale=1,
        title=None,
        layout=None,
        pos=None,
        cmap=None,
        frameon=True,
        min_edge_width=None,
        max_edge_width=None,
        export_to_gexf=False,
        random_state=0):
    node_labels = groups
    if (node_labels is not None
        and isinstance(node_labels, str)
        and node_labels != adata.uns['paga']['groups']):
        raise ValueError('Provide a list of group labels for the PAGA groups {}, not {}.'
                         .format(adata.uns['paga']['groups'], node_labels))
    groups_key = adata.uns['paga']['groups']
    if node_labels is None:
        node_labels = adata.obs[groups_key].cat.categories

    if color is None and groups_key is not None:
        if (groups_key + '_colors' not in adata.uns
            or len(adata.obs[groups_key].cat.categories)
               != len(adata.uns[groups_key + '_colors'])):
            utils.add_colors_for_categorical_sample_annotation(adata, groups_key)
        color = adata.uns[groups_key + '_colors']
        for iname, name in enumerate(adata.obs[groups_key].cat.categories):
            if name in settings.categories_to_ignore: color[iname] = 'grey'

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
    if threshold_solid is not None:
        adjacency_solid[adjacency_solid < threshold_solid] = 0
    nx_g_solid = nx.Graph(adjacency_solid)
    if dashed_edges is not None:
        adjacency_dashed = adata.uns['paga'][dashed_edges].copy()
        if threshold_dashed is not None:
            adjacency_dashed[adjacency_dashed < threshold_dashed] = 0
        nx_g_dashed = nx.Graph(adjacency_dashed)

    # degree of the graph for coloring
    if isinstance(color, str) and color.startswith('degree'):
        # see also tools.paga.paga_degrees
        if color == 'degree_dashed':
            color = [d for _, d in nx_g_dashed.degree(weight='weight')]
        elif color == 'degree_solid':
            color = [d for _, d in nx_g_solid.degree(weight='weight')]
        else:
            raise ValueError('`degree` either "degree_dashed" or "degree_solid".')
        color = (np.array(color) - np.min(color)) / (np.max(color) - np.min(color))

    # plot numeric colors
    colorbar = False
    if isinstance(color, (list, np.ndarray)) and not isinstance(color[0], (str, dict)):
        import matplotlib
        norm = matplotlib.colors.Normalize()
        color = norm(color)
        if cmap is None: cmap = rcParams['image.cmap']
        cmap = matplotlib.cm.get_cmap(cmap)
        color = [cmap(c) for c in color]
        colorbar = True

    if len(color) < len(node_labels):
        print(node_labels, color)
        raise ValueError('`color` list need to be at least as long as `node_labels` list.')

    # node positions from adjacency_solid
    if pos is None:
        if layout is None:
            layout = 'fr'
        # igraph layouts
        if layout != 'eq_tree':
            from ... import utils as sc_utils
            adj_solid_weights = adjacency_solid
            g = sc_utils.get_igraph_from_adjacency(adj_solid_weights)
            if 'rt' in layout:
                g_tree = g
                if solid_edges != 'confidence_tree':
                    adj_tree = adata.uns['paga']['confidence_tree']
                    g_tree = sc_utils.get_igraph_from_adjacency(adj_tree)
                pos_list = g_tree.layout(
                    layout, root=root if isinstance(root, list) else [root],
                    rootlevel=rootlevel).coords
            elif layout == 'circle':
                pos_list = g.layout(layout).coords
            else:
                np.random.seed(random_state)
                init_coords = np.random.random((adjacency_solid.shape[0], 2)).tolist()
                pos_list = g.layout(layout, seed=init_coords, weights='weight').coords
            pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
        # equally-spaced tree
        else:
            nx_g_tree = nx_g_solid
            if solid_edges != 'confidence_tree':
                adj_tree = adata.uns['paga']['confidence_tree']
                nx_g_tree = nx.Graph(adj_tree)
            pos = utils.hierarchy_pos(nx_g_tree, root)
            if len(pos) < adjacency_solid.shape[0]:
                raise ValueError('This is a forest and not a single tree. '
                                 'Try another `layout`, e.g., {\'fr\'}.')
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
    widths = [x[-1]['weight'] for x in nx_g_solid.edges(data=True)]
    widths = base_edge_width * np.array(widths)
    if min_edge_width is not None or max_edge_width is not None:
        widths = np.clip(widths, min_edge_width, max_edge_width)
    nx.draw_networkx_edges(nx_g_solid, pos, ax=ax, width=widths, edge_color='black')

    if export_to_gexf:
        for count, n in enumerate(nx_g_solid.nodes()):
            nx_g_solid.node[count]['label'] = node_labels[count]
            nx_g_solid.node[count]['color'] = color[count]
            nx_g_solid.node[count]['viz'] = {
                'position': {'x': 1000*pos[count][0],
                             'y': 1000*pos[count][1],
                             'z': 0}}
        logg.msg('exporting to {}'.format(settings.writedir + 'paga_graph.gexf'), v=1)
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
    if is_color_like(color[0]):
        ax.scatter(pos_array[:, 0], pos_array[:, 1],
                   c=color, edgecolors='face', s=groups_sizes)
        for count, group in enumerate(node_labels):
            ax.text(pos_array[count, 0], pos_array[count, 1], group,
                verticalalignment='center',
                horizontalalignment='center', size=fontsize)
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
            a = pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
            if not isinstance(color[count], dict):
                raise ValueError('{} is neither a dict of valid matplotlib colors '
                                 'nor a valid matplotlib color.'.format(color[count]))
            color_single = color[count].keys()
            fracs = [color[count][c] for c in color_single]
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
    if colorbar:
        ax1 = pl.axes([0.95, 0.1, 0.03, 0.7])
        cb = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap,
                                              norm=norm)
    return pos_array


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
