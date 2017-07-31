# Author: F. Alex Wolf (http://falexwolf.de)
"""Plotting

Plotting functions for each tool and toplevel plotting functions for AnnData.
"""

import numpy as np
import networkx as nx
from matplotlib import pyplot as pl
from matplotlib.colors import is_color_like
from matplotlib.figure import SubplotParams as sppars
from matplotlib import rcParams

from . import utils
from .. import settings as sett
from .. import logging as logg

from .ann_data import scatter, violin
from .ann_data import ranking
from .utils import matrix
from .utils import timeseries, timeseries_subplot, timeseries_as_heatmap


# ------------------------------------------------------------------------------
# Visualization tools
# ------------------------------------------------------------------------------


def pca(adata, **params):
    """Plot PCA results.

    The parameters are the ones of the scatter plot. Call pca_ranking separately
    if you want to change the default settings.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {'2d', '3d'}, optional (default: '2d')
         Projection of plot.
    legend_loc : str, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    """
    show = params['show'] if 'show' in params else None
    if 'show' in params: del params['show']
    pca_scatter(adata, **params, show=False)
    pca_loadings(adata, show=False)
    pca_variance_ratio(adata, show=show)


def pca_scatter(
        adata,
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in PCA coordinates.

    See parameters of sc.pl.pca(). In addition.

    Parameters
    ----------
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    from ..utils import check_adata
    adata = check_adata(adata, verbosity=-1)
    axs = scatter(
        adata,
        basis='pca',
        color=color,
        groups=groups,
        components=components,
        projection=projection,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False, ax=ax)
    utils.savefig_or_show('pca_scatter', show=show, save=save)
    return axs


def pca_loadings(adata, components=None, show=None, save=None):
    """Rank genes according to contributions to PCs.

    Parameters
    ----------
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    """
    if isinstance(components, str): components = components.split(',')
    keys = ['PC1', 'PC2', 'PC3'] if components is None else ['PC{}'.format(c) for c in components]
    ranking(adata, 'var', keys)
    utils.savefig_or_show('pca_loadings', show=show, save=save)


def pca_variance_ratio(adata, log=False, show=None, save=None):
    """Plot the variance ratio.

    Parameters
    ----------
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    """
    ranking(adata, 'add', 'pca_variance_ratio', labels='PC', log=log)
    utils.savefig_or_show('pca_ranking_variance', show=show, save=save)


def diffmap(
        adata,
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {'2d', '3d'}, optional (default: '2d')
         Projection of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object. Only works if plotting a single component.
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    if components == 'all':
        components_list = ['{},{}'.format(*((i, i+1) if i % 2 == 1 else (i+1, i)))
                      for i in range(1, adata.smp['X_diffmap'].shape[1])]
    else:
        if components is None: components = '1,2' if '2d' in projection else '1,2,3'
        if not isinstance(components, list): components_list = [components]
        else: components_list = components
    for components in components_list:
        axs = scatter(
            adata,
            basis='diffmap',
            color=color,
            groups=groups,
            components=components,
            projection=projection,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            color_map=color_map,
            palette=palette,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False,
            save=False,
            ax=ax)
        writekey = 'diffmap'
        if isinstance(components, list): components = ','.join([str(comp) for comp in components])
        writekey += '_components' + components.replace(',', '')
        if sett.savefigs or (save is not None): utils.savefig(writekey)  # TODO: cleaner
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()
    return axs


def draw_graph(
        adata,
        layout=None,
        color=None,
        groups=None,
        components=None,
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """Scatter in tSNE basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layout : {'fr', 'drl', ...}, optional (default: last computed)
        One of the `draw_graph` layouts, see sc.tl.draw_graph. By default,
        the last computed layout is taken.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    matplotlib.Axes object
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    if layout is None: layout = adata.add['draw_graph_layout'][-1]
    if 'X_draw_graph_' + layout not in adata.smp_keys():
        raise ValueError('Did not find {} in adata.smp. Did you compute layout {}?'
                         .format('draw_graph_' + layout, layout))
    axs = scatter(
        adata,
        basis='draw_graph_' + layout,
        color=color,
        groups=groups,
        components=components,
        projection='2d',
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        show=show,
        save=save,
        ax=ax)
    return axs


def tsne(
        adata,
        color=None,
        groups=None,
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter in tSNE basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    matplotlib.Axes object
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    axs = scatter(
        adata,
        basis='tsne',
        color=color,
        groups=groups,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        show=show,
        save=save,
        ax=ax)
    return axs


# ------------------------------------------------------------------------------
# Subgroup identification and ordering - clustering, pseudotime, branching
# and tree inference tools
# ------------------------------------------------------------------------------

def aga(
        adata,
        root=0,       # aga_graph
        fontsize=8,   # aga_graph
        basis='tsne',
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='on data',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        left_margin=0.05,
        layout_graph=None,
        minimal_realized_attachedness=None,
        attachedness_type='relative',
        show=None,
        save=None):
    """Summary figure for approximate graph abstraction.

    See `sc.pl.aga_scatter` and `sc.pl.aga_graph` for the parameters.

    Also see `sc.pl.aga_path` for more possibilities.
    """
    _, axs = pl.subplots(figsize=(8, 4), ncols=2)
    pl.subplots_adjust(left=left_margin, bottom=0.05)
    aga_scatter(adata,
                color='aga_groups',
                basis=basis,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                ax=axs[0],
                show=False)
    axs[1].set_frame_on(False)
    aga_graph(adata, root=root, fontsize=fontsize, ax=axs[1],
              layout=layout_graph,
              attachedness_type=attachedness_type,
              minimal_realized_attachedness=minimal_realized_attachedness,
              show=False)
    utils.savefig_or_show('aga', show=show, save=save)


def aga_scatter(
        adata,
        basis='tsne',
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """Scatter plot of aga groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    matplotlib.Axes object
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    if color is not None:
        if not isinstance(color, list): color = color.split(',')
    else:
        # no need to add pseudotime, is usually not helpful in satter
        color = ['aga_groups']
    if 'aga_groups_original' in adata.add:
        if 'aga_groups' in color:
            color[color.index('aga_groups')] = adata.add['aga_groups_original']
        else:
            color += [adata.add['aga_groups_original']]
    ax = scatter(adata,
                 basis=basis,
                 color=color,
                 groups=groups,
                 components=components,
                 projection=projection,
                 legend_loc=legend_loc,
                 legend_fontsize=legend_fontsize,
                 color_map=color_map,
                 palette=palette,
                 right_margin=right_margin,
                 size=size,
                 title=title,
                 ax=ax,
                 show=False)
    utils.savefig_or_show('aga_' + basis, show=show, save=save)
    return ax


def aga_attachedness(
        adata,
        attachedness_type='scaled',
        color_map=None,
        show=None,
        save=None):
    """Attachedness of aga groups.
    """
    if attachedness_type == 'scaled':
        attachedness = adata.add['aga_attachedness']
    elif attachedness_type == 'distance':
        attachedness = adata.add['aga_distances']
    elif attachedness_type == 'absolute':
        attachedness = adata.add['aga_attachedness_absolute']
    else:
        raise ValueError('Unkown attachedness_type {}.'.format(attachedness_type))
    adjacency = adata.add['aga_adjacency']
    matrix(attachedness, color_map=color_map, show=False)
    for i in range(adjacency.shape[0]):
        neighbors = adjacency[i].nonzero()[1]
        pl.scatter([i for j in neighbors], neighbors, color='green')
    utils.savefig_or_show('aga_attachedness', show=show, save=save)
    # as a stripplot
    if False:
        pl.figure()
        for i, ds in enumerate(attachedness):
            ds = np.log1p(ds)
            x = [i for j, d in enumerate(ds) if i != j]
            y = [d for j, d in enumerate(ds) if i != j]
            pl.scatter(x, y, color='gray')
            neighbors = adjacency[i]
            pl.scatter([i for j in neighbors],
                       ds[neighbors], color='green')
        pl.show()


def aga_graph(
        adata,
        root=0,
        layout=None,
        colors=None,
        groups=None,
        fontsize=None,
        node_size=1,
        node_size_power=0.5,
        edge_width=1,
        title=None,
        ext='png',
        add_noise_to_node_positions=None,
        left_margin=0.01,
        attachedness_type='relative',
        minimal_realized_attachedness=None,
        force_labels_to_front=False,
        show=None,
        save=None,
        ax=None):
    """Plot the abstracted graph.

    Parameters
    ----------
    attachedness_type : {'relative', 'absolute', 'full'}
        For 'full', use the fully connected graph weighted with the
        attachedness matrix.
    layout : {'simple', 'rt', 'rt_circular', 'circle', ...}
        Plotting layout. 'rt' stands for Reingold Tilford and uses
        the igraph layout function.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    if isinstance(colors, list) and isinstance(colors[0], dict): colors = [colors]
    if colors is None or isinstance(colors, str): colors = [colors]
    if isinstance(groups, list) and isinstance(groups[0], str): groups = [groups]
    if groups is None or isinstance(groups, dict) or isinstance(groups, str): groups = [groups]
    if len(colors) != len(groups):
        print(colors, groups)
        raise ValueError('`colors` and `groups` lists need to have the same length.')
    if title is None or isinstance(title, str): title = [title for name in groups]
    if ax is None:
        # 3.72 is the default figure_width obtained in utils.scatter_base
        # for a single panel when rcParams['figure.figsize'][0] = 4
        figure_width = rcParams['figure.figsize'][0] * len(colors)
        top = 0.93
        fig, axs = pl.subplots(ncols=len(colors),
                               figsize=(figure_width, rcParams['figure.figsize'][1]),
                               subplotpars=sppars(left=left_margin, bottom=0,
                                                  right=0.99, top=top))
    else:
        axs = ax
    if len(colors) == 1: axs = [axs]
    for icolor, color in enumerate(colors):
        _aga_graph_single(
            adata,
            layout=layout,
            root=root,
            colors=color,
            groups=groups[icolor],
            fontsize=fontsize,
            node_size=node_size,
            node_size_power=node_size_power,
            edge_width=edge_width,
            attachedness_type=attachedness_type,
            minimal_realized_attachedness=minimal_realized_attachedness,
            ext=ext,
            ax=axs[icolor],
            title=title[icolor],
            add_noise_to_node_positions=add_noise_to_node_positions,
            force_labels_to_front=force_labels_to_front)
    if ext == 'pdf':
        logg.warn('Be aware that saving as pdf exagerates thin lines.')
    utils.savefig_or_show('aga_graph', show=show, ext=ext, save=save)
    return axs if ax is None else None


def _aga_graph_single(
        adata,
        root=0,
        colors=None,
        groups=None,
        fontsize=None,
        node_size=1,
        node_size_power=0.5,
        edge_width=1,
        title=None,
        ext='pdf',
        ax=None,
        layout=None,
        add_noise_to_node_positions=None,
        minimal_realized_attachedness=None,
        attachedness_type=False,
        draw_edge_labels=False,
        force_labels_to_front=False):
    avail_attachedness_types = {'relative', 'absolute', 'full'}
    if attachedness_type not in avail_attachedness_types:
        raise ValueError('Pass available `attachedness_type`, one of {}.'
                         .format(avail_attachedness_types))
    if colors is None and 'aga_groups_colors_original' in adata.add:
        colors = adata.add['aga_groups_colors_original']
    if groups is None and 'aga_groups_names_original' in adata.add:
        groups = adata.add['aga_groups_names_original']
    elif groups in adata.smp_keys():
        groups = adata.add[groups + '_names']
    elif groups is None:
        groups = adata.add['aga_groups_names']
    if isinstance(root, str) and root in groups:
        root = list(groups).index(root)

    # plot the tree
    if isinstance(adata, nx.Graph):
        G = adata
        colors = ['grey' for n in enumerate(G)]
    else:
        if colors is None:
            if ('aga_groups_colors' not in adata.add
                or len(adata.add['aga_groups_names']) != len(adata.add['aga_groups_colors'])):
                utils.add_colors_for_categorical_sample_annotation(adata, 'aga_groups')
            colors = adata.add['aga_groups_colors']
        for iname, name in enumerate(adata.add['aga_groups_names']):
            if name in sett._ignore_categories: colors[iname] = 'grey'
        nx_g = nx.Graph(adata.add['aga_adjacency'])
    # node positions
    if attachedness_type in {'relative', 'absolute'}:
        if layout is None: layout = 'simple'
        if layout == 'simple':
            pos = utils.hierarchy_pos(nx_g, root)
            if len(pos) < nx_g.number_of_nodes():
                raise ValueError('This is a forest and not a single tree. '
                                 'Try another `layout`, e.g.,  {"fr"}.')
        else:
            from .. import utils as sc_utils
            g = sc_utils.get_igraph_from_adjacency(adata.add['aga_adjacency'])
            if 'rt' in layout:
                pos_list = g.layout(layout, root=[root]).coords
            else:
                pos_list = g.layout(layout).coords
            pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
        pos_array = np.array([pos[n] for count, n in enumerate(nx_g)])
        pos_y_scale = np.max(pos_array[:, 1]) - np.min(pos_array[:, 1])
        if add_noise_to_node_positions:
            np.random.seed(0)
            pos = {n: pos[n] + 0.025*pos_y_scale*2*(np.random.random()-0.5)
                   for n in pos.keys()}
    elif attachedness_type == 'full':
        from .. import utils as sc_utils
        if layout is None: layout = 'fr'
        g = sc_utils.get_igraph_from_adjacency(adata.add['aga_attachedness'])
        if 'rt' in layout:
            pos_list = g.layout(layout, root=[root]).coords
        else:
            pos_list = g.layout(layout).coords
        pos = {n: [p[0], -p[1]] for n, p in enumerate(pos_list)}
    if len(pos) == 1: pos[0] = 0.5, 0.5
    if ax is None:
        fig = pl.figure()
        ax = pl.axes([0.08, 0.08, 0.9, 0.9], frameon=False)
    # edge widths
    from ..tools import aga
    minimal_realized_attachedness = aga.MINIMAL_REALIZED_ATTACHEDNESS if minimal_realized_attachedness is None else minimal_realized_attachedness
    base_edge_width = edge_width * 1.5*rcParams['lines.linewidth']
    if 'aga_attachedness' in adata.add:
        if attachedness_type == 'relative':
            nx_g = nx.Graph(adata.add['aga_attachedness'])
        else:
            nx_g = nx.Graph(adata.add['aga_attachedness_absolute'])
        widths = [base_edge_width*x[-1]['weight'] for x in nx_g.edges(data=True)]
        if attachedness_type in {'relative', 'absolute'}:
            nx.draw_networkx_edges(nx_g, pos, ax=ax, width=widths, edge_color='grey',
                                   style='dashed', alpha=0.5)
            if attachedness_type == 'relative':
                nx_g = nx.Graph(adata.add['aga_adjacency'])
            else:
                nx_g = nx.Graph(adata.add['aga_adjacency_absolute'])
            if minimal_realized_attachedness == aga.MINIMAL_REALIZED_ATTACHEDNESS:
                widths = [base_edge_width*x[-1]['weight'] for x in nx_g.edges(data=True)]
            else:
                widths = [base_edge_width*(x[-1]['weight'] if x[-1]['weight'] != aga.MINIMAL_REALIZED_ATTACHEDNESS else minimal_realized_attachedness)
                          for x in nx_g.edges(data=True)]
            nx.draw_networkx_edges(nx_g, pos, ax=ax, width=widths, edge_color='black')
        else:
            nx.draw_networkx_edges(nx_g, pos, ax=ax, width=widths, edge_color='black')
    else:
        widths = [base_edge_width for x in nx_g.edges()]
        nx.draw_networkx_edges(nx_g, pos, ax=ax, width=widths, edge_color='black')
    # labels
    if draw_edge_labels:
        edge_labels = {}
        for n1, n2, label in nx_g.edges(data=True):
            edge_labels[(n1, n2)] = '{:.3f}'.format(1. / (1 + 10*label['weight']))
        nx.draw_networkx_edge_labels(nx_g, pos, edge_labels=edge_labels, ax=ax, font_size=5)
    trans = ax.transData.transform
    bbox = ax.get_position().get_points()
    ax_x_min = bbox[0, 0]
    ax_x_max = bbox[1, 0]
    ax_y_min = bbox[0, 1]
    ax_y_max = bbox[1, 1]
    ax_len_x = ax_x_max - ax_x_min
    ax_len_y = ax_y_max - ax_y_min
    trans2 = ax.transAxes.inverted().transform
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])
    base_pie_size = 1/(np.sqrt(nx_g.number_of_nodes()) + 10) * node_size
    median_group_size = np.median(adata.add['aga_groups_sizes'])
    for count, n in enumerate(nx_g.nodes_iter()):
        pie_size = base_pie_size
        pie_size *= np.power(adata.add['aga_groups_sizes'][count] / median_group_size,
                             node_size_power)
        xx, yy = trans(pos[n])     # data coordinates
        xa, ya = trans2((xx, yy))  # axis coordinates
        xa = ax_x_min + (xa - pie_size/2) * ax_len_x
        ya = ax_y_min + (ya - pie_size/2) * ax_len_y
        a = pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
        if is_color_like(colors[count]):
            fracs = [100]
            color = [colors[count]]
        elif isinstance(colors[count], dict):
            color = colors[count].keys()
            fracs = [colors[count][c] for c in color]
            if sum(fracs) < 1:
                color = list(color)
                color.append('grey')
                fracs.append(1-sum(fracs))
        else:
            raise ValueError('{} is neither a dict of valid matplotlib colors '
                             'nor a valid matplotlib color.'.format(colors[count]))
        a.pie(fracs, colors=color)
        if not force_labels_to_front and groups is not None:
            a.text(0.5, 0.5, groups[count],
                   verticalalignment='center',
                   horizontalalignment='center',
                   transform=a.transAxes,
                   size=fontsize)
    # TODO: this is a terrible hack, but if we use the solution above, labels
    # get hidden behind pies
    if force_labels_to_front and groups is not None:
        for count, n in enumerate(nx_g.nodes_iter()):
            # all copy and paste from above
            pie_size = base_pie_size
            pie_size *= np.power(adata.add['aga_groups_sizes'][count] / median_group_size,
                                 node_size_power)
            xx, yy = trans(pos[n])     # data coordinates
            xa, ya = trans2((xx, yy))  # axis coordinates
            xa = ax_x_min + (xa - pie_size/2.0000001) * ax_len_x  # make sure a new axis is created
            ya = ax_y_min + (ya - pie_size/2.0000001) * ax_len_y
            a = pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
            a.set_frame_on(False)
            a.set_xticks([])
            a.set_yticks([])
            a.text(0.5, 0.5, groups[count],
                   verticalalignment='center',
                   horizontalalignment='center',
                   transform=a.transAxes, size=fontsize)
    if title is not None: ax.set_title(title)
    return ax


def aga_path(
        adata,
        nodes=[0],
        keys=[0],
        as_heatmap=True,
        color_map=None,
        xlim=[None, None],
        n_avg=1,
        title=None,
        left_margin=None,
        show_left_y_ticks=None,
        ytick_fontsize=None,
        show_nodes_twin=True,
        legend_fontsize=None,
        save=None,
        show=None,
        ax=None):
    """Gene expression changes along paths in the abstracted graph.
    """
    ax_was_none = ax is None
    if show_left_y_ticks is None:
        show_left_y_ticks = False if show_nodes_twin else True

    orig_node_names = []
    if ('aga_groups_names_original' in adata.add
        and adata.add['aga_groups_original'] != 'louvain_groups'):
        orig_node_names = adata.add['aga_groups_names_original']
    else:
        logg.m('did not find field "aga_groups_names_original" in adata.add, '
               'using aga_group integer ids instead', v=4)

    def moving_average(a, n=n_avg):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    ax = pl.gca() if ax is None else ax
    from matplotlib import transforms
    trans = transforms.blended_transform_factory(
        ax.transData, ax.transAxes)
    if as_heatmap:
        X = []
    x_tick_locs = []
    x_tick_labels = []
    for ikey, key in enumerate(keys):
        x = []
        for igroup, group in enumerate(nodes):
            if ikey == 0: x_tick_locs.append(len(x))
            idcs = np.arange(adata.n_smps)[adata.smp['aga_groups'] == str(group)]
            idcs_group = np.argsort(adata.smp['aga_pseudotime'][adata.smp['aga_groups'] == str(group)])
            idcs = idcs[idcs_group]
            if key in adata.smp_keys(): x += list(adata.smp[key][idcs])
            else: x += list(adata[:, key].X[idcs])
        if n_avg > 1:
            old_len_x = len(x)
            x = moving_average(x)
            if ikey == 0: x_tick_locs = len(x)/old_len_x * np.array(x_tick_locs)
        if not as_heatmap:
            ax.plot(x[xlim[0]:xlim[1]], label=key)
        else:
            X.append(x)
        if ikey == 0:
            for igroup, group in enumerate(nodes):
                if len(orig_node_names) > 0 and group not in orig_node_names:
                    label = orig_node_names[int(group)]
                else:
                    label = group
                if not isinstance(label, int):
                    pl.text(x_tick_locs[igroup], -0.05*(igroup+1),
                            label, transform=trans)
                else:
                    x_tick_labels.append(label)
    if as_heatmap:
        img = ax.imshow(np.array(X), aspect='auto', interpolation='nearest',
                        cmap=color_map)
        ax.set_yticks(range(len(X)))
        ax.set_yticklabels(keys, fontsize=ytick_fontsize)
        ax.set_frame_on(False)
        pl.colorbar(img, ax=ax)
        left_margin = 0.2 if left_margin is None else left_margin
        pl.subplots_adjust(left=left_margin)
    else:
        left_margin = 0.4 if left_margin is None else left_margin
        pl.legend(frameon=False, loc='center left',
                  bbox_to_anchor=(-left_margin, 0.5),
                  fontsize=legend_fontsize)
    ax.set_xticks(x_tick_locs)
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel(adata.add['aga_groups_original'] if ('aga_groups_original' in adata.add
                  and adata.add['aga_groups_original'] != 'louvain_groups')
                  else 'aga groups')
    if show_left_y_ticks:
        utils.pimp_axis(pl.gca().get_yaxis())
        pl.ylabel('as indicated on legend')
    elif not as_heatmap:
        pl.yticks([])
        pl.ylabel('as indicated on legend (a.u.)')
    if show_nodes_twin and not as_heatmap:
        pl.twinx()
        x = []
        for g in nodes:
            x += list(adata.smp['aga_groups'][adata.smp['aga_groups'] == str(g)].astype(int))
        if n_avg > 1: x = moving_average(x)
        pl.plot(x[xlim[0]:xlim[1]], '--', color='black')
        label = 'aga groups' + (' / original groups' if len(orig_node_names) > 0 else '')
        pl.ylabel(label)
        utils.pimp_axis(pl.gca().get_yaxis())
    if title is not None: pl.title(title)
    if show is None and not ax_was_none: show = False
    else: show = sett.autoshow if show is None else show
    utils.savefig_or_show('aga_path', show=show, save=save)
    return ax if ax_was_none else None


def dpt(
        adata,
        basis='diffmap',
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None):
    """Plot results of DPT analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'draw_graph_...'}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: None)
        Sample/ cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {'2d', '3d'}, optional (default: '2d')
         Projection of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    show_tree : bool, optional (default: False)
         This shows the inferred tree. For more than a single branching, the
         result is pretty unreliable. Use tool `aga` (Approximate Graph
         Abstraction) instead.
    """
    dpt_scatter(
        adata,
        basis=basis,
        color=color,
        groups=groups,
        components=components,
        projection=projection,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=save)
    colors = ['dpt_pseudotime']
    if len(np.unique(adata.smp['dpt_groups'])) > 1: colors += ['dpt_groups']
    if color is not None:
        if not isinstance(color, list): colors = color.split(',')
        else: colors = color
    dpt_groups_pseudotime(adata, color_map=color_map, show=False, save=save)
    dpt_timeseries(adata, color_map=color_map, show=show, save=save)


def dpt_scatter(
        adata,
        basis='diffmap',
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None):
    """See parameters of sc.pl.dpt().
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    colors = ['dpt_pseudotime']
    if len(np.unique(adata.smp['dpt_groups'])) > 1: colors += ['dpt_groups']
    if color is not None:
        if not isinstance(color, list): colors = color.split(',')
        else: colors = color
    if components == 'all':
        components_list = ['1,2', '1,3', '1,4', '1,5', '2,3', '2,4', '2,5', '3,4', '3,5', '4,5']
    else:
        if components is None:
            components = '1,2' if '2d' in projection else '1,2,3'
        if not isinstance(components, list): components_list = [components]
        else: components_list = components
    for components in components_list:
        axs = scatter(
            adata,
            basis=basis,
            color=colors,
            groups=groups,
            components=components,
            projection=projection,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            color_map=color_map,
            palette=palette,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False)
        writekey = 'dpt_' + basis + '_components' + components.replace(',', '')
        save = False if save is None else save
        if sett.savefigs or save: utils.savefig(writekey)
    utils.savefig_or_show(writekey, show=show, save=False)


def dpt_timeseries(adata, color_map=None, show=None, save=None):
    # only if number of genes is not too high
    if adata.n_vars <= 11:
        # plot time series as gene expression vs time
        timeseries(adata.X[adata.smp['dpt_order_indices']],
                   var_names=adata.var_names,
                   highlightsX=adata.add['dpt_changepoints'],
                   xlim=[0, 1.3*adata.X.shape[0]])
        pl.xlabel('dpt order')
    elif adata.n_vars <= 50:
        # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
        timeseries_as_heatmap(adata.X[adata.smp['dpt_order_indices'], :40],
                              var_names=adata.var_names,
                              highlightsX=adata.add['dpt_changepoints'])
        pl.xlabel('dpt order')
    if adata.n_vars <= 50:
        utils.savefig_or_show('dpt_timeseries', save=save, show=show)


def dpt_groups_pseudotime(adata, color_map=None, palette=None, show=None, save=None):
    """Plot groups and pseudotime."""
    pl.figure()
    pl.subplot(211)
    timeseries_subplot(adata.smp['dpt_groups'],
                       time=adata.smp['dpt_order'],
                       color=adata.smp['dpt_groups'],
                       highlightsX=adata.add['dpt_changepoints'],
                       ylabel='dpt groups',
                       yticks=(np.arange(len(adata.add['dpt_groups_names']), dtype=int)
                                     if len(adata.add['dpt_groups_names']) < 5 else None),
                       palette=palette)
    pl.subplot(212)
    timeseries_subplot(adata.smp['dpt_pseudotime'],
                       time=adata.smp['dpt_order'],
                       color=adata.smp['dpt_pseudotime'],
                       xlabel='dpt order',
                       highlightsX=adata.add['dpt_changepoints'],
                       ylabel='pseudotime',
                       yticks=[0, 1],
                       color_map=color_map)
    utils.savefig_or_show('dpt_groups_pseudotime', save=save, show=show)


def louvain(
        adata,
        basis='tsne',
        color=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Plot results of Louvain clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'draw_graph_...'}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical sample annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {'2d', '3d'}, optional (default: '2d')
         Projection of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    color_map : str (default: `matplotlib.rcParams['image.cmap']`)
         String denoting matplotlib color map.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    from ..utils import check_adata
    adata = check_adata(adata)
    add_color = []
    if color is not None:
        add_color = color if isinstance(color, list) else color.split(',')
    color = ['louvain_groups'] + add_color
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        groups=groups,
        components=components,
        projection=projection,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False)
    utils.savefig_or_show('louvain_' + basis, show=show, save=save)


def rank_genes_groups(adata, groups=None, n_genes=20, fontsize=8, show=None, save=None):
    """Plot ranking of genes.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    groups : str or list of str
        The groups for which to show the gene ranking.
    n_genes : int, optional (default: 20)
        Number of genes to show.
    fontsize : int, optional (default: 8)
        Fontsize for gene names.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    groups_key = adata.add['rank_genes_groups']
    group_names = adata.add['rank_genes_groups_names'] if groups is None else groups
    # one panel for each group
    n_panels = len(group_names)
    # set up the figure
    if n_panels <= 5:
        n_panels_y = 1
        n_panels_x = n_panels
    else:
        n_panels_y = 2
        n_panels_x = int(n_panels/2+0.5)
    from matplotlib import gridspec
    fig = pl.figure(figsize=(n_panels_x * rcParams['figure.figsize'][0],
                             n_panels_y * rcParams['figure.figsize'][1]))
    left = 0.2/n_panels_x
    bottom = 0.13/n_panels_y
    gs = gridspec.GridSpec(nrows=n_panels_y,
                           ncols=n_panels_x,
                           left=left,
                           right=1-(n_panels_x-1)*left-0.01/n_panels_x,
                           bottom=bottom,
                           top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
                           wspace=0.18)

    for count, group_name in enumerate(group_names):
        pl.subplot(gs[count])
        gene_names = adata.add['rank_genes_groups_gene_names'][group_name]
        scores = adata.add['rank_genes_groups_gene_scores'][group_name]
        for ig, g in enumerate(gene_names[:n_genes]):
            pl.text(ig, scores[ig], gene_names[ig],
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center', fontsize=fontsize)
        pl.title(group_name)
        if n_panels <= 5 or count >= n_panels_x:
            pl.xlabel('ranking')
        if count == 0 or count == n_panels_x:
            pl.ylabel('mean of z-score w.r.t. to bulk mean')
        ymin = np.min(scores)
        ymax = np.max(scores)
        ymax += 0.3*(ymax-ymin)
        pl.ylim([ymin, ymax])
        pl.xlim(-0.9, ig+1-0.1)
    writekey = 'rank_genes_groups_' + adata.add['rank_genes_groups']
    utils.savefig_or_show(writekey, show=show, save=save)


def rank_genes_groups_violin(adata, groups=None, n_genes=20, show=None, save=None):
    """Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list of str
        List of valid gene names.
    groups : list of str
        List of valid group names.
    n_genes : int
        Number of genes to show.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    from ..tools import rank_genes_groups
    groups_key = adata.add['rank_genes_groups']
    group_names = adata.add['rank_genes_groups_names'] if groups is None else groups
    group_loop = (group_name for group_name in group_names)
    check_is_computed = True
    for group_name in group_loop:
        keys = []
        gene_names = []
        gene_loop = (gene_item for gene_item in enumerate(adata.add['rank_genes_groups_gene_names'][group_name][:n_genes]))
        for gene_counter, gene_name in gene_loop:
            identifier = rank_genes_groups._build_identifier(
                groups_key, group_name, gene_counter, gene_name)
            if check_is_computed and identifier not in set(adata.smp_keys()):
                raise ValueError('You need to set `compute_distribution=True` in `sc.tl.rank_genes_groups()` if you want to use this visualiztion. '
                                 'You might consider simply using `sc.pl.rank_genes_groups_means()` and `sc.pl.violin(adata_raw, gene_name, group_by=grouping)` instead.')
                check_is_computed = False
            keys.append(identifier)
            gene_names.append(gene_name)
        ax = violin(adata, keys, show=False)
        ax.set_title(group_name)
        ax.set_ylabel('z-score w.r.t. to bulk mean')
        ax.set_xticklabels(gene_names, rotation='vertical')
        writekey = 'rank_genes_groups_' + adata.add['rank_genes_groups'] + '_' + group_name
        utils.savefig_or_show(writekey, show=show, save=save)


def sim(adata, shuffle=False, show=None, save=None):
    """Plot results of simulation.

    Parameters
    ----------
    show : bool, optional (default: None)
         Show the plot.
    shuffle : bool, optional (default: False)
         Shuffle the data.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    from .. import utils as sc_utils
    X = adata.X
    genenames = adata.var_names
    tmax = adata.add['tmax_write']
    n_real = X.shape[0]/tmax
    if not shuffle:
        timeseries(X,
                   var_names=genenames,
                   xlim=[0, 1.25*X.shape[0]],
                   highlightsX=np.arange(tmax, n_real * tmax, tmax),
                   xlabel='realizations / time steps')
        utils.savefig_or_show('sim', save=save, show=show)
    else:
        # shuffled data
        X, rows = sc_utils.subsample(X, seed=1)
        timeseries(X,
                   var_names=genenames,
                   xlim=[0, 1.25*X.shape[0]],
                   highlightsX=np.arange(tmax, n_real * tmax, tmax),
                   xlabel='index (arbitrary order)')
        utils.savefig_or_show('sim_shuffled', save=save, show=show)
