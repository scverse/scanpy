# Author: F. Alex Wolf (http://falexwolf.de)
"""Plotting

Plotting functions for each tool and toplevel plotting functions for AnnData.
"""

import warnings
import numpy as np
import networkx as nx
from ..compat.matplotlib import pyplot as pl
from matplotlib.colors import is_color_like
# general functions
from .toplevel import scatter, violin
from .toplevel import matrix
from .toplevel import timeseries, timeseries_subplot, timeseries_as_heatmap
from .toplevel import ranking, ranking_deprecated
from .toplevel import savefig, savefig_or_show
# preprocessing
from .preprocessing import filter_genes_dispersion
from . import utils
from .. import settings as sett

utils.init_plotting_params()


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
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    """
    show = params['show'] if 'show' in params else None
    if 'show' in params: del params['show']
    pca_scatter(adata, **params, show=False)
    pca_ranking(adata, show)


def pca_scatter(
        adata,
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """See parameters of pl.pca().
    """
    from ..examples import check_adata
    adata = check_adata(adata, verbosity=-1)
    axs = scatter(
        adata,
        basis='pca',
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    savefig_or_show('pca_scatter', show=show)
    return axs


def pca_ranking(adata, comps=None, show=None):
    """Rank genes according to contributions to PCs.
    """
    if isinstance(comps, str): comps = comps.split(',')
    keys = ['PC1', 'PC2', 'PC3'] if comps is None else ['PC{}'.format(c) for c in comps]
    ranking(adata, 'var', keys)
    if sett.savefigs: savefig('pca_ranking_components')
    ranking(adata, 'add', 'pca_variance_ratio', labels='PC')
    savefig_or_show('pca_ranking_variance')


def diffmap(
        adata,
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str or list, optional (default: '1,2')
         String of the form '1,2' or 'all' or list. First component is 1 or '1'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    if comps == 'all':
        comps_list = ['{},{}'.format(*((i, i+1) if i % 2 == 1 else (i+1, i)))
                      for i in range(1, adata.smp['X_diffmap'].shape[1])]
    else:
        if comps is None: comps = '1,2' if '2d' in layout else '1,2,3'
        if not isinstance(comps, list): comps_list = [comps]
        else: comps_list = comps
    for comps in comps_list:
        axs = scatter(
            adata,
            basis='diffmap',
            color=color,
            names=names,
            comps=comps,
            cont=cont,
            layout=layout,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            cmap=cmap,
            pal=pal,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False)
        writekey = 'diffmap'
        if isinstance(comps, list): comps = ','.join([str(comp) for comp in comps])
        writekey += '_comps' + comps.replace(',', '')
        if sett.savefigs: savefig(writekey)
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()
    return axs


def draw_graph(
        adata,
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Scatter in tSNE basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".

    Returns
    -------
    matplotlib.Axes object
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    axs = scatter(
        adata,
        basis='draw_graph',
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    savefig_or_show('draw_graph', show=show)
    return axs


def tsne(
        adata,
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Scatter in tSNE basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".

    Returns
    -------
    matplotlib.Axes object
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    axs = scatter(
        adata,
        basis='tsne',
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    savefig_or_show('tsne', show=show)
    return axs


def spring(
        adata,
        color=None,
        names=None,
        comps='1,2',
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Plot spring scatter plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    if 'X_spring' in adata.smp:
        Y = adata.smp['X_spring']
    else:
        raise ValueError('Need to run tool `spring` before plotting.')
    axs = scatter(
        adata,
        basis='spring',
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        # defined in plotting
        title=title,
        show=False)
    savefig_or_show('spring', show=show)


# ------------------------------------------------------------------------------
# Subgroup identification and ordering - clustering, pseudotime, branching
# and tree inference tools
# ------------------------------------------------------------------------------

def aga(
        adata,
        root=0,       # aga_tree
        fontsize=8,   # aga_tree
        basis='tsne',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='on data',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Summary figure for approximate graph abstraction.

    See aga_scatter and aga_tree for most of the arguments.
    """
    _, axs = pl.subplots(figsize=(8, 4), ncols=2)
    aga_scatter(adata,
                color='aga_groups',
                basis=basis,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                ax=axs[0],
                show=False)
    axs[1].set_frame_on(False)
    aga_tree(adata, root=root, fontsize=fontsize, ax=axs[1], show=False)
    show = sett.autoshow if show is None else show
    savefig_or_show('aga', show=show)


def aga_scatter(
        adata,
        basis='tsne',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        ax=None,
        show=None):
    """See parameters of sc.pl.aga().
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    if color is not None:
        if not isinstance(color, list): color = color.split(',')
    else:
        if 'aga_pseudotime' in adata.smp_keys(): color = ['aga_pseudotime']
        else: color = []
        color += ['aga_groups']
    if 'aga_groups_original' in adata.add:
        if 'aga_groups' in color:
            color[color.index('aga_groups')] = adata.add['aga_groups_original']
        else:
            color += [adata.add['aga_groups_original']]
    if comps == 'all':
        comps_list = ['1,2', '1,3', '1,4', '1,5', '2,3', '2,4', '2,5', '3,4',
                      '3,5', '4,5']
    else:
        if comps is None: comps = '1,2' if '2d' in layout else '1,2,3'
        if not isinstance(comps, list): comps_list = [comps]
        else: comps_list = comps
    ax_was_none = ax is None
    for comps in comps_list:
        ax = scatter(adata,
                     basis=basis,
                     color=color,
                     names=names,
                     comps=comps,
                     cont=cont,
                     layout=layout,
                     legend_loc=legend_loc,
                     legend_fontsize=legend_fontsize,
                     cmap=cmap,
                     pal=pal,
                     right_margin=right_margin,
                     size=size,
                     title=title,
                     ax=ax,
                     show=False)
        writekey = 'aga_' + basis + '_comps' + comps.replace(',', '')
        if sett.savefigs: savefig(writekey)
    if show is None and not ax_was_none: show = False
    else: show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()
    return ax

def aga_attachedness(adata):
    attachedness = adata.add['aga_attachedness']
    adjacency = adata.add['aga_adjacency']
    matrix(attachedness, show=False)
    for i in range(adjacency.shape[0]):
        neighbors = adjacency[i].nonzero()[1]
        pl.scatter([i for j in neighbors], neighbors, color='green')
    pl.show()
    # pl.figure()
    # for i, ds in enumerate(attachedness):
    #     ds = np.log1p(ds)
    #     x = [i for j, d in enumerate(ds) if i != j]
    #     y = [d for j, d in enumerate(ds) if i != j]
    #     pl.scatter(x, y, color='gray')
    #     neighbors = adjacency[i]
    #     pl.scatter([i for j in neighbors],
    #                ds[neighbors], color='green')
    # pl.show()


def aga_tree(
        adata,
        root=0,
        colors=None,
        names=None,
        fontsize=None,
        node_size=1,
        ext='pdf',
        ax=None,
        show=None):
    if colors is None or isinstance(colors, str): colors = [colors]
    if isinstance(colors, list) and isinstance(colors[0], dict): colors = [colors]
    if names is None or isinstance(names, str) or isinstance(names, dict): names = [names]
    if len(colors) != len(names):
        raise ValueError('`colors` and `names` lists need to have the same length.')
    from matplotlib import rcParams
    figure_width = 1.3 * rcParams['figure.figsize'][0] * len(colors)
    fig, axs = pl.subplots(ncols=len(colors),
                           figsize=(figure_width, 1.3*rcParams['figure.figsize'][1]))
    if len(colors) == 1: axs = [axs]
    for icolor, color in enumerate(colors):
        show_color = False if icolor != len(colors)-1 else show
        _aga_tree_single_color(
            adata,
            root=root,
            colors=color,
            names=names[icolor],
            fontsize=fontsize,
            node_size=1,
            ext=ext,
            ax=axs[icolor],
            show=show_color)


def _aga_tree_single_color(
        adata,
        root=0,
        colors=None,
        names=None,
        fontsize=None,
        node_size=1,
        ext='pdf',
        ax=None,
        draw_edge_labels=False,
        show=None):
    from .. import logging as logg
    from matplotlib import rcParams
    if colors is None and 'aga_groups_colors_original' in adata.add:
        colors = adata.add['aga_groups_colors_original']
    if names is None and 'aga_groups_names_original' in adata.add:
        names = adata.add['aga_groups_names_original']
    elif names in adata.smp_keys():
        names = adata.add[names + '_names']
    elif names is None:
        names = adata.add['aga_groups_names']
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
        G = nx.Graph(adata.add['aga_adjacency'])
    try:
        pos = utils.hierarchy_pos(G, root)
        pos_array = np.array([pos[n] for count, n in enumerate(G)])
        pos_y_scale = np.max(pos_array[:, 1]) - np.min(pos_array[:, 1])
        np.random.seed(0)
        pos = {n: pos[n] + 0.025*pos_y_scale*2*(np.random.random()-0.5)
               for n in pos.keys()}
    # print(pos)
    # print(0.01*np.random.random((adata.add['aga_adjacency'].shape[0], 2)))
    # pos += 0.01*np.random.random((adata.add['aga_adjacency'].shape[0], 2))
    except Exception:
        pos = nx.spring_layout(G)
        logg.warn('could not draw tree layout, now using fruchterman-reingold layout')
    if len(pos) == 1: pos[0] = 0.5, 0.5
    ax_was_none = False
    if ax is None:
        fig = pl.figure()
        ax = pl.axes([0.08, 0.08, 0.9, 0.9], frameon=False)
        ax_was_none = True
    # edge widths
    G = nx.Graph(adata.add['aga_attachedness'])
    widths = [1.5*rcParams['lines.linewidth']*x[-1]['weight'] for x in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color='grey',
                           style='dashed', alpha=0.7)
    G = nx.Graph(adata.add['aga_adjacency'])
    widths = [1.5*rcParams['lines.linewidth']*x[-1]['weight'] for x in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color='black')
    # labels
    if draw_edge_labels:
        edge_labels = {}
        for n1, n2, label in G.edges(data=True):
            edge_labels[(n1, n2)] = '{:.3f}'.format(1. / (1 + 10*label['weight']))
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax, font_size=5)
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
    base_pie_size = 1/(np.sqrt(G.number_of_nodes()) + 10) * node_size
    median_group_size = np.median(adata.add['aga_groups_sizes'])
    for count, n in enumerate(G):
        pie_size = base_pie_size
        pie_size *= np.sqrt(adata.add['aga_groups_sizes'][count] / median_group_size)
        xx, yy = trans(pos[n])     # data coordinates
        xa, ya = trans2((xx, yy))  # axis coordinates
        xa = ax_x_min + (xa - pie_size/2) * ax_len_x
        ya = ax_y_min + (ya - pie_size/2) * ax_len_y
        a = pl.axes([xa, ya, pie_size * ax_len_x, pie_size * ax_len_y])
        if is_color_like(colors[count]):
            fracs = [100]
            color = [colors[count]]
        else:
            color = colors[count].keys()
            fracs = [colors[count][c] for c in color]
            if sum(fracs) < 1:
                color = list(color)
                color.append('grey')
                fracs.append(1-sum(fracs))
                # names[count] += '\n?'
        a.pie(fracs, colors=color)
        if names is not None:
            a.text(0.5, 0.5, names[count],
                   verticalalignment='center',
                   horizontalalignment='center',
                   transform=a.transAxes, size=fontsize)
    if show is None and not ax_was_none: show = False
    else: show = sett.autoshow if show is None else show
    savefig_or_show('aga_tree', show, ext=ext)
    return ax if ax_was_none else None


def aga_timeseries(
        adata,
        nodes=[0],
        keys=[0],
        xlim=[None, None],
        n_avg=1,
        left_margin=0.4,
        show_left_y_ticks=None,
        show_nodes_twin=True,
        legend_fontsize=None,
        ax=None,
        show=None):
    ax_was_none = ax is None
    if show_left_y_ticks is None:
        show_left_y_ticks = False if show_nodes_twin else True

    orig_node_names = []
    if 'aga_groups_names_original' in adata.add:
        orig_node_names = adata.add['aga_groups_names_original']
    else:
        logg.warn('did not find field "aga_groups_names_original" in adata.add, '
                  'using aga_group integer ids instead')

    def moving_average(a, n=n_avg):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    ax = pl.gca()
    from matplotlib import transforms
    trans = transforms.blended_transform_factory(
        ax.transData, ax.transAxes)
    for ikey, key in enumerate(keys):
        x = []
        for igroup, group in enumerate(nodes):
            if ikey == 0:
                if len(orig_node_names) > 0 and group not in orig_node_names:
                    label = orig_node_names[int(group)]
                else:
                    label = group
                pl.text(len(x), -0.05*(igroup+1), label, transform=trans)
            idcs = np.arange(adata.n_smps)[adata.smp['aga_groups'] == str(group)]
            idcs_group = np.argsort(adata.smp['aga_pseudotime'][adata.smp['aga_groups'] == str(group)])
            idcs = idcs[idcs_group]
            if key in adata.smp_keys(): x += list(adata.smp[key][idcs])
            else: x += list(adata[:, key].X[idcs])
        if n_avg > 1: x = moving_average(x)
        pl.plot(x[xlim[0]:xlim[1]], label=key)
    pl.legend(frameon=False, loc='center left',
              bbox_to_anchor=(-left_margin, 0.5),
              fontsize=legend_fontsize)
    pl.xticks([])
    if show_left_y_ticks:
        utils.pimp_axis(pl.gca().get_yaxis())
        pl.ylabel('as indicated on legend')
    else:
        pl.yticks([])
        pl.ylabel('as indicated on legend (a.u.)')
    if show_nodes_twin:
        pl.twinx()
        x = []
        for g in nodes:
            x += list(adata.smp['aga_groups'][adata.smp['aga_groups'] == str(g)].astype(int))
        if n_avg > 1: x = moving_average(x)
        pl.plot(x[xlim[0]:xlim[1]], '--', color='black')
        label = 'aga groups' + (' / original groups' if len(orig_node_names) > 0 else '')
        pl.ylabel(label)
    if show is None and not ax_was_none: show = False
    else: show = sett.autoshow if show is None else show
    savefig_or_show('aga_timeseries', show)
    return ax if ax_was_none else None


def aga_sc_tree(adata, root, show=None):
    G = nx.Graph(adata.add['aga_adjacency'])
    node_sets = []
    sorted_aga_groups = adata.smp['aga_groups'][adata.smp['aga_order']]
    for n in adata.add['aga_groups_names']:
        node_sets.append(np.flatnonzero(n == sorted_aga_groups))
    sc_G = utils.hierarchy_sc(G, root, node_sets)
    ax = aga_tree(sc_G, root)
    savefig_or_show('aga_sc_tree', show)
    return ax


def dpt(
        adata,
        basis='diffmap',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Plot results of DPT analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: first annotation)
        Sample/ cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String of the form '1,2' or 'all'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    show_tree : bool, optional (default: False)
         This shows the inferred tree. For more than a single branching, the
         result is pretty unreliable. Use tool `aga` (Approximate Graph
         Abstraction) instead.
    """
    dpt_scatter(
        adata,
        basis=basis,
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    colors = ['dpt_pseudotime']
    if len(np.unique(adata.smp['dpt_groups'])) > 1: colors += ['dpt_groups']
    if color is not None:
        if not isinstance(color, list): colors = color.split(',')
        else: colors = color
    dpt_timeseries(adata, cmap=cmap, show=show)


def dpt_scatter(
        adata,
        basis='diffmap',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """See parameters of sc.pl.dpt().
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    colors = ['dpt_pseudotime']
    if len(np.unique(adata.smp['dpt_groups'])) > 1: colors += ['dpt_groups']
    if color is not None:
        if not isinstance(color, list): colors = color.split(',')
        else: colors = color
    if comps == 'all':
        comps_list = ['1,2', '1,3', '1,4', '1,5', '2,3', '2,4', '2,5', '3,4', '3,5', '4,5']
    else:
        if comps is None:
            comps = '1,2' if '2d' in layout else '1,2,3'
        if not isinstance(comps, list): comps_list = [comps]
        else: comps_list = comps
    # adata.add['highlights'] = (
    #    # list([adata.add['iroot']])
    #    [i for g in adata.add['dpt_groupconnects'] for i in g])
    #  + [adata.add['dpt_grouptips'][i][1]
    #     for i in range(len(adata.add['dpt_grouptips']))
    #     if adata.add['dpt_grouptips'][i][1] != -1])
    #  + [adata.add['dpt_grouptips'][i][0]
    #     for i in range(len(adata.add['dpt_grouptips']))
    #  if adata.add['dpt_grouptips'][i][1] != -1])

    #    adata.add['highlights'] = adata.add['dpt_groups_connects'][adata.add['dpt_groups_connects'].nonzero()].A1
    #    adata.add['highlights'] = {i: '{}, {}->{}'.format(
    #        i,
    #        adata.add['dpt_groups_connects'].nonzero()[0][ii],
    #        adata.add['dpt_groups_connects'].nonzero()[1][ii])
    #        for ii, i in enumerate(adata.add['highlights'])}
    for comps in comps_list:
        axs = scatter(adata,
                      basis=basis,
                      color=colors,
                      names=names,
                      comps=comps,
                      cont=cont,
                      layout=layout,
                      legend_loc=legend_loc,
       legend_fontsize=legend_fontsize,
                      cmap=cmap,
                      pal=pal,
                      right_margin=right_margin,
                      size=size,
                      title=title,
                      show=False)
        writekey = 'dpt_' + basis + '_comps' + comps.replace(',', '')
        if sett.savefigs: savefig(writekey)
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()


def dpt_timeseries(adata, cmap=None, show=None):
    # plot segments and pseudotime
    if True:
        dpt_segments_pseudotime(adata, 'viridis' if cmap is None else cmap)
        # time series plot
        # only if number of genes is not too high
        if adata.X.shape[1] <= 11:
            # plot time series as gene expression vs time
            timeseries(adata.X[adata.smp['dpt_order']],
                             varnames=adata.var_names,
                             highlightsX=adata.add['dpt_changepoints'],
                             xlim=[0, 1.3*adata.X.shape[0]])
            pl.xlabel('dpt order')
            if sett.savefigs: savefig('dpt_vsorder')
        elif adata.X.shape[1] < 50:
            # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
            timeseries_as_heatmap(adata.X[adata.smp['dpt_order'], :40],
                                        varnames=adata.var_names,
                                        highlightsX=adata.add['dpt_changepoints'])
            pl.xlabel('dpt order')
            if sett.savefigs: savefig('dpt_heatmap')
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()


def dpt_segments_pseudotime(adata, cmap=None, pal=None):
    """Plot segments and pseudotime."""
    pl.figure()
    pl.subplot(211)
    timeseries_subplot(adata.smp['dpt_groups'][adata.smp['dpt_order'], np.newaxis],
                             c=adata.smp['dpt_groups'][adata.smp['dpt_order']],
                             highlightsX=adata.add['dpt_changepoints'],
                             ylabel='dpt groups',
                             yticks=(np.arange(len(adata.add['dpt_groups_names']), dtype=int)
                                     if len(adata.add['dpt_groups_names']) < 5 else None),
                             pal=pal)
    pl.subplot(212)
    timeseries_subplot(adata.smp['dpt_pseudotime'][adata.smp['dpt_order'], np.newaxis],
                             c=adata.smp['dpt_pseudotime'][adata.smp['dpt_order']],
                             xlabel='dpt order',
                             highlightsX=adata.add['dpt_changepoints'],
                             ylabel='pseudotime',
                             yticks=[0, 1],
                             cmap=cmap)
    if sett.savefigs: savefig('dpt_segpt')


def dbscan(
        adata,
        basis='tsne',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Plot results of DBSCAN clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    colors = ['dbscan_groups']
    if color is not None: colors += color.split(',')
    axs = scatter(
        adata,
        basis=basis,
        color=colors,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    savefig_or_show('dbscan_' + basis)


def louvain(
        adata,
        basis='tsne',
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        pal=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Plot results of Louvain clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: None)
         Colors cycle to use for categorical groups.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    title : str, optional (default: None)
         Provide title for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    add_color = []
    if color is not None:
        add_color = color if isinstance(color, list) else color.split(',')
    color = ['louvain_groups'] + add_color
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        names=names,
        comps=comps,
        cont=cont,
        layout=layout,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        cmap=cmap,
        pal=pal,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False)
    savefig_or_show('louvain_' + basis)


def paths(
        adata,
        basis='diffmap',
        dist_threshold=None,
        single_panel=True,
        color=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        cmap=None,
        right_margin=None,
        size=None,
        title=None,
        show=None):
    """Plot paths.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    dist_threshold : float
        Distance threshold to decide what still belongs in the path.
    single_panel : bool (default: True)
        If False, show separate panel for each group.
    color : string or list of strings, optional (default: first annotation)
        Keys for sample/cell annotation either as list or string "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names)
        Restrict to a few categories in categorical sample annotation.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    cmap : str (default: continuous: inferno/ categorical: finite palette)
         String denoting matplotlib color map.
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    names = None if names is None else names.split(',') if isinstance(names, str) else names
    # from ..tools.paths import process_dists_from_paths
    # process_dists_from_paths(adata, dist_threshold)

    color_base = ['paths_groups']
    if color is not None:
        if isinstance(color, list):
            color_base += color
        else:
            color_base += [color]
    adata.add['highlights'] = [adata.add['iroot']]

    # add continuous distance coloring
    if single_panel:
        for iname, name in enumerate(adata.add['paths_groups_names']):
            if names is None or (names is not None and name in names):
                title = 'dist_from_path_' + name
                adata.smp[title] = adata.add['paths_dists_from_paths'][iname]
                # color_base.append(title)
                adata.add['highlights'] += [adata.add['paths_groups_fateidx'][iname]]

        axs = scatter(
            adata,
            basis=basis,
            color=color_base,
            names=names,
            comps=comps,
            cont=cont,
            layout=layout,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            cmap=cmap,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False)
        writekey = 'paths_' + basis + '_' + adata.add['paths_type']
        if sett.savefigs: savefig(writekey)
    else:
        for iname, name in enumerate(adata.add['paths_groups_names']):
            if names is None or (names is not None and name in names):
                title = 'dist_from_path_' + name
                adata.smp[title] = adata.add['paths_dists_from_paths'][iname]
                # color_base.append(title)
                adata.add['highlights'] = ([adata.add['iroot']]
                                       + [adata.add['paths_groups_fateidx'][iname]])
            axs = scatter(
                adata,
                basis=basis,
                color=color_base,
                names=[name],
                comps=comps,
                cont=cont,
                layout=layout,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                cmap=cmap,
                right_margin=right_margin,
                size=size,
                title=title,
                show=False)
            del color_base[-1]
            writekey = 'paths_' + basis + '_' + adata.add['paths_type'] + '_' + name
            if sett.savefigs: savefig(writekey)
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()


def diffrank(adata, n_genes=20, show=None):
    """Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_genes : int
        Number of genes to show.
    """
    ranking_deprecated(adata, toolkey='diffrank', n_genes=n_genes)
    writekey = 'diffrank_' + adata.add['diffrank_groups']
    savefig_or_show(writekey)


def tgdyn_simple(adata, n_genes=10, show=None):
    """Plot analysis of single-cell dynamics on the graph.

    Parameters
    ----------
    adata : dict
        Dictionary returned by get_data.
    n_genes : int
    """
    plot_tgdyn_simple_ranking(adata, n_genes, show)


def tgdyn_simple_ranking(adata, n_genes=10, show=None):
    """Plot ranking.

    TODO
    ----
    Replace with call to plotting.plot_ranking.

    Parameters
    ----------
    dtgdyn : dict
        Dictionary returned by tgdyn.
    adata : dict
        Dictionary returned by get_data.
    """
    n_panels = adata.add['tgdyn_simple_groups_ids_bigenough'].shape[0]

    # find minimum velocity to set y-axis limit
    ymin = 1
    for igroup in adata.add['tgdyn_simple_groups_ids_bigenough']:
        genes = adata.add['tgdyn_simple_genes_sorted'][igroup, :n_genes]
        ymin = np.min([ymin,
                       np.min(np.abs(adata.add['tgdyn_simple_vs_norm'][igroup, genes]))])

    # determine n_panels in x and y direction
    if n_panels <= 5:
        n_panels_y = 1
        n_panels_x = n_panels
    else:
        n_panels_y = 2
        n_panels_x = int(n_panels/2+0.5)

    # do the actual plotting
    fig = pl.figure(figsize=(n_panels_x*4, n_panels_y*4))
    pl.subplots_adjust(left=0.17/n_panels_x, right=0.99, bottom=0.13)
    count = 1
    for igroup in adata.add['tgdyn_simple_groups_ids_bigenough']:
        # generate the panel
        fig.add_subplot(n_panels_y, n_panels_x, count)
        # get the velocity to plot
        v = adata.add['tgdyn_simple_vs_norm'][igroup]
        # loop over the top-ranked genes
        for ig, g in enumerate(adata.add['tgdyn_simple_genes_sorted'][igroup, :n_genes]):
            marker = r'\leftarrow' if v[g] < 0 else r'\rightarrow'
            color = 'red' if v[g] < 0 else 'green'
            pl.text(ig,
                    np.abs(v[g]),
                    r'$ ' + marker + '$ ' + adata.var_names[g],
                    color=color,
                    rotation='vertical',
                    verticalalignment='bottom',
                    horizontalalignment='center',
                    fontsize=8)
        title = adata.add['tgdyn_simple_groups'] + ' ' + str(adata.add['tgdyn_simple_groups_names'][igroup])
        pl.title(title)
        pl.xlim(-0.9, ig+1-0.1)
        pl.ylim(-0.02+ymin, 1.15)
        if n_panels <= 5 or count > n_panels_x:
            pl.xlabel('ranking')
        if count == 1 or count == n_panels_x + 1:
            pl.ylabel('|velocity$_{gene}$|/max$_{genes}$|velocity$_{gene}$|')
        else:
            pl.yticks([])
        count += 1

    writekey = 'tgdyn_simple_' + adata.add['tgdyn_simple_groups']
    savefig_or_show(writekey)


def sim(adata, params=None, show=None):
    """Plot results of simulation.
    """
    from .. import utils as sc_utils
    X = adata.X
    genenames = adata.var_names
    tmax = adata.add['tmax_write']
    n_real = X.shape[0]/tmax
    timeseries(X,
               varnames=genenames,
               xlim=[0, 1.25*X.shape[0]],
               highlightsX=np.arange(tmax, n_real * tmax, tmax),
               xlabel='realizations / time steps')
    if sett.savefigs: savefig('sim')
    # shuffled data
    X, rows = sc_utils.subsample(X, seed=1)
    timeseries(X,
               varnames=genenames,
               xlim=[0, 1.25*X.shape[0]],
               highlightsX=np.arange(tmax, n_real * tmax, tmax),
               xlabel='index (arbitrary order)')
    savefig_or_show('sim_shuffled')
