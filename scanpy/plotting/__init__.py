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
    savefig_or_show('tsne')
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


def ega_tree(adata, root=0, colors=None, names=None, show=None, fontsize=None):
    # plot the tree
    if isinstance(adata, nx.Graph):
        G = adata
        colors = ['grey' for n in enumerate(G)]
    else:
        if colors is None:
            if ('ega_groups_colors' not in adata.add
                or len(adata.add['ega_groups_names']) != len(adata.add['ega_groups_colors'])):
                utils.add_colors_for_categorical_sample_annotation(adata, 'ega_groups')
            colors = adata.add['ega_groups_colors']
        else: colors = colors
        if names is None:
            names = {i: n for i, n in enumerate(adata.add['ega_groups_names'])}
        for iname, name in enumerate(adata.add['ega_groups_names']):
            if name in sett._ignore_categories: colors[iname] = 'grey'
        G = nx.Graph(adata.add['ega_groups_adjacency'])
    pos = utils.hierarchy_pos(G, root)
    # pos = nx.spring_layout(G)
    if len(pos) == 1: pos[0] = 0.5, 0.5
    fig = pl.figure()
    ax = pl.axes([0.08, 0.08, 0.9, 0.9], frameon=False)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edges(G, pos, ax=ax)
    edge_labels = {}
    for n1, n2, label in G.edges(data=True):
        edge_labels[(n1, n2)] = '{:.3f}'.format(label['weight'])
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax, font_size=5)
    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform
    pl.xticks([])
    pl.yticks([])
    piesize = 1/(np.sqrt(G.number_of_nodes()) + 10)
    p2 = piesize/2.0
    for n_cnt, n in enumerate(G):
        xx, yy = trans(pos[n])     # figure coordinates
        xa, ya = trans2((xx, yy))  # normalized coordinates
        a = pl.axes([xa-p2, ya-p2, piesize, piesize])
        if is_color_like(colors[n_cnt]):
            fracs = [100]
            color = [colors[n_cnt]]
        else:
            color = colors[n_cnt].keys()
            fracs = [colors[n_cnt][c] for c in color]
            if sum(fracs) < 1:
                color = list(color)
                color.append('grey')
                fracs.append(1-sum(fracs))
                # names[n_cnt] += '\n?'
        a.pie(fracs, colors=color)
        if names is not None:
            a.text(0.5, 0.5, names[n_cnt],
                   verticalalignment='center',
                   horizontalalignment='center',
                   transform=a.transAxes, size=fontsize)
    savefig_or_show('ega_tree', show)
    return ax


def ega_sc_tree(adata, root, show=None):
    G = nx.Graph(adata.add['ega_groups_adjacency'])
    node_sets = []
    sorted_ega_groups = adata.smp['ega_groups'][adata.smp['ega_order']]
    for n in adata.add['ega_groups_names']:
        node_sets.append(np.flatnonzero(n == sorted_ega_groups))
    # print(node_sets)
    sc_G = utils.hierarchy_sc(G, root, node_sets)
    ax = ega_tree(sc_G, root)
    savefig_or_show('ega_sc_tree', show)
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
        show_tree=False,
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
         result is pretty unreliable. Use tool `ega` (Extremal Graph
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
    if 'dpt_groups' in colors and show_tree: dpt_tree(adata, show=False)
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


def dpt_tree(adata, root=0, colors=None, names=None, show=None, fontsize=None):
    # plot the tree
    if isinstance(adata, nx.Graph):
        G = adata
        colors = ['grey' for n in enumerate(G)]
    else:
        if colors is None:
            if ('dpt_groups_colors' not in adata.add
                or len(adata.add['dpt_groups_names']) != len(adata.add['dpt_groups_colors'])):
                utils.add_colors_for_categorical_sample_annotation(adata, 'dpt_groups')
            colors = adata.add['dpt_groups_colors']
        else: colors = colors
        if names is None:
            names = {i: n for i, n in enumerate(adata.add['dpt_groups_names'])}
        for iname, name in enumerate(adata.add['dpt_groups_names']):
            if name in sett._ignore_categories: colors[iname] = 'grey'
        G = nx.Graph(adata.add['dpt_groups_adjacency'])
    pos = utils.hierarchy_pos(G, root)
    # pos = nx.spring_layout(G)
    if len(pos) == 1: pos[0] = 0.5, 0.5
    fig = pl.figure()
    ax = pl.axes([0.08, 0.08, 0.9, 0.9], frameon=False)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edges(G, pos, ax=ax)  #, edge_labels=labels)
    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform
    pl.xticks([])
    pl.yticks([])
    piesize = 1/(np.sqrt(G.number_of_nodes()) + 10)
    p2 = piesize/2.0
    for n_cnt, n in enumerate(G):
        xx, yy = trans(pos[n])     # figure coordinates
        xa, ya = trans2((xx, yy))  # normalized coordinates
        a = pl.axes([xa-p2, ya-p2, piesize, piesize])
        if is_color_like(colors[n_cnt]):
            fracs = [100]
            color = [colors[n_cnt]]
        else:
            color = colors[n_cnt].keys()
            fracs = [colors[n_cnt][c] for c in color]
            if sum(fracs) < 1:
                color = list(color)
                color.append('grey')
                fracs.append(1-sum(fracs))
                names[n_cnt] += '\n?'
        a.pie(fracs, colors=color)
        if names is not None:
            a.text(0.5, 0.5, names[n_cnt],
                   verticalalignment='center',
                   horizontalalignment='center',
                   transform=a.transAxes, size=fontsize)
    savefig_or_show('dpt_tree', show)
    return ax


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
