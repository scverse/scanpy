# Author: F. Alex Wolf (http://falexwolf.de)
"""Plotting

Plotting functions for each tool and toplevel plotting functions for AnnData.
"""

import numpy as np
from ..compat.matplotlib import pyplot as pl
from .toplevel import scatter, violin
from .toplevel import timeseries, timeseries_subplot, timeseries_as_heatmap
from .toplevel import ranking
from .toplevel import savefig, savefig_or_show
from . import utils
from .. import sett

utils.init_plotting_params()


def pca(adata,
        smp=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legendloc='right margin',
        cmap=None,
        pal=None,
        right_margin=None,
        size=3,
        titles=None,
        show=None):
    """
    Plot PCA scatter.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    smps = scatter(adata,
                   basis='pca',
                   color=smp,
                   names=names,
                   comps=comps,
                   cont=cont,
                   layout=layout,
                   legendloc=legendloc,
                   cmap=cmap,
                   pal=pal,
                   right_margin=right_margin,
                   size=size,
                   titles=titles)
    writekey = sett.basekey + '_pca' + sett.plotsuffix
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(writekey)
    elif show: pl.show()


def diffmap(adata,
            smp=None,
            names=None,
            comps=None,
            cont=None,
            layout='2d',
            legendloc='right margin',
            cmap=None,
            pal=None,
            right_margin=None,
            size=3,
            titles=None,
            show=None):
    """Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str or list, optional (default: '1,2')
         String of the form '1,2' or 'all' or list. First component is 1 or '1'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    if comps == 'all':
        comps_list = ['1,2', '1,3', '1,4', '1,5', '2,3', '2,4', '2,5', '3,4', '3,5', '4,5']
    else:
        if comps is None:
            comps = '1,2' if '2d' in layout else '1,2,3'
        comps_list = [comps]
    for comps in comps_list:
        scatter(adata,
                basis='diffmap',
                color=smp,
                names=names,
                comps=comps,
                cont=cont,
                layout=layout,
                legendloc=legendloc,
                cmap=cmap,
                pal=pal,
                right_margin=right_margin,
                size=size,
                titles=titles)
        writekey = sett.basekey + '_diffmap'
        if isinstance(comps, list):
            comps = ','.join([str(comp) for comp in comps])
        writekey += sett.plotsuffix + '_comps' + comps.replace(',', '')
        if sett.savefigs: savefig(writekey)
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()


def tsne(adata,
         smp=None,
         names=None,
         comps=None,
         cont=None,
         layout='2d',
         legendloc='right margin',
         cmap=None,
         pal=None,
         right_margin=None,
         size=3,
         titles=None,
         show=None):
    """Scatter in tSNE basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    scatter(adata,
            basis='tsne',
            color=smp,
            names=names,
            comps=comps,
            cont=cont,
            layout=layout,
            legendloc=legendloc,
            cmap=cmap,
            pal=pal,
            right_margin=right_margin,
            size=size,
            titles=titles)
    writekey = sett.basekey + '_tsne' + sett.plotsuffix
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(writekey)
    elif show: pl.show()


def spring(adata,
           smp=None,
           names=None,
           comps='1,2',
           cont=None,
           layout='2d',
           legendloc='right margin',
           cmap=None,
           pal=None,
           right_margin=None,
           size=3,
           show=None):
    """Plot spring scatter plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: 0.2)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    Y = adata.smp['X_spring']
    if True:
        smps = scatter(adata,
                       basis='spring',
                       color=smp,
                       names=names,
                       comps=comps,
                       cont=cont,
                       layout=layout,
                       legendloc=legendloc,
                       cmap=cmap,
                       pal=pal,
                       right_margin=right_margin,
                       size=size,
                       # defined in plotting
                       titles=['Fruchterman-Reingold step: 12'])
        writekey = sett.basekey + '_spring'
        writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
        show = sett.autoshow if show is None else show
        if sett.savefigs: savefig(writekey)
        elif show: pl.show()
    else:
        Adj = dspring['Adj']
        istep = dspring['istep']
        # TODO: don't save the adjacency matrix!!!
        import scanpy as sc
        sc.write(dspring['writekey']+'_step{:02}'.format(istep), dspring)
        # compute the next steps
        istep_init = istep + 1
        add_steps = params['add_steps']
        del params['add_steps']
        for istep in istep_init + np.arange(add_steps, dtype=int):
            sett.mt(0, 'compute Fruchterman-Reingold layout: step', istep)
            Y = fruchterman_reingold_layout(Adj, Yinit=Y, iterations=step_size)
            sett.mt(0, 'finished computation')
            _plot({'Y': Y}, adata, istep, **params)
        # save state of Y to outfile
        dspring['Y'] = Y
        dspring['istep'] = istep
        sc.write(dspring['writekey'], dspring)


def dpt(adata,
        basis='diffmap',
        smp=None,
        names=None,
        comps=None,
        cont=None,
        layout='2d',
        legendloc='right margin',
        cmap=None,
        pal=None,
        right_margin=None,
        size=3,
        show=None):
    """
    Plot results of DPT analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    smp : str, optional (default: first annotation)
        Sample/ cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    comps : str, optional (default: '1,2')
         String of the form '1,2' or 'all'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    # scatter plot
    smps = ['dpt_pseudotime']
    if len(np.unique(adata.smp['dpt_groups'])) > 1:
        smps += ['dpt_groups']
    adata.add['highlights'] = (list([adata.add['iroot']])   # also plot the tip cell indices
                               + [adata.add['dpt_segtips'][i][1] for i in range(len(adata.add['dpt_segtips']))
                               if adata.add['dpt_segtips'][i][1] != -1])
    if smp is not None:
        smps += smp.split(',')
    if comps == 'all':
        comps_list = ['1,2', '1,3', '1,4', '1,5', '2,3', '2,4', '2,5', '3,4', '3,5', '4,5']
    else:
        if comps is None:
            comps = '1,2' if '2d' in layout else '1,2,3'
        comps_list = [comps]
    for comps in comps_list:
        smps = scatter(adata,
                       basis=basis,
                       color=smps,
                       names=names,
                       comps=comps,
                       cont=cont,
                       layout=layout,
                       legendloc=legendloc,
                       cmap=cmap,
                       pal=pal,
                       right_margin=right_margin,
                       size=size)
        writekey = sett.basekey + '_dpt_' + basis
        writekey += (sett.plotsuffix + '_comps' + comps.replace(',', ''))
        if sett.savefigs: savefig(writekey)
    # plot segments and pseudotime
    dpt_segments_pseudotime(adata, 'viridis' if cmap is None else cmap)
    # time series plot
    # only if number of genes is not too high
    X = adata.X
    writekey = sett.basekey + '_' + 'dpt' + sett.plotsuffix
    if X.shape[1] <= 11:
        # plot time series as gene expression vs time
        timeseries(X[adata.smp['dpt_order']],
                         varnames=adata.var_names,
                         highlightsX=adata.add['dpt_changepoints'],
                         xlim=[0, 1.3*X.shape[0]])
        pl.xlabel('dpt order')
        if sett.savefigs: savefig(writekey + '_vsorder')
    elif X.shape[1] < 50:
        # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
        timeseries_as_heatmap(X[adata.smp['dpt_order'], :40],
                                    varnames=adata.var_names,
                                    highlightsX=adata.add['dpt_changepoints'])
        pl.xlabel('dpt order')
        if sett.savefigs: savefig(writekey + '_heatmap')
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
    writekey = sett.basekey + '_' + 'dpt' + sett.plotsuffix
    if sett.savefigs: savefig(writekey + '_segpt')


def dbscan(adata,
           basis='tsne',
           smp=None,
           names=None,
           comps=None,
           cont=None,
           layout='2d',
           legendloc='right margin',
           cmap=None,
           pal=None,
           right_margin=None,
           size=3,
           titles=None,
           show=None):
    """
    Plot results of DBSCAN clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : {'diffmap', 'pca', 'tsne', 'spring'}
        Choose the basis in which to plot.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    smps = ['dbscan_groups']
    if smp is not None:
        smps += smp.split(',')
    smps = scatter(adata,
                   basis=basis,
                   color=smps,
                   names=names,
                   comps=comps,
                   cont=cont,
                   layout=layout,
                   legendloc=legendloc,
                   cmap=cmap,
                   pal=pal,
                   right_margin=right_margin,
                   size=size,
                   titles=titles)
    writekey = sett.basekey + '_dbscan_' + basis
    writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(writekey)
    elif show: pl.show()


def paths(adata,
          basis='diffmap',
          dist_threshold=None,
          single_panel=True,
          smp=None,
          names=None,
          comps=None,
          cont=None,
          layout='2d',
          legendloc='right margin',
          cmap=None,
          right_margin=None,
          size=3,
          titles=None,
          show=None):
    """Plot paths

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    dist_threshold : float
        Distance threshold to decide what still belongs in the path.
    single_panel : bool (default: True)
        If False, show separate panel for each group.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: continuous: inferno/ categorical: finite palette)
         String denoting matplotlib color map.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    names = None if names is None else names.split(',') if isinstance(names, str) else names
    from ..tools.paths import process_dists_from_paths
    process_dists_from_paths(adata, dist_threshold)

    smp_base = ['paths_groups']
    if smp is not None:
        smp_base += [smp]
    adata.add['highlights'] = [adata.add['iroot']]

    # add continuous distance coloring
    if single_panel:
        for iname, name in enumerate(adata.add['paths_groups_names']):
            if names is None or (names is not None and name in names):
                title = 'dist_from_path_' + name
                adata.smp[title] = adata.add['paths_dists_from_paths'][iname]
                smp_base.append(title)
                adata.add['highlights'] += [adata.add['paths_groups_fateidx'][iname]]

        smps = scatter(adata,
                       basis=basis,
                       color=smp_base,
                       names=names,
                       comps=comps,
                       cont=cont,
                       layout=layout,
                       legendloc=legendloc,
                       cmap=cmap,
                       right_margin=right_margin,
                       size=size,
                       titles=titles)
        writekey = sett.basekey + '_paths_' + basis + '_' + adata.add['paths_type']
        writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
        if sett.savefigs: savefig(writekey)
    else:
        for iname, name in enumerate(adata.add['paths_groups_names']):
            if names is None or (names is not None and name in names):
                title = 'dist_from_path_' + name
                adata.smp[title] = adata.add['paths_dists_from_paths'][iname]
                smp_base.append(title)
                adata.add['highlights'] = ([adata.add['iroot']]
                                       + [adata.add['paths_groups_fateidx'][iname]])
            smps = scatter(adata,
                           basis=basis,
                           color=smp_base,
                           names=[name],
                           comps=comps,
                           cont=cont,
                           layout=layout,
                           legendloc=legendloc,
                           cmap=cmap,
                           right_margin=right_margin,
                           size=size,
                           titles=titles)
            del smp_base[-1]
            writekey = sett.basekey + '_paths_' + basis + '_' + adata.add['paths_type']
            writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + '_' + name + sett.plotsuffix
            if sett.savefigs: savefig(writekey)
    show = sett.autoshow if show is None else show
    if not sett.savefigs and show: pl.show()


def diffrank(adata, n_genes=20, show=None):
    """
    Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_genes : int
        Number of genes to show.
    """
    ranking(adata, toolkey='diffrank', n_genes=n_genes)
    writekey = sett.basekey + '_diffrank_' + adata.add['diffrank_groups'] + sett.plotsuffix
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(writekey)
    elif show: pl.show()


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

    writekey = sett.basekey + '_tgdyn_simple_' + adata.add['tgdyn_simple_groups'] + sett.plotsuffix
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(writekey + '_ranking')
    elif show: pl.show()


def sim(adata, params=None, show=None):
    """
    Plot results of simulation.
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
    if sett.savefigs: savefig(sett.basekey + '_sim')
    # shuffled data
    X, rows = sc_utils.subsample(X, seed=1)
    timeseries(X,
               varnames=genenames,
               xlim=[0, 1.25*X.shape[0]],
               highlightsX=np.arange(tmax, n_real * tmax, tmax),
               xlabel='index (arbitrary order)')
    show = sett.autoshow if show is None else show
    if sett.savefigs: savefig(sett.basekey + '_sim_shuffled')
    elif show: pl.show()
