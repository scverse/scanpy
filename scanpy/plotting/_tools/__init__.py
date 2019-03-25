import numpy as np
import pandas as pd
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams, cm, colors
from anndata import AnnData
from typing import Union, Optional

from .. import _utils as utils
from ...utils import doc_params, sanitize_anndata
from ... import logging as logg
from .._anndata import scatter, ranking
from .._utils import timeseries, timeseries_subplot, timeseries_as_heatmap
from .._docs import doc_scatter_bulk, doc_show_save_ax
from .scatterplots import pca, plot_scatter
from matplotlib.colors import Colormap

# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------


@doc_params(scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def pca_overview(adata, **params):
    """\
    Plot PCA results.

    The parameters are the ones of the scatter plot. Call pca_ranking separately
    if you want to change the default settings.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: `None`)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {scatter_bulk}
    show : bool, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    """
    show = params['show'] if 'show' in params else None
    if 'show' in params: del params['show']
    scatterplots.pca(adata, **params, show=False)
    pca_loadings(adata, show=False)
    pca_variance_ratio(adata, show=show)


# backwards compat
pca_scatter = pca


def pca_loadings(adata, components=None, show=None, save=None):
    """Rank genes according to contributions to PCs.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    components : str or list of integers, optional
        For example, ``'1,2,3'`` means ``[1, 2, 3]``, first, second, third
        principal component.
    show : bool, optional (default: `None`)
        Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    """
    if components is None: components = [1, 2, 3]
    elif isinstance(components, str): components = components.split(',')
    components = np.array(components) - 1
    ranking(adata, 'varm', 'PCs', indices=components)
    utils.savefig_or_show('pca_loadings', show=show, save=save)


def pca_variance_ratio(adata, n_pcs=30, log=False, show=None, save=None):
    """Plot the variance ratio.

    Parameters
    ----------
    n_pcs : `int`, optional (default: `30`)
         Number of PCs to show.
    log : `bool`, optional (default: `False`)
         Plot on logarithmic scale..
    show : `bool`, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    """
    ranking(adata, 'uns', 'variance_ratio', n_points=n_pcs, dictionary='pca', labels='PC', log=log)
    utils.savefig_or_show('pca_variance_ratio', show=show, save=save)


# ------------------------------------------------------------------------------
# Subgroup identification and ordering - clustering, pseudotime, branching
# and tree inference tools
# ------------------------------------------------------------------------------


def dpt_timeseries(adata, color_map=None, show=None, save=None, as_heatmap=True):
    """Heatmap of pseudotime series.

    Parameters
    ----------
    as_heatmap : bool (default: False)
        Plot the timeseries as heatmap.
    """
    if adata.n_vars > 100:
        logg.warn('Plotting more than 100 genes might take some while,'
                  'consider selecting only highly variable genes, for example.')
    # only if number of genes is not too high
    if as_heatmap:
        # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
        timeseries_as_heatmap(adata.X[adata.obs['dpt_order_indices'].values],
                              var_names=adata.var_names,
                              highlightsX=adata.uns['dpt_changepoints'],
                              color_map=color_map)
    else:
        # plot time series as gene expression vs time
        timeseries(adata.X[adata.obs['dpt_order_indices'].values],
                   var_names=adata.var_names,
                   highlightsX=adata.uns['dpt_changepoints'],
                   xlim=[0, 1.3*adata.X.shape[0]])
    pl.xlabel('dpt order')
    utils.savefig_or_show('dpt_timeseries', save=save, show=show)


def dpt_groups_pseudotime(adata, color_map=None, palette=None, show=None, save=None):
    """Plot groups and pseudotime."""
    pl.figure()
    pl.subplot(211)
    timeseries_subplot(adata.obs['dpt_groups'].cat.codes,
                       time=adata.obs['dpt_order'].values,
                       color=np.asarray(adata.obs['dpt_groups']),
                       highlightsX=adata.uns['dpt_changepoints'],
                       ylabel='dpt groups',
                       yticks=(np.arange(len(adata.obs['dpt_groups'].cat.categories), dtype=int)
                                     if len(adata.obs['dpt_groups'].cat.categories) < 5 else None),
                       palette=palette)
    pl.subplot(212)
    timeseries_subplot(adata.obs['dpt_pseudotime'].values,
                       time=adata.obs['dpt_order'].values,
                       color=adata.obs['dpt_pseudotime'].values,
                       xlabel='dpt order',
                       highlightsX=adata.uns['dpt_changepoints'],
                       ylabel='pseudotime',
                       yticks=[0, 1],
                       color_map=color_map)
    utils.savefig_or_show('dpt_groups_pseudotime', save=save, show=show)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups(adata, groups=None, n_genes=20, gene_symbols=None, key=None, fontsize=8,
                      ncols=4, sharey=True, show=None, save=None, ax=None, **kwds):
    """\
    Plot ranking of genes.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    gene_symbols : `str`
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names`.
    n_genes : `int`, optional (default: 20)
        Number of genes to show.
    fontsize : `int`, optional (default: 8)
        Fontsize for gene names.
    ncols : `int`, optional (default: 4)
        Number of panels shown per row.
    sharey: `bool`, optional (default: True)
        Controls if the y-axis of each panels should be shared. But passing
        `sharey=False`, each panel has its own y-axis range.
    {show_save_ax}
    """
    if 'n_panels_per_row' in kwds:  n_panels_per_row  = kwds['n_panels_per_row']
    else: n_panels_per_row = ncols
    if key is None: key = 'rank_genes_groups'
    reference = str(adata.uns[key]['params']['reference'])
    group_names = (adata.uns[key]['names'].dtype.names
                   if groups is None else groups)
    # one panel for each group
    # set up the figure
    n_panels_x = min(n_panels_per_row, len(group_names))
    n_panels_y = np.ceil(len(group_names) / n_panels_x).astype(int)

    from matplotlib import gridspec
    fig = pl.figure(figsize=(n_panels_x * rcParams['figure.figsize'][0],
                             n_panels_y * rcParams['figure.figsize'][1]))
    gs = gridspec.GridSpec(nrows=n_panels_y,
                           ncols=n_panels_x,
                           wspace=0.22,
                           hspace=0.3)

    ax0 = None
    ymin = np.Inf
    ymax = -np.Inf
    for count, group_name in enumerate(group_names):
        if sharey is True:
            if ax0 is None:
                ax = fig.add_subplot(gs[count])
                ax0 = ax
            else:
                ax = fig.add_subplot(gs[count], sharey=ax0)
        else:
            ax = fig.add_subplot(gs[count])

        gene_names = adata.uns[key]['names'][group_name]
        scores = adata.uns[key]['scores'][group_name]
        for ig, g in enumerate(gene_names[:n_genes]):
            gene_name = gene_names[ig]
            if adata.raw is not None and adata.uns[key]['params']['use_raw']:
                ax.text(
                    ig, scores[ig],
                    gene_name if gene_symbols is None else adata.raw.var[gene_symbols][gene_name],
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center', fontsize=fontsize)
            else:
                ax.text(
                    ig, scores[ig],
                    gene_name if gene_symbols is None else adata.var[gene_symbols][gene_name],
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center', fontsize=fontsize)
        ax.set_title('{} vs. {}'.format(group_name, reference))
        if count >= n_panels_x * (n_panels_y - 1):
            ax.set_xlabel('ranking')

        # print the 'score' label only on the first panel per row.
        if count % n_panels_x == 0:
            ax.set_ylabel('score')

        ax.set_xlim(-0.9, ig + 1-0.1)

        if sharey is True:
            ymin = min(ymin, np.min(scores))
            ymax = max(ymax, np.max(scores))
        else:
            ymin = np.min(scores)
            ymax = np.max(scores)
            ymax += 0.3*(np.max(scores)-np.min(scores))
            ax.set_ylim(ymin, ymax)

    if sharey is True:
        ymax += 0.3*(ymax-ymin)
        ax.set_ylim(ymin, ymax)

    writekey = ('rank_genes_groups_'
                + str(adata.uns[key]['params']['groupby']))
    utils.savefig_or_show(writekey, show=show, save=save)


@doc_params(show_save_ax=doc_show_save_ax)
def _rank_genes_groups_plot(adata, plot_type='heatmap', groups=None,
                            n_genes=10, groupby=None, key=None,
                            show=None, save=None, **kwds):
    """\
    Plot ranking of genes using the specified plot type

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.
    {show_save_ax}
    """
    if key is None:
        key = 'rank_genes_groups'

    if 'dendrogram' not in kwds:
        kwds['dendrogram'] = True
    if groupby is None:
        groupby = str(adata.uns[key]['params']['groupby'])
    group_names = (adata.uns[key]['names'].dtype.names
                   if groups is None else groups)

    gene_names = []
    start = 0
    group_positions = []
    group_names_valid = []
    for group in group_names:
        # get all genes that are 'not-nan'
        genes_list = [gene for gene in adata.uns[key]['names'][group] if not pd.isnull(gene)][:n_genes]
        if len(genes_list) == 0:
            logg.warn("No genes found for group {}".format(group))
            continue
        gene_names.extend(genes_list)
        end = start + len(genes_list)
        group_positions.append((start, end -1))
        group_names_valid.append(group)
        start = end

    group_names = group_names_valid
    if plot_type == 'dotplot':
        from .._anndata import dotplot
        dotplot(adata, gene_names, groupby, var_group_labels=group_names,
                var_group_positions=group_positions, show=show, save=save, **kwds)

    elif plot_type == 'heatmap':
        from .._anndata import heatmap
        heatmap(adata, gene_names, groupby, var_group_labels=group_names,
                var_group_positions=group_positions, show=show, save=save, **kwds)

    elif plot_type == 'stacked_violin':
        from .._anndata import stacked_violin
        return stacked_violin(adata, gene_names, groupby, var_group_labels=group_names,
                       var_group_positions=group_positions, show=show, save=save, **kwds)

    elif plot_type == 'tracksplot':
        from .._anndata import tracksplot
        return tracksplot(adata, gene_names, groupby, var_group_labels=group_names,
                       var_group_positions=group_positions, show=show, save=save, **kwds)

    elif plot_type == 'matrixplot':
        from .._anndata import matrixplot
        matrixplot(adata, gene_names, groupby, var_group_labels=group_names,
                   var_group_positions=group_positions, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_heatmap(adata, groups=None, n_genes=10, groupby=None, key=None,
                              show=None, save=None, **kwds):
    """\
    Plot ranking of genes using heatmap plot (see `scanpy.api.pl.heatmap`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see `scanpy.api.pl.heatmap`).
    key : `str`
        Key used to store the ranking results in `adata.uns`.
    **kwds : keyword arguments
        Are passed to `scanpy.api.pl.heatmap`.
    {show_save_ax}
    """

    _rank_genes_groups_plot(adata, plot_type='heatmap', groups=groups, n_genes=n_genes,
                            groupby=groupby, key=key, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_tracksplot(adata, groups=None, n_genes=10, groupby=None, key=None,
                              show=None, save=None, **kwds):
    """\
    Plot ranking of genes using heatmap plot (see `scanpy.api.pl.heatmap`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see `scanpy.api.pl.heatmap`).
    key : `str`
        Key used to store the ranking results in `adata.uns`.
    **kwds : keyword arguments
        Are passed to `scanpy.api.pl.tracksplot`.
    {show_save_ax}
    """

    _rank_genes_groups_plot(adata, plot_type='tracksplot', groups=groups, n_genes=n_genes,
                            groupby=groupby, key=key, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_dotplot(adata, groups=None, n_genes=10, groupby=None, key=None,
                              show=None, save=None, **kwds):
    """\
    Plot ranking of genes using dotplot plot (see `scanpy.api.pl.dotplot`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see `scanpy.api.pl.dotplot`).
    key : `str`
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `scanpy.api.pl.dotplot`.
    """

    _rank_genes_groups_plot(adata, plot_type='dotplot', groups=groups, n_genes=n_genes,
                            groupby=groupby, key=key, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_stacked_violin(adata, groups=None, n_genes=10, groupby=None, key=None,
                                     show=None, save=None, **kwds):
    """\
    Plot ranking of genes using stacked_violin plot (see `scanpy.api.pl.stacked_violin`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see `scanpy.api.pl.stacked_violin`).
    key : `str`
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `scanpy.api.pl.stacked_violin`.
    """

    _rank_genes_groups_plot(adata, plot_type='stacked_violin', groups=groups, n_genes=n_genes,
                            groupby=groupby, key=key, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_matrixplot(adata, groups=None, n_genes=10, groupby=None, key=None,
                                 show=None, save=None, **kwds):
    """\
    Plot ranking of genes using matrixplot plot (see `scanpy.api.pl.matrixplot`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the gene ranking.
    n_genes : `int`, optional (default: 10)
        Number of genes to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see `scanpy.api.pl.matrixplot`).
    key : `str`
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `scanpy.api.pl.matrixplot`.
    """

    _rank_genes_groups_plot(adata, plot_type='matrixplot', groups=groups, n_genes=n_genes,
                            groupby=groupby, key=key, show=show, save=save, **kwds)


@doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_violin(
        adata, groups=None, n_genes=20,
        gene_names=None, gene_symbols=None,
        use_raw=None,
        key=None,
        split=True,
        scale='width',
        strip=True, jitter=True, size=1,
        ax=None, show=None, save=None):
    """\
    Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : list of `str`, optional (default: `None`)
        List of group names.
    n_genes : `int`, optional (default: 20)
        Number of genes to show. Is ignored if `gene_names` is passed.
    gene_names : `None` or list of `str` (default: `None`)
        List of genes to plot. Is only useful if interested in a custom gene list,
        which is not the result of :func:`scanpy.api.tl.rank_genes_groups`.
    gene_symbols : `str`, optional (default: `None`)
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names` displayed in the plot.
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present. Defaults to the value that
        was used in :func:`~scanpy.api.tl.rank_genes_groups`.
    split : `bool`, optional (default: `True`)
        Whether to split the violins or not.
    scale : `str`, optional (default: 'width')
        See `seaborn.violinplot`.
    strip : `bool`, optional (default: `True`)
        Show a strip plot on top of the violin plot.
    jitter : `int`, `float`, `bool`, optional (default: `True`)
        If set to 0, no points are drawn. See `seaborn.stripplot`.
    size : `int`, optional (default: 1)
        Size of the jitter points.
    {show_save_ax}
    """
    if key is None:
        key = 'rank_genes_groups'
    groups_key = str(adata.uns[key]['params']['groupby'])
    if use_raw is None:
        use_raw = bool(adata.uns[key]['params']['use_raw'])
    reference = str(adata.uns[key]['params']['reference'])
    groups_names = (adata.uns[key]['names'].dtype.names
                    if groups is None else groups)
    if isinstance(groups_names, str): groups_names = [groups_names]
    axs = []
    for group_name in groups_names:
        if gene_names is None:
            gene_names = adata.uns[
                key]['names'][group_name][:n_genes]
        df = pd.DataFrame()
        new_gene_names = []
        for g in gene_names:
            if adata.raw is not None and use_raw:
                X_col = adata.raw[:, g].X
                if gene_symbols:
                    g = adata.raw.var[gene_symbols][g]
            else:
                X_col = adata[:, g].X
                if gene_symbols:
                    g = adata.var[gene_symbols][g]
            if issparse(X_col): X_col = X_col.toarray().flatten()
            new_gene_names.append(g)
            df[g] = X_col
        df['hue'] = adata.obs[groups_key].astype(str).values
        if reference == 'rest':
            df.loc[df['hue'] != group_name, 'hue'] = 'rest'
        else:
            df.loc[~df['hue'].isin([group_name, reference]), 'hue'] = np.nan
        df['hue'] = df['hue'].astype('category')
        df_tidy = pd.melt(df, id_vars='hue', value_vars=new_gene_names)
        x = 'variable'
        y = 'value'
        hue_order = [group_name, reference]
        import seaborn as sns
        _ax = sns.violinplot(x=x, y=y, data=df_tidy, inner=None,
                             hue_order=hue_order, hue='hue', split=split,
                             scale=scale, orient='vertical', ax=ax)
        if strip:
            _ax = sns.stripplot(x=x, y=y, data=df_tidy,
                                hue='hue', dodge=True, hue_order=hue_order,
                                jitter=jitter, color='black', size=size, ax=_ax)
        _ax.set_xlabel('genes')
        _ax.set_title('{} vs. {}'.format(group_name, reference))
        _ax.legend_.remove()
        _ax.set_ylabel('expression')
        _ax.set_xticklabels(new_gene_names, rotation='vertical')
        writekey = ('rank_genes_groups_'
                    + str(adata.uns[key]['params']['groupby'])
                    + '_' + group_name)
        utils.savefig_or_show(writekey, show=show, save=save)
        axs.append(_ax)
    if show == False: return axs


def sim(adata, tmax_realization=None, as_heatmap=False, shuffle=False,
        show=None, save=None):
    """Plot results of simulation.

    Parameters
    ----------
    as_heatmap : bool (default: False)
        Plot the timeseries as heatmap.
    tmax_realization : int or None (default: False)
        Number of observations in one realization of the time series. The data matrix
        adata.X consists in concatenated realizations.
    shuffle : bool, optional (default: False)
        Shuffle the data.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    show : bool, optional (default: `None`)
        Show the plot, do not return axis.
    """
    from ... import utils as sc_utils
    if tmax_realization is not None: tmax = tmax_realization
    elif 'tmax_write' in adata.uns: tmax = adata.uns['tmax_write']
    else: tmax = adata.n_obs
    n_realizations = adata.n_obs/tmax
    if not shuffle:
        if not as_heatmap:
            timeseries(adata.X,
                       var_names=adata.var_names,
                       xlim=[0, 1.25*adata.n_obs],
                       highlightsX=np.arange(tmax, n_realizations*tmax, tmax),
                       xlabel='realizations')
        else:
            # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
            timeseries_as_heatmap(adata.X,
                                  var_names=adata.var_names,
                                  highlightsX=np.arange(tmax, n_realizations*tmax, tmax))
        pl.xticks(np.arange(0, n_realizations*tmax, tmax),
                  np.arange(n_realizations).astype(int) + 1)
        utils.savefig_or_show('sim', save=save, show=show)
    else:
        # shuffled data
        X = adata.X
        X, rows = sc_utils.subsample(X, seed=1)
        timeseries(X,
                   var_names=adata.var_names,
                   xlim=[0, 1.25*adata.n_obs],
                   highlightsX=np.arange(tmax, n_realizations*tmax, tmax),
                   xlabel='index (arbitrary order)')
        utils.savefig_or_show('sim_shuffled', save=save, show=show)


@doc_params(show_save_ax=doc_show_save_ax)
def embedding_density(
    adata: AnnData,
    basis: str,
    key: str,
    *,
    group: Optional[str] = None,
    color_map: Union[Colormap, str] = 'YlOrRd',
    bg_dotsize: Optional[int] = 80,
    fg_dotsize:  Optional[int] = 180,
    vmax:  Optional[int] = 1,
    vmin:  Optional[int] = 0,
    save: Union[bool, str, None] = None,
    **kwargs
):
    """Plot the density of cells in an embedding (per condition)

    Plots the gaussian kernel density estimates (over condition) from the
    `sc.tl.embedding_density()` output.

    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    basis
        The embedding over which the density was calculated. This embedded
        representation should be found in `adata.obsm['X_[basis]']``.
    key
        Name of the `.obs` covariate that contains the density estimates
    group
        The category in the categorical observation annotation to be plotted.
        For example, 'G1' in the cell cycle 'phase' covariate.
    color_map
        Matplolib color map to use for density plotting.
    bg_dotsize
        Dot size for background data points not in the `group`.
    fg_dotsize
        Dot size for foreground data points in the `group`.
    vmax
        Density that corresponds to color bar maximum.
    vmin
        Density that corresponds to color bar minimum.
    {show_save_ax}

    Examples
    --------
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.umap(adata)
    >>> sc.tl.embedding_density(adata, basis='umap', groupby='phase')
    >>> sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase', 
    ...                         group='G1')
    >>> sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase', 
    ...                         group='S')
    """
    sanitize_anndata(adata)
    
    # Test user inputs
    basis = basis.lower()

    if basis == 'fa':
        basis = 'draw_graph_fa'

    if 'X_'+basis not in adata.obsm_keys():
        raise ValueError('Cannot find the embedded representation `adata.obsm[X_{!r}]`. '
                         'Compute the embedding first.'.format(basis))

    if key not in adata.obs:
        raise ValueError('Please run `sc.tl.embedding_density()` first and specify the correct key.')

    if key+'_params' not in adata.uns:
        raise ValueError('Please run `sc.tl.embedding_density()` first and specify the correct key.')

    if 'components' in kwargs:
        logg.warn('Components were specified, but will be ignored. Only the '
                  'components used to calculate the density can be plotted.')
        del kwargs['components']

    components = adata.uns[key+'_params']['components']
    groupby = adata.uns[key+'_params']['covariate']

    if (group is None) and (groupby is not None):
        raise ValueError('Densities were calculated over an `.obs` covariate. '
                         'Please specify a group from this covariate to plot.')

    if (group is not None) and (group not in adata.obs[groupby].cat.categories):
        raise ValueError('Please specify a group from the `.obs` category over which the density '
                         'was calculated.')

    if (np.min(adata.obs[key]) < 0) or (np.max(adata.obs[key]) > 1):
        raise ValueError('Densities should be scaled between 0 and 1.')
    
    # Define plotting data
    dens_values = -np.ones(adata.n_obs)
    dot_sizes = np.ones(adata.n_obs)*bg_dotsize

    if group is not None:
        group_mask = (adata.obs[groupby] == group)
        dens_values[group_mask] = adata.obs[key][group_mask]
        dot_sizes[group_mask] = np.ones(sum(group_mask))*fg_dotsize

    else:
        dens_values = adata.obs[key]
        dot_sizes = np.ones(adata.n_obs)*fg_dotsize

    # Make the color map
    if isinstance(color_map, str):
        cmap = cm.get_cmap(color_map)
    else:
        cmap = color_map

    #norm = colors.Normalize(vmin=-1, vmax=1)
    adata_vis = adata.copy()
    adata_vis.obs['Density'] = dens_values

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap.set_over('black')
    cmap.set_under('lightgray')
    
    # Ensure title is blank as default
    if 'title' not in kwargs:
        title=""
    else:
        title = kwargs.pop('title')

    # Plot the graph
    return plot_scatter(adata_vis, basis, components=components, color='Density',
                        color_map=cmap, norm=norm, size=dot_sizes, vmax=vmax,
                        vmin=vmin, save=save, title=title, **kwargs)
