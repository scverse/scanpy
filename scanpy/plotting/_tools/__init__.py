import collections.abc as cabc
from copy import copy
import numpy as np
import pandas as pd
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib import pyplot as pl
from matplotlib import rcParams, cm
from anndata import AnnData
from typing import Union, Optional, List, Sequence, Iterable, Mapping

from .._utils import savefig_or_show
from ..._utils import _doc_params, sanitize_anndata, subsample
from ..._compat import Literal
from ... import logging as logg
from .._anndata import ranking
from .._utils import timeseries, timeseries_subplot, timeseries_as_heatmap
from ..._settings import settings
from .._docs import (
    doc_scatter_embedding,
    doc_show_save_ax,
    doc_rank_genes_groups_plot_args,
    doc_rank_genes_groups_values_to_plot,
    doc_vbound_percentile,
    doc_panels,
)
from ...get import rank_genes_groups_df
from .scatterplots import pca, embedding, _panel_grid
from matplotlib.colors import Colormap
from scanpy.get import obs_df

# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------


@_doc_params(scatter_bulk=doc_scatter_embedding, show_save_ax=doc_show_save_ax)
def pca_overview(adata: AnnData, **params):
    """\
    Plot PCA results.

    The parameters are the ones of the scatter plot. Call pca_ranking separately
    if you want to change the default settings.

    Parameters
    ----------
    adata
        Annotated data matrix.
    color
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw
        Use `raw` attribute of `adata` if present.
    {scatter_bulk}
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.
    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.pl.pca_overview(adata, color="louvain")

    .. currentmodule:: scanpy

    See also
    --------
    tl.pca
    pp.pca
    """
    show = params['show'] if 'show' in params else None
    if 'show' in params:
        del params['show']
    pca(adata, **params, show=False)
    pca_loadings(adata, show=False)
    pca_variance_ratio(adata, show=show)


# backwards compat
pca_scatter = pca


def pca_loadings(
    adata: AnnData,
    components: Union[str, Sequence[int], None] = None,
    include_lowest: bool = True,
    n_points: Union[int, None] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
):
    """\
    Rank genes according to contributions to PCs.

    Parameters
    ----------
    adata
        Annotated data matrix.
    components
        For example, ``'1,2,3'`` means ``[1, 2, 3]``, first, second, third
        principal component.
    include_lowest
        Whether to show the variables with both highest and lowest loadings.
    show
        Show the plot, do not return axis.
    n_points
        Number of variables to plot for each component.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()

    Show first 3 components loadings

    .. plot::
        :context: close-figs

        sc.pl.pca_loadings(adata, components = '1,2,3')


    """
    if components is None:
        components = [1, 2, 3]
    elif isinstance(components, str):
        components = [int(x) for x in components.split(',')]
    components = np.array(components) - 1

    if np.any(components < 0):
        raise ValueError("Component indices must be greater than zero.")

    if n_points is None:
        n_points = min(30, adata.n_vars)
    elif adata.n_vars < n_points:
        raise ValueError(
            f"Tried to plot {n_points} variables, but passed anndata only has {adata.n_vars}."
        )

    ranking(
        adata,
        'varm',
        'PCs',
        n_points=n_points,
        indices=components,
        include_lowest=include_lowest,
    )
    savefig_or_show('pca_loadings', show=show, save=save)


def pca_variance_ratio(
    adata: AnnData,
    n_pcs: int = 30,
    log: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
):
    """\
    Plot the variance ratio.

    Parameters
    ----------
    n_pcs
         Number of PCs to show.
    log
         Plot on logarithmic scale..
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    """
    ranking(
        adata,
        'uns',
        'variance_ratio',
        n_points=n_pcs,
        dictionary='pca',
        labels='PC',
        log=log,
    )
    savefig_or_show('pca_variance_ratio', show=show, save=save)


# ------------------------------------------------------------------------------
# Subgroup identification and ordering – clustering, pseudotime, branching
# and tree inference tools
# ------------------------------------------------------------------------------


def dpt_timeseries(
    adata: AnnData,
    color_map: Union[str, Colormap] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    as_heatmap: bool = True,
):
    """\
    Heatmap of pseudotime series.

    Parameters
    ----------
    as_heatmap
        Plot the timeseries as heatmap.
    """
    if adata.n_vars > 100:
        logg.warning(
            'Plotting more than 100 genes might take some while, '
            'consider selecting only highly variable genes, for example.'
        )
    # only if number of genes is not too high
    if as_heatmap:
        # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
        timeseries_as_heatmap(
            adata.X[adata.obs['dpt_order_indices'].values],
            var_names=adata.var_names,
            highlights_x=adata.uns['dpt_changepoints'],
            color_map=color_map,
        )
    else:
        # plot time series as gene expression vs time
        timeseries(
            adata.X[adata.obs['dpt_order_indices'].values],
            var_names=adata.var_names,
            highlights_x=adata.uns['dpt_changepoints'],
            xlim=[0, 1.3 * adata.X.shape[0]],
        )
    pl.xlabel('dpt order')
    savefig_or_show('dpt_timeseries', save=save, show=show)


def dpt_groups_pseudotime(
    adata: AnnData,
    color_map: Union[str, Colormap, None] = None,
    palette: Union[Sequence[str], Cycler, None] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
):
    """Plot groups and pseudotime."""
    _, (ax_grp, ax_ord) = pl.subplots(2, 1)
    timeseries_subplot(
        adata.obs['dpt_groups'].cat.codes,
        time=adata.obs['dpt_order'].values,
        color=np.asarray(adata.obs['dpt_groups']),
        highlights_x=adata.uns['dpt_changepoints'],
        ylabel='dpt groups',
        yticks=(
            np.arange(len(adata.obs['dpt_groups'].cat.categories), dtype=int)
            if len(adata.obs['dpt_groups'].cat.categories) < 5
            else None
        ),
        palette=palette,
        ax=ax_grp,
    )
    timeseries_subplot(
        adata.obs['dpt_pseudotime'].values,
        time=adata.obs['dpt_order'].values,
        color=adata.obs['dpt_pseudotime'].values,
        xlabel='dpt order',
        highlights_x=adata.uns['dpt_changepoints'],
        ylabel='pseudotime',
        yticks=[0, 1],
        color_map=color_map,
        ax=ax_ord,
    )
    savefig_or_show('dpt_groups_pseudotime', save=save, show=show)


@_doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: int = 20,
    gene_symbols: Optional[str] = None,
    key: Optional[str] = 'rank_genes_groups',
    fontsize: int = 8,
    ncols: int = 4,
    sharey: bool = True,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    ax: Optional[Axes] = None,
    **kwds,
):
    """\
    Plot ranking of genes.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the gene ranking.
    gene_symbols
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names`.
    n_genes
        Number of genes to show.
    fontsize
        Fontsize for gene names.
    ncols
        Number of panels shown per row.
    sharey
        Controls if the y-axis of each panels should be shared. But passing
        `sharey=False`, each panel has its own y-axis range.
    {show_save_ax}


    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.rank_genes_groups(adata)


    Plot top 10 genes (default 20 genes)

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups(adata, n_genes=10)

    .. currentmodule:: scanpy

    See also
    --------
    tl.rank_genes_groups

    """
    if 'n_panels_per_row' in kwds:
        n_panels_per_row = kwds['n_panels_per_row']
    else:
        n_panels_per_row = ncols
    if n_genes < 1:
        raise NotImplementedError(
            "Specifying a negative number for n_genes has not been implemented for "
            f"this plot. Received n_genes={n_genes}."
        )

    reference = str(adata.uns[key]['params']['reference'])
    group_names = adata.uns[key]['names'].dtype.names if groups is None else groups
    # one panel for each group
    # set up the figure
    n_panels_x = min(n_panels_per_row, len(group_names))
    n_panels_y = np.ceil(len(group_names) / n_panels_x).astype(int)

    from matplotlib import gridspec

    fig = pl.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0],
            n_panels_y * rcParams['figure.figsize'][1],
        )
    )
    gs = gridspec.GridSpec(nrows=n_panels_y, ncols=n_panels_x, wspace=0.22, hspace=0.3)

    ax0 = None
    ymin = np.Inf
    ymax = -np.Inf
    for count, group_name in enumerate(group_names):
        gene_names = adata.uns[key]['names'][group_name][:n_genes]
        scores = adata.uns[key]['scores'][group_name][:n_genes]

        # Setting up axis, calculating y bounds
        if sharey:
            ymin = min(ymin, np.min(scores))
            ymax = max(ymax, np.max(scores))

            if ax0 is None:
                ax = fig.add_subplot(gs[count])
                ax0 = ax
            else:
                ax = fig.add_subplot(gs[count], sharey=ax0)
        else:
            ymin = np.min(scores)
            ymax = np.max(scores)
            ymax += 0.3 * (ymax - ymin)

            ax = fig.add_subplot(gs[count])
            ax.set_ylim(ymin, ymax)

        ax.set_xlim(-0.9, n_genes - 0.1)

        # Mapping to gene_symbols
        if gene_symbols is not None:
            if adata.raw is not None and adata.uns[key]['params']['use_raw']:
                gene_names = adata.raw.var[gene_symbols][gene_names]
            else:
                gene_names = adata.var[gene_symbols][gene_names]

        # Making labels
        for ig, gene_name in enumerate(gene_names):
            ax.text(
                ig,
                scores[ig],
                gene_name,
                rotation='vertical',
                verticalalignment='bottom',
                horizontalalignment='center',
                fontsize=fontsize,
            )

        ax.set_title('{} vs. {}'.format(group_name, reference))
        if count >= n_panels_x * (n_panels_y - 1):
            ax.set_xlabel('ranking')

        # print the 'score' label only on the first panel per row.
        if count % n_panels_x == 0:
            ax.set_ylabel('score')

    if sharey is True:
        ymax += 0.3 * (ymax - ymin)
        ax.set_ylim(ymin, ymax)

    writekey = f"rank_genes_groups_{adata.uns[key]['params']['groupby']}"
    savefig_or_show(writekey, show=show, save=save)


def _fig_show_save_or_axes(plot_obj, return_fig, show, save):
    """
    Decides what to return
    """
    if return_fig:
        return plot_obj
    else:
        plot_obj.make_figure()
        savefig_or_show(plot_obj.DEFAULT_SAVE_PREFIX, show=show, save=save)
        show = settings.autoshow if show is None else show
        if show is False:
            return plot_obj.get_axes()


def _rank_genes_groups_plot(
    adata: AnnData,
    plot_type: str = 'heatmap',
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    values_to_plot: Optional[str] = None,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    min_logfoldchange: Optional[float] = None,
    key: Optional[str] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    return_fig: Optional[bool] = False,
    gene_symbols: Optional[str] = None,
    **kwds,
):
    """\
    Common function to call the different rank_genes_groups_* plots
    """
    if var_names is not None and n_genes is not None:
        raise ValueError(
            "The arguments n_genes and var_names are mutually exclusive. Please "
            "select only one."
        )

    if var_names is None and n_genes is None:
        # set n_genes = 10 as default when none of the options is given
        n_genes = 10

    if key is None:
        key = 'rank_genes_groups'

    if groupby is None:
        groupby = str(adata.uns[key]['params']['groupby'])
    group_names = adata.uns[key]['names'].dtype.names if groups is None else groups

    if var_names is not None:
        if isinstance(var_names, Mapping):
            # get a single list of all gene names in the dictionary
            var_names_list = sum([list(x) for x in var_names.values()], [])
        elif isinstance(var_names, str):
            var_names_list = [var_names]
        else:
            var_names_list = var_names
    else:
        # dict in which each group is the key and the n_genes are the values
        var_names = {}
        var_names_list = []
        for group in group_names:
            df = rank_genes_groups_df(
                adata,
                group,
                key=key,
                gene_symbols=gene_symbols,
                log2fc_min=min_logfoldchange,
            )

            if gene_symbols is not None:
                df['names'] = df[gene_symbols]

            genes_list = df.names[df.names.notnull()].tolist()

            if len(genes_list) == 0:
                logg.warning(f'No genes found for group {group}')
                continue
            if n_genes < 0:
                genes_list = genes_list[n_genes:]
            else:
                genes_list = genes_list[:n_genes]
            var_names[group] = genes_list
            var_names_list.extend(genes_list)

    # by default add dendrogram to plots
    kwds.setdefault('dendrogram', True)

    if plot_type in ['dotplot', 'matrixplot']:
        # these two types of plots can also
        # show score, logfoldchange and pvalues, in general any value from rank
        # genes groups
        title = None
        values_df = None
        if values_to_plot is not None:
            values_df = _get_values_to_plot(
                adata,
                values_to_plot,
                var_names_list,
                key=key,
                gene_symbols=gene_symbols,
            )
            title = values_to_plot
            if values_to_plot == 'logfoldchanges':
                title = 'log fold change'
            else:
                title = values_to_plot.replace("_", " ").replace('pvals', 'p-value')

        if plot_type == 'dotplot':
            from .._dotplot import dotplot

            _pl = dotplot(
                adata,
                var_names,
                groupby,
                dot_color_df=values_df,
                return_fig=True,
                gene_symbols=gene_symbols,
                **kwds,
            )
            if title is not None and 'colorbar_title' not in kwds:
                _pl.legend(colorbar_title=title)
        elif plot_type == 'matrixplot':
            from .._matrixplot import matrixplot

            _pl = matrixplot(
                adata,
                var_names,
                groupby,
                values_df=values_df,
                return_fig=True,
                gene_symbols=gene_symbols,
                **kwds,
            )

            if title is not None and 'colorbar_title' not in kwds:
                _pl.legend(title=title)

        return _fig_show_save_or_axes(_pl, return_fig, show, save)

    elif plot_type == 'stacked_violin':
        from .._stacked_violin import stacked_violin

        _pl = stacked_violin(
            adata,
            var_names,
            groupby,
            return_fig=True,
            gene_symbols=gene_symbols,
            **kwds,
        )
        return _fig_show_save_or_axes(_pl, return_fig, show, save)
    elif plot_type == 'heatmap':
        from .._anndata import heatmap

        return heatmap(
            adata,
            var_names,
            groupby,
            show=show,
            save=save,
            gene_symbols=gene_symbols,
            **kwds,
        )

    elif plot_type == 'tracksplot':
        from .._anndata import tracksplot

        return tracksplot(
            adata,
            var_names,
            groupby,
            show=show,
            save=save,
            gene_symbols=gene_symbols,
            **kwds,
        )


@_doc_params(params=doc_rank_genes_groups_plot_args, show_save_ax=doc_show_save_ax)
def rank_genes_groups_heatmap(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    gene_symbols: Optional[str] = None,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    min_logfoldchange: Optional[float] = None,
    key: str = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    **kwds,
):
    """\
    Plot ranking of genes using heatmap plot (see :func:`~scanpy.pl.heatmap`)

    Parameters
    ----------
    {params}
    {show_save_ax}
    **kwds
        Are passed to :func:`~scanpy.pl.heatmap`.
    {show_save_ax}

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.rank_genes_groups(adata, 'bulk_labels')
        sc.pl.rank_genes_groups_heatmap(adata)

    Show gene names per group on the heatmap

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_heatmap(adata, show_gene_labels=True)

    Plot top 5 genes per group (default 10 genes)

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True)

    .. currentmodule:: scanpy

    See also
    --------
    tl.rank_genes_groups
    tl.dendrogram
    """
    return _rank_genes_groups_plot(
        adata,
        plot_type='heatmap',
        groups=groups,
        n_genes=n_genes,
        gene_symbols=gene_symbols,
        groupby=groupby,
        var_names=var_names,
        key=key,
        min_logfoldchange=min_logfoldchange,
        show=show,
        save=save,
        **kwds,
    )


@_doc_params(params=doc_rank_genes_groups_plot_args, show_save_ax=doc_show_save_ax)
def rank_genes_groups_tracksplot(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    gene_symbols: Optional[str] = None,
    min_logfoldchange: Optional[float] = None,
    key: Optional[str] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    **kwds,
):
    """\
    Plot ranking of genes using heatmap plot (see :func:`~scanpy.pl.heatmap`)

    Parameters
    ----------
    {params}
    {show_save_ax}
    **kwds
        Are passed to :func:`~scanpy.pl.tracksplot`.
    {show_save_ax}

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.rank_genes_groups(adata, 'bulk_labels')
        sc.pl.rank_genes_groups_tracksplot(adata)
    """

    return _rank_genes_groups_plot(
        adata,
        plot_type='tracksplot',
        groups=groups,
        n_genes=n_genes,
        var_names=var_names,
        gene_symbols=gene_symbols,
        groupby=groupby,
        key=key,
        min_logfoldchange=min_logfoldchange,
        show=show,
        save=save,
        **kwds,
    )


@_doc_params(
    params=doc_rank_genes_groups_plot_args,
    vals_to_plot=doc_rank_genes_groups_values_to_plot,
    show_save_ax=doc_show_save_ax,
)
def rank_genes_groups_dotplot(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    values_to_plot: Optional[
        Literal[
            'scores',
            'logfoldchanges',
            'pvals',
            'pvals_adj',
            'log10_pvals',
            'log10_pvals_adj',
        ]
    ] = None,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    gene_symbols: Optional[str] = None,
    min_logfoldchange: Optional[float] = None,
    key: Optional[str] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    return_fig: Optional[bool] = False,
    **kwds,
):
    """\
    Plot ranking of genes using dotplot plot (see :func:`~scanpy.pl.dotplot`)

    Parameters
    ----------
    {params}
    {vals_to_plot}
    {show_save_ax}
    return_fig
        Returns :class:`DotPlot` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    **kwds
        Are passed to :func:`~scanpy.pl.dotplot`.

    Returns
    -------
    If `return_fig` is `True`, returns a :class:`DotPlot` object,
    else if `show` is false, return axes dict

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.rank_genes_groups(adata, 'bulk_labels', n_genes=adata.raw.shape[1])

    Plot top 2 genes per group.

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_dotplot(adata,n_genes=2)

    Plot with scaled expressions for easier identification of differences.

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_dotplot(adata, n_genes=2, standard_scale='var')

    Plot `logfoldchanges` instead of gene expression. In this case a diverging colormap
    like `bwr` or `seismic` works better. To center the colormap in zero, the minimum
    and maximum values to plot are set to -4 and 4 respectively.
    Also, only genes with a log fold change of 3 or more are shown.

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=4,
            values_to_plot="logfoldchanges", cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change'
        )

    Also, the last genes can be plotted. This can be useful to identify genes
    that are lowly expressed in a group. For this `n_genes=-4` is used

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=-4,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )

    A list specific genes can be given to check their log fold change. If a
    dictionary, the dictionary keys will be added as labels in the plot.

    .. plot::
        :context: close-figs

        var_names = {{'T-cell': ['CD3D', 'CD3E', 'IL32'],
                      'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
                      'myeloid': ['CST3', 'LYZ'] }}
        sc.pl.rank_genes_groups_dotplot(
            adata,
            var_names=var_names,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )

    .. currentmodule:: scanpy

    See also
    --------
    tl.rank_genes_groups
    """
    return _rank_genes_groups_plot(
        adata,
        plot_type='dotplot',
        groups=groups,
        n_genes=n_genes,
        groupby=groupby,
        values_to_plot=values_to_plot,
        var_names=var_names,
        gene_symbols=gene_symbols,
        key=key,
        min_logfoldchange=min_logfoldchange,
        show=show,
        save=save,
        return_fig=return_fig,
        **kwds,
    )


@_doc_params(params=doc_rank_genes_groups_plot_args, show_save_ax=doc_show_save_ax)
def rank_genes_groups_stacked_violin(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    gene_symbols: Optional[str] = None,
    *,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    min_logfoldchange: Optional[float] = None,
    key: Optional[str] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    return_fig: Optional[bool] = False,
    **kwds,
):
    """\
    Plot ranking of genes using stacked_violin plot
    (see :func:`~scanpy.pl.stacked_violin`)

    Parameters
    ----------
    {params}
    {show_save_ax}
    return_fig
        Returns :class:`StackedViolin` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    **kwds
        Are passed to :func:`~scanpy.pl.stacked_violin`.

    Returns
    -------
    If `return_fig` is `True`, returns a :class:`StackedViolin` object,
    else if `show` is false, return axes dict

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels')

    >>> sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=4,
    ... min_logfoldchange=4, figsize=(8,6))

    """

    return _rank_genes_groups_plot(
        adata,
        plot_type='stacked_violin',
        groups=groups,
        n_genes=n_genes,
        gene_symbols=gene_symbols,
        groupby=groupby,
        var_names=var_names,
        key=key,
        min_logfoldchange=min_logfoldchange,
        show=show,
        save=save,
        return_fig=return_fig,
        **kwds,
    )


@_doc_params(
    params=doc_rank_genes_groups_plot_args,
    vals_to_plot=doc_rank_genes_groups_values_to_plot,
    show_save_ax=doc_show_save_ax,
)
def rank_genes_groups_matrixplot(
    adata: AnnData,
    groups: Union[str, Sequence[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None,
    values_to_plot: Optional[
        Literal[
            'scores',
            'logfoldchanges',
            'pvals',
            'pvals_adj',
            'log10_pvals',
            'log10_pvals_adj',
        ]
    ] = None,
    var_names: Optional[Union[Sequence[str], Mapping[str, Sequence[str]]]] = None,
    gene_symbols: Optional[str] = None,
    min_logfoldchange: Optional[float] = None,
    key: Optional[str] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
    return_fig: Optional[bool] = False,
    **kwds,
):
    """\
    Plot ranking of genes using matrixplot plot (see :func:`~scanpy.pl.matrixplot`)

    Parameters
    ----------
    {params}
    {vals_to_plot}
    {show_save_ax}
    return_fig
        Returns :class:`MatrixPlot` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    **kwds
        Are passed to :func:`~scanpy.pl.matrixplot`.

    Returns
    -------
    If `return_fig` is `True`, returns a :class:`MatrixPlot` object,
    else if `show` is false, return axes dict

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.rank_genes_groups(adata, 'bulk_labels', n_genes=adata.raw.shape[1])

    Plot `logfoldchanges` instead of gene expression. In this case a diverging colormap
    like `bwr` or `seismic` works better. To center the colormap in zero, the minimum
    and maximum values to plot are set to -4 and 4 respectively.
    Also, only genes with a log fold change of 3 or more are shown.


    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_matrixplot(
            adata,
            n_genes=4,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )

    Also, the last genes can be plotted. This can be useful to identify genes
    that are lowly expressed in a group. For this `n_genes=-4` is used

    .. plot::
        :context: close-figs

        sc.pl.rank_genes_groups_matrixplot(
            adata,
            n_genes=-4,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )

    A list specific genes can be given to check their log fold change. If a
    dictionary, the dictionary keys will be added as labels in the plot.

    .. plot::
        :context: close-figs

        var_names = {{"T-cell": ['CD3D', 'CD3E', 'IL32'],
                      'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
                      'myeloid': ['CST3', 'LYZ'] }}
        sc.pl.rank_genes_groups_matrixplot(
            adata,
            var_names=var_names,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )
    """

    return _rank_genes_groups_plot(
        adata,
        plot_type='matrixplot',
        groups=groups,
        n_genes=n_genes,
        groupby=groupby,
        values_to_plot=values_to_plot,
        var_names=var_names,
        gene_symbols=gene_symbols,
        key=key,
        min_logfoldchange=min_logfoldchange,
        show=show,
        save=save,
        return_fig=return_fig,
        **kwds,
    )


@_doc_params(show_save_ax=doc_show_save_ax)
def rank_genes_groups_violin(
    adata: AnnData,
    groups: Optional[Sequence[str]] = None,
    n_genes: int = 20,
    gene_names: Optional[Iterable[str]] = None,
    gene_symbols: Optional[str] = None,
    use_raw: Optional[bool] = None,
    key: Optional[str] = None,
    split: bool = True,
    scale: str = 'width',
    strip: bool = True,
    jitter: Union[int, float, bool] = True,
    size: int = 1,
    ax: Optional[Axes] = None,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
):
    """\
    Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        List of group names.
    n_genes
        Number of genes to show. Is ignored if `gene_names` is passed.
    gene_names
        List of genes to plot. Is only useful if interested in a custom gene list,
        which is not the result of :func:`scanpy.tl.rank_genes_groups`.
    gene_symbols
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names` displayed in the plot.
    use_raw
        Use `raw` attribute of `adata` if present. Defaults to the value that
        was used in :func:`~scanpy.tl.rank_genes_groups`.
    split
        Whether to split the violins or not.
    scale
        See :func:`~seaborn.violinplot`.
    strip
        Show a strip plot on top of the violin plot.
    jitter
        If set to 0, no points are drawn. See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    {show_save_ax}
    """
    if key is None:
        key = 'rank_genes_groups'
    groups_key = str(adata.uns[key]['params']['groupby'])
    if use_raw is None:
        use_raw = bool(adata.uns[key]['params']['use_raw'])
    reference = str(adata.uns[key]['params']['reference'])
    groups_names = adata.uns[key]['names'].dtype.names if groups is None else groups
    if isinstance(groups_names, str):
        groups_names = [groups_names]
    axs = []
    for group_name in groups_names:
        if gene_names is None:
            _gene_names = adata.uns[key]['names'][group_name][:n_genes]
        else:
            _gene_names = gene_names
        if isinstance(_gene_names, np.ndarray):
            _gene_names = _gene_names.tolist()
        df = obs_df(adata, _gene_names, use_raw=use_raw, gene_symbols=gene_symbols)
        new_gene_names = df.columns
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

        _ax = sns.violinplot(
            x=x,
            y=y,
            data=df_tidy,
            inner=None,
            hue_order=hue_order,
            hue='hue',
            split=split,
            scale=scale,
            orient='vertical',
            ax=ax,
        )
        if strip:
            _ax = sns.stripplot(
                x=x,
                y=y,
                data=df_tidy,
                hue='hue',
                dodge=True,
                hue_order=hue_order,
                jitter=jitter,
                color='black',
                size=size,
                ax=_ax,
            )
        _ax.set_xlabel('genes')
        _ax.set_title('{} vs. {}'.format(group_name, reference))
        _ax.legend_.remove()
        _ax.set_ylabel('expression')
        _ax.set_xticklabels(new_gene_names, rotation='vertical')
        writekey = (
            f"rank_genes_groups_"
            f"{adata.uns[key]['params']['groupby']}_"
            f"{group_name}"
        )
        savefig_or_show(writekey, show=show, save=save)
        axs.append(_ax)
    if show is False:
        return axs


def sim(
    adata,
    tmax_realization: Optional[int] = None,
    as_heatmap: bool = False,
    shuffle: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
):
    """\
    Plot results of simulation.

    Parameters
    ----------
    tmax_realization
        Number of observations in one realization of the time series. The data matrix
        adata.X consists in concatenated realizations.
    as_heatmap
        Plot the timeseries as heatmap.
    shuffle
        Shuffle the data.
    show
        Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.
    """
    if tmax_realization is not None:
        tmax = tmax_realization
    elif 'tmax_write' in adata.uns:
        tmax = adata.uns['tmax_write']
    else:
        tmax = adata.n_obs
    n_realizations = adata.n_obs / tmax
    if not shuffle:
        if not as_heatmap:
            timeseries(
                adata.X,
                var_names=adata.var_names,
                xlim=[0, 1.25 * adata.n_obs],
                highlights_x=np.arange(tmax, n_realizations * tmax, tmax),
                xlabel='realizations',
            )
        else:
            # plot time series as heatmap, as in Haghverdi et al. (2016), Fig. 1d
            timeseries_as_heatmap(
                adata.X,
                var_names=adata.var_names,
                highlights_x=np.arange(tmax, n_realizations * tmax, tmax),
            )
        pl.xticks(
            np.arange(0, n_realizations * tmax, tmax),
            np.arange(n_realizations).astype(int) + 1,
        )
        savefig_or_show('sim', save=save, show=show)
    else:
        # shuffled data
        X = adata.X
        X, rows = subsample(X, seed=1)
        timeseries(
            X,
            var_names=adata.var_names,
            xlim=[0, 1.25 * adata.n_obs],
            highlights_x=np.arange(tmax, n_realizations * tmax, tmax),
            xlabel='index (arbitrary order)',
        )
        savefig_or_show('sim_shuffled', save=save, show=show)


@_doc_params(
    vminmax=doc_vbound_percentile, panels=doc_panels, show_save_ax=doc_show_save_ax
)
def embedding_density(
    adata: AnnData,
    # on purpose, there is no asterisk here (for backward compat)
    basis: str = 'umap',  # was positional before 1.4.5
    key: Optional[str] = None,  # was positional before 1.4.5
    groupby: Optional[str] = None,
    group: Optional[Union[str, List[str], None]] = 'all',
    color_map: Union[Colormap, str] = 'YlOrRd',
    bg_dotsize: Optional[int] = 80,
    fg_dotsize: Optional[int] = 180,
    vmax: Optional[int] = 1,
    vmin: Optional[int] = 0,
    vcenter: Optional[int] = None,
    norm: Optional[Normalize] = None,
    ncols: Optional[int] = 4,
    hspace: Optional[float] = 0.25,
    wspace: Optional[None] = None,
    title: str = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Axes] = None,
    return_fig: Optional[bool] = None,
    **kwargs,
) -> Union[Figure, Axes, None]:
    """\
    Plot the density of cells in an embedding (per condition).

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
        Name of the `.obs` covariate that contains the density estimates. Alternatively, pass `groupby`.
    groupby
        Name of the condition used in `tl.embedding_density`. Alternatively, pass `key`.
    group
        The category in the categorical observation annotation to be plotted.
        For example, 'G1' in the cell cycle 'phase' covariate. If all categories
        are to be plotted use group='all' (default), If multiple categories
        want to be plotted use a list (e.g.: ['G1', 'S']. If the overall density
        wants to be ploted set group to 'None'.
    color_map
        Matplolib color map to use for density plotting.
    bg_dotsize
        Dot size for background data points not in the `group`.
    fg_dotsize
        Dot size for foreground data points in the `group`.
    {vminmax}
    {panels}
    {show_save_ax}

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.umap(adata)
        sc.tl.embedding_density(adata, basis='umap', groupby='phase')

    Plot all categories be default

    .. plot::
        :context: close-figs

        sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase')

    Plot selected categories

    .. plot::
        :context: close-figs

        sc.pl.embedding_density(
            adata,
            basis='umap',
            key='umap_density_phase',
            group=['G1', 'S'],
        )

    .. currentmodule:: scanpy

    See also
    --------
    tl.embedding_density
    """
    sanitize_anndata(adata)

    # Test user inputs
    basis = basis.lower()

    if basis == 'fa':
        basis = 'draw_graph_fa'

    if key is not None and groupby is not None:
        raise ValueError('either pass key or groupby but not both')

    if key is None:
        key = 'umap_density'
    if groupby is not None:
        key += f'_{groupby}'

    if f'X_{basis}' not in adata.obsm_keys():
        raise ValueError(
            f'Cannot find the embedded representation `adata.obsm[X_{basis!r}]`. '
            'Compute the embedding first.'
        )

    if key not in adata.obs or f'{key}_params' not in adata.uns:
        raise ValueError(
            'Please run `sc.tl.embedding_density()` first '
            'and specify the correct key.'
        )

    if 'components' in kwargs:
        logg.warning(
            'Components were specified, but will be ignored. Only the '
            'components used to calculate the density can be plotted.'
        )
        del kwargs['components']

    components = adata.uns[f'{key}_params']['components']
    groupby = adata.uns[f'{key}_params']['covariate']

    # turn group into a list if needed
    if group == 'all':
        if groupby is None:
            group = None
        else:
            group = list(adata.obs[groupby].cat.categories)
    elif isinstance(group, str):
        group = [group]

    if group is None and groupby is not None:
        raise ValueError(
            'Densities were calculated over an `.obs` covariate. '
            'Please specify a group from this covariate to plot.'
        )

    if group is not None and groupby is None:
        logg.warning(
            "value of 'group' is ignored because densities "
            "were not calculated for an `.obs` covariate."
        )
        group = None

    if np.min(adata.obs[key]) < 0 or np.max(adata.obs[key]) > 1:
        raise ValueError('Densities should be scaled between 0 and 1.')

    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02

    # Make the color map
    if isinstance(color_map, str):
        color_map = copy(cm.get_cmap(color_map))

    color_map.set_over('black')
    color_map.set_under('lightgray')
    # a name to store the density values is needed. To avoid
    # overwriting a user name a new random name is created
    while True:
        col_id = np.random.randint(1000, 10000)
        density_col_name = f'_tmp_embedding_density_column_{col_id}_'
        if density_col_name not in adata.obs.columns:
            break

    # if group is set, then plot it using multiple panels
    # (even if only one group is set)
    if (
        group is not None
        and not isinstance(group, str)
        and isinstance(group, cabc.Sequence)
    ):
        if ax is not None:
            raise ValueError("Can only specify `ax` if no `group` sequence is given.")
        fig, gs = _panel_grid(hspace, wspace, ncols, len(group))

        axs = []
        for count, group_name in enumerate(group):
            if group_name not in adata.obs[groupby].cat.categories:
                raise ValueError(
                    'Please specify a group from the `.obs` category '
                    'over which the density was calculated. '
                    f'Invalid group name: {group_name}'
                )

            ax = pl.subplot(gs[count])
            # Define plotting data
            dot_sizes = np.ones(adata.n_obs) * bg_dotsize
            group_mask = adata.obs[groupby] == group_name
            dens_values = -np.ones(adata.n_obs)
            dens_values[group_mask] = adata.obs[key][group_mask]
            adata.obs[density_col_name] = dens_values
            dot_sizes[group_mask] = np.ones(sum(group_mask)) * fg_dotsize

            if title is None:
                _title = group_name
            else:
                _title = title

            ax = embedding(
                adata,
                basis,
                dimensions=np.array(components) - 1,  # Saved with 1 based indexing
                color=density_col_name,
                color_map=color_map,
                size=dot_sizes,
                vmax=vmax,
                vmin=vmin,
                vcenter=vcenter,
                norm=norm,
                save=False,
                title=_title,
                ax=ax,
                show=False,
                **kwargs,
            )
            axs.append(ax)

        ax = axs
    else:
        dens_values = adata.obs[key]
        dot_sizes = np.ones(adata.n_obs) * fg_dotsize

        adata.obs[density_col_name] = dens_values

        # Ensure title is blank as default
        if title is None:
            title = group if group is not None else ""

        # Plot the graph
        fig_or_ax = embedding(
            adata,
            basis,
            dimensions=np.array(components) - 1,  # Saved with 1 based indexing
            color=density_col_name,
            color_map=color_map,
            size=dot_sizes,
            vmax=vmax,
            vmin=vmin,
            vcenter=vcenter,
            norm=norm,
            save=False,
            show=False,
            title=title,
            ax=ax,
            return_fig=return_fig,
            **kwargs,
        )
        if return_fig:
            fig = fig_or_ax
        else:
            ax = fig_or_ax

    # remove temporary column name
    adata.obs = adata.obs.drop(columns=[density_col_name])

    if return_fig:
        return fig
    savefig_or_show(f"{key}_", show=show, save=save)
    if show is False:
        return ax


def _get_values_to_plot(
    adata,
    values_to_plot: Literal[
        'scores',
        'logfoldchanges',
        'pvals',
        'pvals_adj',
        'log10_pvals',
        'log10_pvals_adj',
    ],
    gene_names: Sequence[str],
    groups: Optional[Sequence[str]] = None,
    key: Optional[str] = 'rank_genes_groups',
    gene_symbols: Optional[str] = None,
):
    """
    If rank_genes_groups has been called, this function
    prepares a dataframe containing scores, pvalues, logfoldchange etc to be plotted
    as dotplot or matrixplot.

    The dataframe index are the given groups and the columns are the gene_names

    used by rank_genes_groups_dotplot

    Parameters
    ----------
    adata
    values_to_plot
        name of the value to plot
    gene_names
        gene names
    groups
        groupby categories
    key
        adata.uns key where the rank_genes_groups is stored.
        By default 'rank_genes_groups'
    gene_symbols
        Key for field in .var that stores gene symbols.
    Returns
    -------
    pandas DataFrame index=groups, columns=gene_names

    """
    valid_options = [
        'scores',
        'logfoldchanges',
        'pvals',
        'pvals_adj',
        'log10_pvals',
        'log10_pvals_adj',
    ]
    if values_to_plot not in valid_options:
        raise ValueError(
            f"given value_to_plot: '{values_to_plot}' is not valid. Valid options are {valid_options}"
        )

    values_df = None
    check_done = False
    if groups is None:
        groups = adata.uns[key]['names'].dtype.names
    if values_to_plot is not None:

        df_list = []
        for group in groups:
            df = rank_genes_groups_df(adata, group, key=key, gene_symbols=gene_symbols)
            if gene_symbols is not None:
                df['names'] = df[gene_symbols]
            # check that all genes are present in the df as sc.tl.rank_genes_groups
            # can be called with only top genes
            if not check_done:
                if df.shape[0] < adata.shape[1]:
                    message = (
                        "Please run `sc.tl.rank_genes_groups` with "
                        "'n_genes=adata.shape[1]' to save all gene "
                        f"scores. Currently, only {df.shape[0]} "
                        "are found"
                    )
                    logg.error(message)
                    raise ValueError(message)
            df['group'] = group
            df_list.append(df)

        values_df = pd.concat(df_list)
        if values_to_plot.startswith('log10'):
            column = values_to_plot.replace('log10_', '')
        else:
            column = values_to_plot
        values_df = pd.pivot(
            values_df, index='names', columns='group', values=column
        ).fillna(1)

        if values_to_plot in ['log10_pvals', 'log10_pvals_adj']:
            values_df = -1 * np.log10(values_df)

        values_df = values_df.loc[gene_names].T

    return values_df
