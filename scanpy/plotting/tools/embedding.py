from matplotlib import pyplot as pl
from .. import utils
from ..anndata import scatter
from ...utils import doc_params
from ... import settings
from ..docs import doc_adata_color_etc, doc_edges_arrows, doc_scatter_bulk, doc_show_save_ax


@doc_params(adata_color_etc=doc_adata_color_etc, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def diffmap(
        adata,
        color=None,
        use_raw=None,
        sort_order=True,
        alpha=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        frameon=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """\
    Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    if components == 'all':
        components_list = ['{},{}'.format(*((i, i+1) if i % 2 == 1 else (i+1, i)))
            for i in range(1, adata.obsm['X_diffmap'].shape[1])]
    else:
        if components is None: components = '1,2' if '2d' in projection else '1,2,3'
        if not isinstance(components, list): components_list = [components]
        else: components_list = components
    for components in components_list:
        axs = scatter(
            adata,
            basis='diffmap',
            color=color,
            use_raw=use_raw,
            sort_order=sort_order,
            alpha=alpha,
            groups=groups,
            components=components,
            projection=projection,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            legend_fontweight=legend_fontweight,
            color_map=color_map,
            palette=palette,
            frameon=frameon,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False,
            save=False,
            ax=ax)
        writekey = 'diffmap'
        if isinstance(components, list): components = ','.join(
            [str(comp) for comp in components])
        writekey += ('_components' + components.replace(',', '')
                     + (save if isinstance(save, str) else ''))
        if settings.autosave or (save is not None):
            utils.savefig(writekey)
    show = settings.autoshow if show is None else show
    if show: pl.show()
    if show == False: return axs


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def draw_graph(
        adata,
        layout=None,
        color=None,
        use_raw=None,
        edges=False,
        edges_width=0.1,
        edges_color='grey',
        arrows=False,
        arrows_kwds=None,
        sort_order=True,
        alpha=None,
        groups=None,
        components=None,
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        frameon=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """\
    Scatter plot in graph-drawing basis.

    Parameters
    ----------
    {adata_color_etc}
    layout : {{'fa', 'fr', 'drl', ...}}, optional (default: last computed)
        One of the `draw_graph` layouts, see
        :func:`~scanpy.api.tl.draw_graph`. By default, the last computed layout
        is used.
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    if layout is None: layout = str(adata.uns['draw_graph']['params']['layout'])
    basis = 'draw_graph_' + layout
    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError('Did not find {} in adata.obs. Did you compute layout {}?'
                         .format('draw_graph_' + layout, layout))
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        use_raw=use_raw,
        sort_order=sort_order,
        alpha=alpha,
        groups=groups,
        components=components,
        projection='2d',
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False,
        ax=ax)
    if edges: utils.plot_edges(axs, adata, basis, edges_width, edges_color)
    if arrows: utils.plot_arrows(axs, adata, basis, arrows_kwds)
    utils.savefig_or_show(basis, show=show, save=save)
    if show == False: return axs


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def tsne(
        adata,
        color=None,
        use_raw=None,
        edges=False,
        edges_width=0.1,
        edges_color='grey',
        arrows=False,
        arrows_kwds=None,
        sort_order=True,
        alpha=None,
        groups=None,
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        frameon=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """\
    Scatter plot in tSNE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    basis = 'tsne'
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        use_raw=use_raw,
        sort_order=sort_order,
        alpha=alpha,
        groups=groups,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False,
        ax=ax)
    if edges: utils.plot_edges(axs, adata, basis, edges_width, edges_color)
    if arrows: utils.plot_arrows(axs, adata, basis, arrows_kwds)
    utils.savefig_or_show(basis, show=show, save=save)
    if show == False: return axs


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def umap(
        adata,
        color=None,
        use_raw=None,
        edges=False,
        edges_width=0.1,
        edges_color='grey',
        arrows=False,
        arrows_kwds=None,
        sort_order=True,
        alpha=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        frameon=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """\
    Scatter plot in UMAP basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    basis = 'umap'
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        use_raw=use_raw,
        sort_order=sort_order,
        alpha=alpha,
        groups=groups,
        components=components,
        projection=projection,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False,
        ax=ax)
    if edges: utils.plot_edges(axs, adata, basis, edges_width, edges_color)
    if arrows: utils.plot_arrows(axs, adata, basis, arrows_kwds)
    utils.savefig_or_show(basis, show=show, save=save)
    if show == False: return axs


@doc_params(adata_color_etc=doc_adata_color_etc, edges_arrows=doc_edges_arrows, scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def phate(
        adata,
        color=None,
        use_raw=None,
        edges=False,
        edges_width=0.1,
        edges_color='grey',
        arrows=False,
        arrows_kwds=None,
        sort_order=True,
        alpha=None,
        groups=None,
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        frameon=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """\
    Scatter plot in PHATE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False`, a list of `matplotlib.Axis` objects. Every second element
    corresponds to the 'right margin' drawing area for color bars and legends.

    Examples
    --------
    >>> import scanpy.api as sc
    >>> import phate
    >>> data, branches = phate.tree.gen_dla(n_dim=100,
                                            n_branch=20,
                                            branch_length=100)
    >>> data.shape
    (2000, 100)
    >>> adata = sc.AnnData(data)
    >>> adata.obs['branches'] = branches
    >>> sc.tl.phate(adata, k=5, a=20, t=150)
    >>> adata.obsm['X_phate'].shape
    (2000, 2)
    >>> sc.pl.phate(adata,
                    color='branches',
                    color_map='tab20')
    """
    basis = 'phate'
    axs = scatter(
        adata,
        basis=basis,
        color=color,
        use_raw=use_raw,
        sort_order=sort_order,
        alpha=alpha,
        groups=groups,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        frameon=frameon,
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False,
        ax=ax)
    if edges:
        utils.plot_edges(axs, adata, basis, edges_width, edges_color)
    if arrows:
        utils.plot_arrows(axs, adata, basis, arrows_kwds)
    utils.savefig_or_show(basis, show=show, save=save)
    if show == False:
        return axs
