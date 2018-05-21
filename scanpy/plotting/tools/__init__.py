"""Plotting

Plotting functions for each tool and toplevel plotting functions for AnnData.
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams

from .. import utils
from ...utils import doc_params
from ... import settings
from ... import logging as logg
from ..anndata import scatter, ranking
from ..utils import timeseries, timeseries_subplot, timeseries_as_heatmap
from ..utils import doc_edges_arrows


# ------------------------------------------------------------------------------
# Visualization tools
# ------------------------------------------------------------------------------


def pca(adata, **params):
    """Plot PCA results.

    The parameters are the ones of the scatter plot. Call pca_ranking separately
    if you want to change the default settings.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    """
    show = params['show'] if 'show' in params else None
    if 'show' in params: del params['show']
    pca_scatter(adata, **params, show=False)
    pca_loadings(adata, show=False)
    pca_variance_ratio(adata, show=show)


def pca_scatter(
        adata,
        color=None,
        use_raw=True,
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """Scatter plot in PCA coordinates.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    axs = scatter(
        adata,
        basis='pca',
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
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False, ax=ax)
    utils.savefig_or_show('pca_scatter', show=show, save=save)
    if show == False: return axs


def pca_loadings(adata, components=None, show=None, save=None):
    """Rank genes according to contributions to PCs.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    components : str or list of integers, optional
        For example, ``'1,2,3'`` means ``[1, 2, 3]``, first, second, third
        principal component.
    show : bool, optional (default: None)
        Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    """
    if components is None: components = [1, 2, 3]
    elif isinstance(components, str): components = components.split(',')
    components = np.array(components) - 1
    ranking(adata, 'varm', 'PCs', indices=components)
    utils.savefig_or_show('pca_loadings', show=show, save=save)


def pca_variance_ratio(adata, log=False, show=None, save=None):
    """Plot the variance ratio.

    Parameters
    ----------
    show : bool, optional (default: None)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    """
    ranking(adata, 'uns', 'variance_ratio', dictionary='pca', labels='PC', log=log)
    utils.savefig_or_show('pca_variance_ratio', show=show, save=save)


def diffmap(
        adata,
        color=None,
        use_raw=True,
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object. Only works if plotting a single component.

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
    if not settings.autosave and show: pl.show()
    if show == False: return axs


@doc_params(edges_arrows=doc_edges_arrows)
def draw_graph(
        adata,
        layout=None,
        color=None,
        use_raw=True,
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
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    layout : {{'fr', 'drl', ...}}, optional (default: last computed)
        One of the `draw_graph` layouts, see
        :func:`~scanpy.api.tl.draw_graph`. By default, the last computed layout
        is used.
    color : `str` or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {edges_arrows}
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.

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


@doc_params(edges_arrows=doc_edges_arrows)
def tsne(
        adata,
        color=None,
        use_raw=True,
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in tSNE basis.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {edges_arrows}
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
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
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.

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


@doc_params(edges_arrows=doc_edges_arrows)
def umap(
        adata,
        color=None,
        use_raw=True,
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in UMAP basis.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {{edges_arrows}}
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {{'2d', '3d'}}, optional (default: '2d')
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.

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


@doc_params(edges_arrows=doc_edges_arrows)
def phate(
        adata,
        color=None,
        use_raw=True,
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Scatter plot in PHATE basis.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {edges_arrows}
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
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
    right_margin : float or list of floats (default: None)
         Adjust the width of the space right of each plotting panel.
    size : float (default: None)
         Point size.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.

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


# ------------------------------------------------------------------------------
# Subgroup identification and ordering - clustering, pseudotime, branching
# and tree inference tools
# ------------------------------------------------------------------------------


def dpt(
        adata,
        basis='diffmap',
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None):
    """Plot results of DPT analysis.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    basis : {`'diffmap'`, `'pca'`, `'tsne'`, `'draw_graph_...'`}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: None)
        Observation/ cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.
    show_tree : bool, optional (default: False)
         This shows the inferred tree. For more than a single branching, the
         result is pretty unreliable. Use tool `aga` (Approximate Graph
         Abstraction) instead.
    """
    colors = ['dpt_pseudotime']
    if len(np.unique(adata.obs['dpt_groups'].values)) > 1: colors += ['dpt_groups']
    if color is not None: colors = color
    dpt_scatter(
        adata,
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
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=save)
    dpt_groups_pseudotime(adata, color_map=color_map, show=False, save=save)
    dpt_timeseries(adata, color_map=color_map, show=show, save=save)


def dpt_scatter(
        adata,
        basis='diffmap',
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None):
    """Scatter plot of DPT results.

    See parameters of sc.pl.dpt().
    """

    colors = ['dpt_pseudotime']
    if len(np.unique(adata.obs['dpt_groups'].values)) > 1: colors += ['dpt_groups']
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
            legend_fontweight=legend_fontweight,
            color_map=color_map,
            palette=palette,
            right_margin=right_margin,
            size=size,
            title=title,
            show=False)
        writekey = 'dpt_' + basis + '_components' + components.replace(',', '')
        utils.savefig_or_show(writekey, show=show, save=save)


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


def louvain(
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
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None, ax=None):
    """Plot results of Louvain clustering.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    basis : {`'diffmap'`, `'pca'`, `'tsne'`, `'draw_graph_...'`}
        Choose the basis in which to plot.
    color : string or list of strings, optional (default: None)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    groups : str, optional (default: all groups)
        Restrict to a few categories in categorical observation annotation.
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
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : matplotlib.Axes
         A matplotlib axes object.
    """
    add_color = []
    if color is not None:
        add_color = color if isinstance(color, list) else color.split(',')
    color = ['louvain_groups'] + add_color
    axs = scatter(
        adata,
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
        right_margin=right_margin,
        size=size,
        title=title,
        show=False,
        save=False)
    utils.savefig_or_show('louvain_' + basis, show=show, save=save)


def rank_genes_groups(adata, groups=None, n_genes=20, gene_symbols=None, key=None, fontsize=8, show=None, save=None, ext=None):
    """Plot ranking of genes.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
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
    show : `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : `matplotlib.Axes`, optional (default: `None`)
        A `matplotlib.Axes` object.
    """
    if key is None:
        key = 'rank_genes_groups'
    groups_key = str(adata.uns[key]['params']['groupby'])
    reference = str(adata.uns[key]['params']['reference'])
    group_names = (adata.uns[key]['names'].dtype.names
                   if groups is None else groups)
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
        gene_names = adata.uns[key]['names'][group_name]
        scores = adata.uns[key]['scores'][group_name]
        for ig, g in enumerate(gene_names[:n_genes]):
            gene_name = gene_names[ig]
            pl.text(ig, scores[ig], gene_name if gene_symbols is None else adata.var[gene_symbols][gene_name],
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center', fontsize=fontsize)
        pl.title('{} vs. {}'.format(group_name, reference))
        if n_panels <= 5 or count >= n_panels_x: pl.xlabel('ranking')
        if count == 0 or count == n_panels_x: pl.ylabel('score')
        ymin = np.min(scores)
        ymax = np.max(scores)
        ymax += 0.3*(ymax-ymin)
        pl.ylim([ymin, ymax])
        pl.xlim(-0.9, ig+1-0.1)
    writekey = ('rank_genes_groups_'
                + str(adata.uns[key]['params']['groupby']))
    utils.savefig_or_show(writekey, show=show, save=save)


def rank_genes_groups_violin(adata, groups=None, n_genes=20,
                             gene_names=None, gene_symbols=None,
                             use_raw=None,
                             key=None,
                             split=True,
                             scale='width',
                             strip=True, jitter=True, size=1,
                             ax=None, show=None, save=None):
    """Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    groups : list of `str`, optional (default: `None`)
        List of group names.
    n_genes : `int`, optional (default: 20)
        Number of genes to show. Not used if gene_names is not None.
    gene_names : list of `str`
        List of genes to plot.
    gene_symbols : `str`
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names`.
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present. Defaults to the value that
        was used in :func:`~scanpy.api.tl.rank_genes_groups`.
    split : `bool`, optional (default: `True`)
        Whether to split the violins or not.
    scale : `str` (default: 'width')
        See `seaborn.violinplot`.
    strip : `bool` (default: `True`)
        Show a strip plot on top of the violin plot.
    jitter : `int`, `float`, `bool`, optional (default: `True`)
        If set to 0, no points are drawn. See `seaborn.stripplot`.
    size : `int`, optional (default: 1)
        Size of the jitter points.
    show : `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    ax : `matplotlib.Axes`, optional (default: `None`)
        A `matplotlib.Axes` object.
    """
    if key is None:
        key = 'rank_genes_groups'
    from ..tools import rank_genes_groups
    groups_key = str(adata.uns[key]['params']['groupby'])
    if use_raw is None:
        use_raw = bool(adata.uns[key]['params']['use_raw'])
    reference = str(adata.uns[key]['params']['reference'])
    groups_names = (adata.uns[key]['names'].dtype.names
                    if groups is None else groups)
    if isinstance(groups_names, str): groups_names = [groups_names]
    for group_name in groups_names:
        if gene_names is None:
            gene_names = adata.uns[
                key]['names'][group_name][:n_genes]
        # make a "hue" option!
        df = pd.DataFrame()
        new_gene_names = []
        for g in gene_names:
            if adata.raw is not None and use_raw:
                X_col = adata.raw[:, g].X
            else:
                X_col = adata[:, g].X
            if issparse(X_col): X_col = X_col.toarray().flatten()
            new_gene_names.append(
                g if gene_symbols is None else adata.var[gene_symbols][g])
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
        ax = sns.violinplot(x=x, y=y, data=df_tidy, inner=None,
                            hue_order=hue_order, hue='hue', split=split,
                            scale=scale, orient='vertical', ax=ax)
        if strip:
            ax = sns.stripplot(x=x, y=y, data=df_tidy,
                               hue='hue', dodge=True, hue_order=hue_order,
                               jitter=jitter, color='black', size=size, ax=ax)
        ax.set_xlabel('genes')
        ax.set_title('{} vs. {}'.format(group_name, reference))
        ax.legend_.remove()
        ax.set_ylabel('expression')
        ax.set_xticklabels(gene_names, rotation='vertical')
        writekey = ('rank_genes_groups_'
                    + str(adata.uns[key]['params']['groupby'])
                    + '_' + group_name)
        utils.savefig_or_show(writekey, show=show, save=save)


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
    show : bool, optional (default: None)
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
