"""Plotting functions for AnnData.
"""

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import is_color_like
import seaborn as sns

from .. import settings
from . import utils
from .utils import scatter_base, scatter_group
from ..utils import sanitize_anndata

def scatter(
        adata,
        x=None,
        y=None,
        color='grey',
        use_raw=True,
        sort_order=True,
        alpha=None,
        basis=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        legend_fontweight=None,
        color_map=None,
        palette=None,
        right_margin=None,
        left_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """Scatter plot.

    Color with sample annotation (`color in adata.obs_keys()`) or gene
    expression (`color in adata.var_names`).

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    x : str or None
        x coordinate.
    y : str or None
        y coordinate.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    sort_order : `bool`, optional (default: `True`)
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    basis : {'pca', 'tsne', 'diffmap', 'draw_graph_fr', etc.}
        String that denotes a plotting tool that computed coordinates.
    groups : str, optional (default: all groups in color)
        Allows to restrict categories in sample annotation to a subset.
    components : str or list of str, optional (default: '1,2')
         String of the form '1,2' or ['1,2', '2,3'].
    projection : {'2d', '3d'}, optional (default: '2d')
         Projection of plot.
    legend_loc : str, optional (default: 'right margin')
         Location of legend, either 'on data', 'right margin' or valid keywords
         for matplotlib.legend.
    legend_fontsize : int (default: None)
         Legend font size.
    legend_fontweight : int (default: None)
         Legend font weight.
    color_map : str (default: 'viridis')
         String denoting matplotlib color map for continuous coloring.
    palette : list of str (default: None)
         Colors to use for plotting groups (categorical annotation).
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: None)
         Point size. Sample-number dependent by default.
    title : str, optional (default: None)
         Provide title for panels either as `["title1", "title2", ...]` or
         `"title1,title2,..."`.
    show : bool, optional (default: None)
         Show the plot.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    A list of matplotlib.Axis objects.
    """
    sanitize_anndata(adata)
    if components is None: components = '1,2' if '2d' in projection else '1,2,3'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1
    title = None if title is None else title.split(',') if isinstance(title, str) else title
    keys = ['grey'] if color is None else color.split(',') if isinstance(color, str) else color
    groups = None if groups is None else groups.split(',') if isinstance(groups, str) else groups
    highlights = adata.uns['highlights'] if 'highlights' in adata.uns else []
    if basis is not None:
        try:
            Y = adata.obsm['X_' + basis][:, components]
        except KeyError:
            raise KeyError('compute coordinates using visualization tool {} first'
                           .format(basis))
    elif x is not None and y is not None:
        x_arr = adata._get_obs_array(x)
        y_arr = adata._get_obs_array(y)
        Y = np.c_[x_arr[:, None], y_arr[:, None]]
    else:
        raise ValueError('Either provide keys for a `basis` or for `x` and `y`.')

    if size is None:
        n = Y.shape[0]
        size = 120000 / n

    if legend_loc == 'on data' and legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']
    elif legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']

    palette_was_none = False
    if palette is None: palette_was_none = True
    if isinstance(palette, list):
        if not is_color_like(palette[0]):
            palettes = palette
        else:
            palettes = [palette]
    else:
        palettes = [palette for i in range(len(keys))]
    for i, palette in enumerate(palettes):
        palettes[i] = utils.default_palette(palette)

    if basis is not None:
        component_name = ('DC' if basis == 'diffmap'
                          else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis
                          else 'tSNE' if basis == 'tsne'
                          else 'PC' if basis == 'pca'
                          else 'Spring' if basis == 'spring'
                          else None)
    else:
        component_name = None
    axis_labels = (x, y) if component_name is None else None
    show_ticks = True if component_name is None else False

    # the actual color ids, e.g. 'grey' or '#109482'
    color_ids = [None if not is_color_like(key)
                 else key for key in keys]
    categoricals = []
    colorbars = []
    for ikey, key in enumerate(keys):
        if color_ids[ikey] is not None:
            c = color_ids[ikey]
            continuous = True
            categorical = False
            colorbars.append(False)
        else:
            c = 'white' if projection == '2d' else 'white'
            categorical = False
            continuous = False
            # test whether we have categorial or continuous annotation
            if key in adata.obs_keys():
                if is_categorical_dtype(adata.obs[key]):
                    categorical = True
                else:
                    continuous = True
                    c = adata.obs[key]
            # coloring according to gene expression
            elif (use_raw
                  and adata.raw is not None
                  and key in adata.raw.var_names):
                c = adata.raw[:, key].X
                continuous = True
            elif key in adata.var_names:
                c = adata[:, key].X
                continuous = True
            else:
                raise ValueError('"' + key + '" is invalid!'
                                 + ' specify valid sample annotation, one of '
                                 + str(adata.obs_keys()) + ' or a gene name '
                                 + str(adata.var_names))
            colorbars.append(True if continuous else False)
        if categorical: categoricals.append(ikey)
        color_ids[ikey] = c

    if right_margin is None and len(categoricals) > 0:
        if legend_loc == 'right margin': right_margin = 0.5
    if title is None and keys[0] is not None:
        title = [key.replace('_', ' ') if not is_color_like(key) else '' for key in keys]

    axs = scatter_base(Y,
                       title=title,
                       alpha=alpha,
                       component_name=component_name,
                       axis_labels=axis_labels,
                       component_indexnames=components + 1,
                       projection=projection,
                       colors=color_ids,
                       highlights=highlights,
                       colorbars=colorbars,
                       right_margin=right_margin,
                       left_margin=left_margin,
                       sizes=[size for c in keys],
                       color_map=color_map,
                       show_ticks=show_ticks,
                       ax=ax)

    def add_centroid(centroids, name, Y, mask):
        Y_mask = Y[mask]
        if Y_mask.shape[0] == 0: return
        median = np.median(Y_mask, axis=0)
        i = np.argmin(np.sum(np.abs(Y_mask - median), axis=1))
        centroids[name] = Y_mask[i]

    for i, ikey in enumerate(categoricals):
        palette = palettes[i]
        key = keys[ikey]
        if (not key + '_colors' in adata.uns
            or not palette_was_none
            or len(adata.obs[key].cat.categories) != len(adata.uns[key + '_colors'])):
            utils.add_colors_for_categorical_sample_annotation(adata, key, palette)
        # actually plot the groups
        mask_remaining = np.ones(Y.shape[0], dtype=bool)
        centroids = {}
        if groups is None:
            for iname, name in enumerate(adata.obs[key].cat.categories):
                if name not in settings.categories_to_ignore:
                    mask = scatter_group(axs[ikey], key, iname,
                                         adata, Y, projection, size=size, alpha=alpha)
                    mask_remaining[mask] = False
                    if legend_loc == 'on data': add_centroid(centroids, name, Y, mask)
        else:
            for name in groups:
                if name not in set(adata.obs[key].cat.categories):
                    raise ValueError('"' + name + '" is invalid!'
                                     + ' specify valid name, one of '
                                     + str(adata.obs[key].cat.categories))
                else:
                    iname = np.flatnonzero(adata.obs[key].cat.categories.values == name)[0]
                    mask = scatter_group(axs[ikey], key, iname,
                                         adata, Y, projection, size=size, alpha=alpha)
                    if legend_loc == 'on data': add_centroid(centroids, name, Y, mask)
                    mask_remaining[mask] = False
        if mask_remaining.sum() > 0:
            data = [Y[mask_remaining, 0], Y[mask_remaining, 1]]
            if projection == '3d': data.append(Y[mask_remaining, 2])
            axs[ikey].scatter(*data, marker='.', c='grey', s=size,
                                    edgecolors='none', zorder=-1)
        legend = None
        if legend_loc == 'on data':
            for name, pos in centroids.items():
                axs[ikey].text(pos[0], pos[1], name,
                                     weight=legend_fontweight,
                                     verticalalignment='center',
                                     horizontalalignment='center',
                                     fontsize=legend_fontsize)
        elif legend_loc == 'right margin':
            legend = axs[ikey].legend(frameon=False, loc='center left',
                                            bbox_to_anchor=(1, 0.5),
                                            ncol=(1 if len(adata.obs[key].cat.categories) <= 14
                                                  else 2 if len(adata.obs[key].cat.categories) <= 30 else 3),
                                            fontsize=legend_fontsize)
        elif legend_loc != 'none':
            legend = axs[ikey].legend(frameon=False, loc=legend_loc,
                                            fontsize=legend_fontsize)
        if legend is not None:
            for handle in legend.legendHandles: handle.set_sizes([300.0])
    utils.savefig_or_show('scatter' if basis is None else basis, show=show, save=save)
    return axs


def ranking(adata, attr, keys, indices=None,
            labels=None, color='black', n_points=30,
            log=False):
    """Plot rankings.

    See, for example, how this is used in pl.pca_ranking.

    Parameters
    ----------
    adata : AnnData
        The data.
    attr : {'var', 'obs', 'add', 'varm', 'obsm'}
        The attribute of AnnData that contains the score.
    keys : str or list of str
        The scores to look up an array from the attribute of adata.

    Returns
    -------
    Returns matplotlib gridspec with access to the axes.
    """
    if isinstance(keys, str) and indices is not None:
        scores = getattr(adata, attr)[keys][:, indices]
        keys = ['{}{}'.format(keys[:-1], i+1) for i in indices]
    else:
        scores = getattr(adata, attr)[keys]
    n_panels = len(keys) if isinstance(keys, list) else 1
    if n_panels == 1: scores, keys = scores[:, None], [keys]
    if log: scores = np.log(scores)
    if labels is None:
        labels = adata.var_names if attr in {'var', 'varm'} else np.arange(scores.shape[0]).astype(str)
    if isinstance(labels, str):
        labels = [labels + str(i+1) for i in range(scores.shape[0])]
    from matplotlib import gridspec
    if n_panels <= 5: n_rows, n_cols = 1, n_panels
    else: n_rows, n_cols = 2, int(n_panels/2 + 0.5)
    fig = pl.figure(figsize=(n_cols * rcParams['figure.figsize'][0],
                             n_rows * rcParams['figure.figsize'][1]))
    left, bottom = 0.2/n_cols, 0.13/n_rows
    gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols, wspace=0.2,
                           left=left, bottom=bottom,
                           right=1-(n_cols-1)*left-0.01/n_cols,
                           top=1-(n_rows-1)*bottom-0.1/n_rows)
    for iscore, score in enumerate(scores.T):
        pl.subplot(gs[iscore])
        indices = np.argsort(score)[::-1][:n_points+1]
        for ig, g in enumerate(indices):
            pl.text(ig, score[g], labels[g], color=color,
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center', fontsize=8)
        pl.title(keys[iscore].replace('_', ' '))
        if n_panels <= 5 or count > n_cols: pl.xlabel('ranking')
        pl.xlim(-0.9, ig + 0.9)
        score_min, score_max = np.min(score[indices]), np.max(score[indices])
        pl.ylim((0.95 if score_min > 0 else 1.05) * score_min,
                (1.05 if score_max > 0 else 0.95) * score_max)
    return gs


def violin(adata, keys, group_by=None, use_raw=True, jitter=True, size=1, scale='width',
           order=None, multi_panel=False, show=None, save=None, ax=None):
    """Violin plot.

    Wraps `seaborn.violinplot` for :class:`~scanpy.api.AnnData`.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    keys : str
        Keys for accessing variables of adata.X or fields of adata.obs.
    group_by : str, optional
        The key of the sample grouping to consider.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    multi_panel : bool, optional
        Show fields in multiple panels. Returns a seaborn FacetGrid in that case.
    jitter : `float` or `bool`, optional (default: `True`)
        See `seaborn.stripplot`.
    size : int, optional (default: 1)
        Size of the jitter points.
    order : list of str, optional (default: `True`)
        Order in which to show the categories.
    scale : {'area', 'count', 'width'}, optional (default: 'width')
        The method used to scale the width of each violin. If 'area', each
        violin will have the same area. If 'count', the width of the violins
        will be scaled by the number of observations in that bin. If 'width',
        each violin will have the same width.
    show : bool, optional (default: `None`)
         Show the plot.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    ax : `matplotlib.Axes`
         A `matplotlib.Axes` object.

    Returns
    -------
    A `matplotlib.Axes` object if `ax` is `None` else `None`.
    """
    sanitize_anndata(adata)
    if group_by is not None and isinstance(keys, list):
        raise ValueError('Pass a single key as string if using `group_by`.')
    if isinstance(keys, str): keys = [keys]
    obs_keys = False
    for key in keys:
        if key in adata.obs_keys(): obs_keys = True
        if obs_keys and key not in set(adata.obs_keys()):
            raise ValueError(
                'Either use sample keys or variable names, but do not mix. '
                'Did not find {} in adata.obs_keys().'.format(key))
    if obs_keys:
        obs_df = adata.obs
    else:
        if group_by is None: obs_df = pd.DataFrame()
        else: obs_df = adata.obs.copy()
        for key in keys:
            if adata.raw is not None and use_raw:
                X_col = adata.raw[:, key].X
            else:
                X_col = adata[:, key].X
            obs_df[key] = X_col
    if group_by is None:
        obs_tidy = pd.melt(obs_df, value_vars=keys)
        x = 'variable'
        y = 'value'
    else:
        obs_tidy = obs_df
        x = group_by
        y = keys[0]
    if multi_panel:
        sns.set_style('whitegrid')
        g = sns.FacetGrid(obs_tidy, col=x, sharey=False)
        g = g.map(sns.violinplot, y, inner=None, orient='vertical', scale=scale)
        g = g.map(sns.stripplot, y, orient='vertical', jitter=jitter, size=size,
                     color='black').set_titles(
                         col_template='{col_name}').set_xlabels('')
        ax = g
    else:
        ax = sns.violinplot(x=x, y=y, data=obs_tidy, inner=None, order=order,
                            orient='vertical', scale=scale, ax=ax)
        ax = sns.stripplot(x=x, y=y, data=obs_tidy, order=order,
                           jitter=jitter, color='black', size=size, ax=ax)
        ax.set_xlabel('' if group_by is None else group_by.replace('_', ' '))
    utils.savefig_or_show('violin', show=show, save=save)
    return ax


def clustermap(
        adata, obs_keys=None, use_raw=True, show=None, save=None, **kwargs):
    """Hierarchically-clustered heatmap [Waskom16]_.

    Wraps `seaborn.clustermap <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_ for :class:`~scanpy.api.AnnData`.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    obs_keys : `str`
        Categorical annotation to plot with a different color map.
        Currently, only a single key is supported.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    show : bool, optional (default: `None`)
         Show the plot.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on \{'.pdf', '.png', '.svg'\}.
    **kwargs : keyword arguments
        Keyword arguments passed to `seaborn.clustermap <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_.

    Returns
    -------
    If `not show`, a `seaborn.ClusterGrid` object.

    Notes
    -----
    The returned object has a savefig() method that should be used if you want
    to save the figure object without clipping the dendrograms.

    To access the reordered row indices, use:
    clustergrid.dendrogram_row.reordered_ind

    Column indices, use: clustergrid.dendrogram_col.reordered_ind

    Examples
    --------
    Soon to come with figures. In the meanwile, see
    https://seaborn.pydata.org/generated/seaborn.clustermap.html.

    >>> import scanpy.api as sc
    >>> adata = sc.datasets.krumsiek11()
    >>> sc.pl.clustermap(adata, obs_keys='cell_type')
    """
    if not isinstance(obs_keys, (str, None)):
        raise ValueError('Currently, only a single key is supported.')
    sanitize_anndata(adata)
    X = adata.raw.X if use_raw and adata.raw is not None else adata.X
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    if obs_keys is not None:
        row_colors = adata.obs[obs_keys]
        utils.add_colors_for_categorical_sample_annotation(adata, obs_keys)
        # do this more efficiently... just a quick solution
        lut = dict(zip(
            row_colors.cat.categories,
            adata.uns[obs_keys + '_colors']))
        row_colors = adata.obs[obs_keys].map(lut)
        g = sns.clustermap(df, row_colors=row_colors, **kwargs)
    else:
        g = sns.clustermap(df, **kwargs)
    show = settings.autoshow if show is None else show
    if show: pl.show()
    else: return g
