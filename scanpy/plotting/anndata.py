"""Plotting functions for AnnData.
"""

import os
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import is_color_like
import seaborn as sns

from .. import settings
from .. import logging as logg
from . import utils
from .utils import scatter_base, scatter_group, setup_axes
from ..utils import sanitize_anndata, doc_params
from .utils import doc_scatter_bulk, doc_show_save_ax

VALID_LEGENDLOCS = {
    'none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left',
    'lower left', 'lower right', 'right', 'center left', 'center right',
    'lower center', 'upper center', 'center'
}


@doc_params(scatter_bulk=doc_scatter_bulk, show_save_ax=doc_show_save_ax)
def scatter(
        adata,
        x=None,
        y=None,
        color=None,
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
    """\
    Scatter plot along observations or variables axes.

    Color the plot using annotations of observations (`.obs`), variables
    (`.var`) or expression of genes (`.var_names`).

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    x : `str` or `None`
        x coordinate.
    y : `str` or `None`
        y coordinate.
    color : string or list of strings, optional (default: `None`)
        Keys for observation/cell or variable/gene annotation
        `[\'ann1\', \'ann2\']`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    basis : {{'pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr', etc.}}
        String that denotes a plotting tool that computed coordinates.
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a `matplotlib.Axis` or a list of it.
    """
    if basis is not None:
        axs = _scatter_obs(
            adata=adata,
            x=x,
            y=y,
            color=color,
            use_raw=use_raw,
            sort_order=sort_order,
            alpha=alpha,
            basis=basis,
            groups=groups,
            components=components,
            projection=projection,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            legend_fontweight=legend_fontweight,
            color_map=color_map,
            palette=palette,
            right_margin=right_margin,
            left_margin=left_margin,
            size=size,
            title=title,
            show=show,
            save=save,
            ax=ax)

    elif x is not None and y is not None:
        if x in adata.obs_keys() and y in adata.obs_keys() and color not in adata.var_keys():
            axs = _scatter_obs(
                adata=adata,
                x=x,
                y=y,
                color=color,
                use_raw=use_raw,
                sort_order=sort_order,
                alpha=alpha,
                basis=basis,
                groups=groups,
                components=components,
                projection=projection,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                legend_fontweight=legend_fontweight,
                color_map=color_map,
                palette=palette,
                right_margin=right_margin,
                left_margin=left_margin,
                size=size,
                title=title,
                show=show,
                save=save,
                ax=ax)

        elif x in adata.var_keys() and y in adata.var_keys() and color not in adata.obs_keys():
            axs = _scatter_var(
                adata=adata,
                x=x,
                y=y,
                color=color,
                use_raw=use_raw,
                sort_order=sort_order,
                alpha=alpha,
                basis=basis,
                groups=groups,
                components=components,
                projection=projection,
                legend_loc=legend_loc,
                legend_fontsize=legend_fontsize,
                legend_fontweight=legend_fontweight,
                color_map=color_map,
                palette=palette,
                right_margin=right_margin,
                left_margin=left_margin,
                size=size,
                title=title,
                show=show,
                save=save,
                ax=ax)
        else:
            raise ValueError(
                '`x`, `y`, and potential `color` inputs must all come from either `.obs` or `.var`')
    else:
        raise ValueError('Either provide a `basis` or `x` and `y`.')

    return axs


def _scatter_var(
        adata,
        x=None,
        y=None,
        color=None,
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

    adata_T = adata.T

    axs = _scatter_obs(
        adata=adata_T,
        x=x,
        y=y,
        color=color,
        use_raw=use_raw,
        sort_order=sort_order,
        alpha=alpha,
        basis=basis,
        groups=groups,
        components=components,
        projection=projection,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        color_map=color_map,
        palette=palette,
        right_margin=right_margin,
        left_margin=left_margin,
        size=size,
        title=title,
        show=show,
        save=save,
        ax=ax)

    # store .uns annotations that were added to the new adata object
    adata.uns = adata_T.uns

    return axs


def _scatter_obs(
        adata,
        x=None,
        y=None,
        color=None,
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
    """See docstring of scatter."""
    sanitize_anndata(adata)
    if legend_loc not in VALID_LEGENDLOCS:
        raise ValueError(
            'Invalid `legend_loc`, need to be one of: {}.'.format(VALID_LEGENDLOCS))
    if components is None: components = '1,2' if '2d' in projection else '1,2,3'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1
    keys = ['grey'] if color is None else [color] if isinstance(color, str) else color
    if title is not None and isinstance(title, str):
        title = [title]
    highlights = adata.uns['highlights'] if 'highlights' in adata.uns else []
    if basis is not None:
        try:
            # ignore the '0th' diffusion component
            if basis == 'diffmap': components += 1
            Y = adata.obsm['X_' + basis][:, components]
            # correct the component vector for use in labeling etc.
            if basis == 'diffmap': components -= 1
        except KeyError:
            raise KeyError('compute coordinates using visualization tool {} first'
                           .format(basis))
    elif x is not None and y is not None:
        x_arr = adata._get_obs_array(x)
        y_arr = adata._get_obs_array(y)
        Y = np.c_[x_arr[:, None], y_arr[:, None]]
    else:
        raise ValueError('Either provide a `basis` or `x` and `y`.')

    if size is None:
        n = Y.shape[0]
        size = 120000 / n

    if legend_loc.startswith('on data') and legend_fontsize is None:
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
        component_name = (
            'DC' if basis == 'diffmap'
            else 'tSNE' if basis == 'tsne'
            else 'UMAP' if basis == 'umap'
            else 'PC' if basis == 'pca'
            else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis
            else basis)
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
                raise ValueError(
                    'key \'{}\' is invalid! pass valid observation annotation, '
                    'one of {} or a gene name {}'
                    .format(key, adata.obs_keys(), adata.var_names))
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
        utils.add_colors_for_categorical_sample_annotation(
            adata, key, palette, force_update_colors=not palette_was_none)
        # actually plot the groups
        mask_remaining = np.ones(Y.shape[0], dtype=bool)
        centroids = {}
        if groups is None:
            for iname, name in enumerate(adata.obs[key].cat.categories):
                if name not in settings.categories_to_ignore:
                    mask = scatter_group(axs[ikey], key, iname,
                                         adata, Y, projection, size=size, alpha=alpha)
                    mask_remaining[mask] = False
                    if legend_loc.startswith('on data'): add_centroid(centroids, name, Y, mask)
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
                    if legend_loc.startswith('on data'): add_centroid(centroids, name, Y, mask)
                    mask_remaining[mask] = False
        if mask_remaining.sum() > 0:
            data = [Y[mask_remaining, 0], Y[mask_remaining, 1]]
            if projection == '3d': data.append(Y[mask_remaining, 2])
            axs[ikey].scatter(*data, marker='.', c='lightgrey', s=size,
                                    edgecolors='none', zorder=-1)
        legend = None
        if legend_loc.startswith('on data'):
            if legend_fontweight is None:
                legend_fontweight = 'bold'
            for name, pos in centroids.items():
                axs[ikey].text(pos[0], pos[1], name,
                               weight=legend_fontweight,
                               verticalalignment='center',
                               horizontalalignment='center',
                               fontsize=legend_fontsize)
            all_pos = np.zeros((len(adata.obs[key].cat.categories), 2))
            for iname, name in enumerate(adata.obs[key].cat.categories):
                if name in centroids:
                    all_pos[iname] = centroids[name]
                else:
                    all_pos[iname] = [np.nan, np.nan]
            utils._tmp_cluster_pos = all_pos
            if legend_loc == 'on data export':
                filename = settings.writedir + 'pos.csv'
                logg.msg('exporting label positions to {}'.format(filename), v=1)
                if settings.writedir != '' and not os.path.exists(settings.writedir):
                    os.makedirs(settings.writedir)
                np.savetxt(filename, all_pos, delimiter=',')
        elif legend_loc == 'right margin':
            legend = axs[ikey].legend(
                frameon=False, loc='center left',
                bbox_to_anchor=(1, 0.5),
                ncol=(1 if len(adata.obs[key].cat.categories) <= 14
                      else 2 if len(adata.obs[key].cat.categories) <= 30 else 3),
                fontsize=legend_fontsize)
        elif legend_loc != 'none':
            legend = axs[ikey].legend(
                frameon=False, loc=legend_loc, fontsize=legend_fontsize)
        if legend is not None:
            for handle in legend.legendHandles: handle.set_sizes([300.0])
    utils.savefig_or_show('scatter' if basis is None else basis, show=show, save=save)
    if show == False: return axs if len(keys) > 1 else axs[0]


def ranking(adata, attr, keys, dictionary=None, indices=None,
            labels=None, color='black', n_points=30,
            log=False, show=None):
    """Plot rankings.

    See, for example, how this is used in pl.pca_ranking.

    Parameters
    ----------
    adata : AnnData
        The data.
    attr : {'var', 'obs', 'uns', 'varm', 'obsm'}
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
        if dictionary is None:
            scores = getattr(adata, attr)[keys]
        else:
            scores = getattr(adata, attr)[dictionary][keys]
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
    if show == False: return gs


@doc_params(show_save_ax=doc_show_save_ax)
def violin(adata, keys, groupby=None, log=False, use_raw=True, stripplot=True, jitter=True,
           size=1, scale='width', order=None, multi_panel=None, show=None,
           xlabel='', rotation=None, save=None, ax=None, multi_panel_figsize=None,
           multi_panel_swap_axes=False, **kwargs):
    """\
    Violin plot [Waskom16]_.

    Wraps `seaborn.violinplot` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    keys : `str` or list of `str`
        Keys for accessing variables of `.var_names` or fields of `.obs`.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider.
    log : `bool`, optional (default: `False`)
        Plot on logarithmic axis.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    multi_panel : `bool`, optional (default: `False`)
        Display keys in multiple panels also when `groupby is not None`.
    stripplot : `bool` optional (default: `True`)
        Add a stripplot on top of the violin plot.
        See `seaborn.stripplot`.
    jitter : `float` or `bool`, optional (default: `True`)
        Add jitter to the stripplot (only when stripplot is True)
        See `seaborn.stripplot`.
    size : int, optional (default: 1)
        Size of the jitter points.
    order : list of str, optional (default: `True`)
        Order in which to show the categories.
    scale : {{'area', 'count', 'width'}}, optional (default: 'width')
        The method used to scale the width of each violin. If 'area', each
        violin will have the same area. If 'count', the width of the violins
        will be scaled by the number of observations in that bin. If 'width',
        each violin will have the same width.
    xlabel : `str`, optional (default: `''`)
        Label of the x axis. Defaults to `groupby` if `rotation` is `None`,
        otherwise, no label is shown.
    rotation : `float`, optional (default: `None`)
        Rotation of xtick labels.
    multi_panel_figsize : (float, float), optional (default: None)
        Figure size when multi_panel = True. Otherwise the rcParam['figure.figsize] value is used.
        Format is (width, height)
    multi_panel_swap_axes: `bool`, optional (default: `False`)
         By default, in multi_panel, the y axis contains the `keys` and the x axis the group by categories.
         By setting `multi_panel_swap_axes` then y are the group categories and x the `keys`.
    {show_save_ax}
    **kwargs : keyword arguments
        Are passed to `seaborn.violinplot`.

    Returns
    -------
    A `matplotlib.Axes` object if `ax` is `None` else `None`.
    """
    sanitize_anndata(adata)
    if isinstance(keys, str): keys = [keys]
    obs_keys = False
    for key in keys:
        if key in adata.obs_keys(): obs_keys = True
        if obs_keys and key not in set(adata.obs_keys()):
            raise ValueError(
                'Either use observation keys or variable names, but do not mix. '
                'Did not find {} in adata.obs_keys().'.format(key))
    if obs_keys:
        obs_df = adata.obs
    else:
        if groupby is None: obs_df = pd.DataFrame()
        else: obs_df = pd.DataFrame(adata.obs[groupby])
        for key in keys:
            if adata.raw is not None and use_raw:
                X_col = adata.raw[:, key].X
            else:
                X_col = adata[:, key].X
            obs_df[key] = X_col
    if groupby is None:
        obs_tidy = pd.melt(obs_df, value_vars=keys)
        x = 'variable'
        ys = ['value']
    else:
        obs_tidy = obs_df
        x = groupby
        ys = keys
    if multi_panel:
        if groupby is None and len(ys) == 1:
            # This is a quick and dirty way for adapting scales across several
            # keys if groupby is None.
            y = ys[0]
            g = sns.FacetGrid(obs_tidy, col=x, col_order=keys, sharey=False)
            # don't really know why this gives a warning without passing `order`
            g = g.map(sns.violinplot, y, inner=None, orient='vertical',
                      scale=scale, order=keys, **kwargs)
            if stripplot:
                g = g.map(sns.stripplot, y, orient='vertical', jitter=jitter, size=size, order=keys,
                          color='black')
            if log:
                g.set(yscale='log')
            g.set_titles(col_template='{col_name}').set_xlabels('')
            if rotation is not None:
                for ax in g.axes[0]:
                    ax.tick_params(labelrotation=rotation)
        else:
            # Make a very compact plot in which the y and x axis are shared.
            # The image is composed of individual plots stacked on top of each
            # other. Each subplot contains and individual violin plot where
            # x = categories in `groupby` and y is each of the keys provided.
            # If multi_panel_swap_axes is True, then x and y are swapped.
            # An example is: keys = marker genes, groupby = louvain clusters.

            categories = adata.obs[groupby].cat.categories
            if multi_panel_swap_axes is False:
                if multi_panel_figsize is None:
                    height = len(ys) * 0.6 + 3
                    width = len(categories) * 0.2 + 1
                else:
                    width, height = multi_panel_figsize
                fig, axs = pl.subplots(nrows=len(ys), ncols=1, sharex=True, sharey=True,
                                       figsize=(width, height))
                for idx, y in enumerate(ys):
                    if len(ys) > 1:
                        ax = axs[idx]
                    else:
                        ax = axs

                    ax = sns.violinplot(x, y=y, data=obs_tidy, inner=None, order=order,
                                        orient='vertical', scale=scale, ax=ax, **kwargs)
                    if stripplot:
                        ax = sns.stripplot(x, y=y, data=obs_tidy, order=order,
                                           jitter=jitter, color='black', size=size, ax=ax)

                    ax.set_ylabel(y, rotation=0, fontsize=11, labelpad=8, ha='right')
                    # remove the grids because in such a compact plot are unnecessary
                    ax.grid(False)
                    ax.tick_params(axis='y', left=False, right=True, labelright=True, labelleft=False)

                    # remove the xticks labels except for the last processed plot (first from bottom-up).
                    # Because the plots share the x axis it is redundant and less compact to plot the
                    # axis for each plot
                    if idx < len(ys) - 1:
                        ax.set_xticklabels([])
                    if log:
                        ax.set_yscale('log')
                    if rotation is not None:
                        ax.tick_params(labelrotation=rotation)
            else:
                if multi_panel_figsize is None:
                    height = len(categories) * 0.6 + 3
                    width = len(ys) * 0.2 + 1
                else:
                    height, width = multi_panel_figsize
                fig, axs = pl.subplots(nrows=len(categories), ncols=1, sharex=True, sharey=True,
                                       figsize=(width, height))
                for idx, category in enumerate(categories):
                    df = pd.melt(obs_tidy[obs_tidy[groupby] == category], value_vars=keys)
                    ax = axs[idx]
                    ax = sns.violinplot('variable', y='value', data=df, inner=None, order=order,
                                        orient='vertical', scale=scale, ax=ax, **kwargs)

                    if stripplot:
                        ax = sns.stripplot('variable', y='value', data=df, order=order,
                                           jitter=jitter, color='black', size=size, ax=ax)

                    ax.set_ylabel(category, rotation=0, fontsize=11, labelpad=8, ha='right')
                    # remove the grids because in such a compact plot are unnecessary
                    ax.grid(False)
                    ax.tick_params(axis='y', left=False, right=True, labelright=True, labelleft=False)

                    # remove the xticks labels except for the last processed plot (first from bottom-up).
                    # Because the plots share the x axis it is redundant and less compact to plot the
                    # axis for each plot
                    if idx < len(categories) - 1:
                        ax.set_xticklabels([])
                    if log:
                        ax.set_yscale('log')
                    if rotation is not None:
                        ax.tick_params(labelrotation=rotation)
            # remove the spacing between subplots
            pl.subplots_adjust(wspace=0, hspace=0)

    else:
        if ax is None:
            axs, _, _, _ = setup_axes(
                ax=ax, panels=['x'] if groupby is None else keys, show_ticks=True, right_margin=0.3)
        for ax, y in zip(axs, ys):
            ax = sns.violinplot(x, y=y, data=obs_tidy, inner=None, order=order,
                                orient='vertical', scale=scale, ax=ax, **kwargs)
            if stripplot:
                ax = sns.stripplot(x, y=y, data=obs_tidy, order=order,
                                   jitter=jitter, color='black', size=size, ax=ax)
            if xlabel == '' and groupby is not None and rotation is None:
                xlabel = groupby.replace('_', ' ')
            ax.set_xlabel(xlabel)
            if log:
                ax.set_yscale('log')
            if rotation is not None:
                ax.tick_params(labelrotation=rotation)
    utils.savefig_or_show('violin', show=show, save=save)
    if show == False: return axs[0] if len(axs) == 1 else axs


@doc_params(show_save_ax=doc_show_save_ax)
def clustermap(
        adata, obs_keys=None, use_raw=True, show=None, save=None, **kwargs):
    """\
    Hierarchically-clustered heatmap [Waskom16]_.

    Wraps `seaborn.clustermap
    <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`__ for
    :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_keys : `str`
        Categorical annotation to plot with a different color map.
        Currently, only a single key is supported.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {show_save_ax}
    **kwargs : keyword arguments
        Keyword arguments passed to `seaborn.clustermap
        <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`__.

    Returns
    -------
    If `show == False`, a `seaborn.ClusterGrid` object.

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
    if not isinstance(obs_keys, (str, type(None))):
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


@doc_params(show_save_ax=doc_show_save_ax)
def heatmap(adata, var_names, groupby=None, use_raw=True, log=False, num_categories=7,
            show=None, save=None, figsize=None, **kwargs):
    """\
    Heatmap of the expression values of set of genes..

    If `groupby` is given, the heatmap is ordered by the respective group. For
    example, a list of marker genes can be plotted, ordered by clustering. If
    the `groupby` observation annotation is not categorical the observation
    annotation is turned into a categorical by binning the data into the number
    especified in `num_categories`.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names : `str` or list of `str`
        `var_names` should be a valid subset of  `adata.var_names`.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories`.
    log : `bool`, optional (default: `False`)
        Use the log of the values
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    figsize : (float, float), optional (default: None)
        Figure size (width, height). If not set, the figure width is set based on the
        number of  `var_names` and the height is set to 10.
    {show_save_ax}
    **kwargs : keyword arguments
        Are passed to `seaborn.heatmap`.

    Returns
    -------
    A list of `matplotlib.Axes` where the first ax is the groupby categories
    colorcode, the second axis is the heatmap and the third axis is the
    colorbar.
    """
    from scipy.sparse import issparse
    sanitize_anndata(adata)
    if isinstance(var_names, str):
        var_names = [var_names]
    if groupby is not None:
        if groupby not in adata.obs_keys():
            raise ValueError('groupby has to be a valid observation. Given value: {}, '
                             'valid observations: {}'.format(groupby, adata.obs_keys()))

    if use_raw:
        matrix = adata.raw[:, var_names].X
    else:
        matrix = adata[:, var_names].X

    if issparse(matrix):
        matrix = matrix.toarray()
    if log:
        matrix = np.log1p(matrix)

    obs_tidy = pd.DataFrame(matrix, columns=var_names)
    if groupby is None:
        groupby = ''
        categorical = pd.Series(np.repeat('', len(obs_tidy))).astype('category')
    else:
        if not is_categorical_dtype(adata.obs[groupby]):
            # if the groupby column is not categorical, turn it into one
            # by subdividing into 7 categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]
    obs_tidy.set_index(categorical, groupby, inplace=True)
    categories = obs_tidy.index.categories

    if figsize is None:
        height = 10
        heatmap_width = len(var_names) * 0.18
        width = heatmap_width + 3  # +3 to account for the colorbar and labels
    else:
        width, height = figsize
    ax_frac2width = 0.25
    fig, axs = pl.subplots(
        nrows=1, ncols=3, sharey=False,
        figsize=(width, height),
        gridspec_kw={'width_ratios': [ax_frac2width, width, ax_frac2width]})
    groupby_ax = axs[0]
    heatmap_ax = axs[1]
    heatmap_cbar_ax = axs[2]

    if groupby:
        obs_tidy = obs_tidy.sort_index()

    # determine groupby label positions
    value_sum = 0
    ticks = []
    labels = []
    label2code = {}
    for code, (label, value) in enumerate(obs_tidy.index.value_counts(sort=False).iteritems()):
        ticks.append(value_sum + (value / 2))
        labels.append(label)
        value_sum += value
        label2code[label] = code

    groupby_ax.imshow(np.matrix([label2code[lab] for lab in obs_tidy.index]).T, aspect='auto')
    if len(categories) > 1:
        groupby_ax.set_yticks(ticks)
        groupby_ax.set_yticklabels(labels)

    # remove y ticks
    groupby_ax.tick_params(axis='y', left=False)
    # remove x ticks and labels
    groupby_ax.tick_params(axis='x', bottom=False, labelbottom=False)

    # remove surrounding lines
    groupby_ax.spines['right'].set_visible(False)
    groupby_ax.spines['top'].set_visible(False)
    groupby_ax.spines['left'].set_visible(False)
    groupby_ax.spines['bottom'].set_visible(False)

    groupby_ax.set_ylabel(groupby)
    groupby_ax.grid(False)

    sns.heatmap(obs_tidy, yticklabels='none', ax=heatmap_ax, cbar_ax=heatmap_cbar_ax, **kwargs)
    heatmap_ax.tick_params(axis='y', left=False, labelleft=False)
    heatmap_ax.set_ylabel('')
    heatmap_ax.set_xticks(np.arange(len(var_names)) + 0.5)
    heatmap_ax.set_xticklabels(var_names)
    pl.subplots_adjust(wspace=0.03, hspace=0.01)
    utils.savefig_or_show('heatmap', show=show, save=save)

    return axs


def dotplot(adata, var_names, groupby=None, use_raw=True, log=False, num_categories=7,
            figsize=None, show=None, save=None, **kwargs):
    """Makes a 'dot plot' of the expression values of `var_names`.
    For each var_name and each groupby category a dot is plotted.
    Each dot represents two values: mean expression within
    each category (visualized by color) and fraction of cells expressing the
    var_name in the category. (visualized by the size of the dot).
    If groupby is not given, the dotplot assumes that all data belongs to a single
    category. A gene is not considered expressed if the expression value in the adata
    (or adata.raw) is equal to zero.
    For example, for each marker gene, the mean value and the percentage of cells
    expressing the gene can be visualized for each cluster.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    var_names : `str` or list of `str`
        var_names should be a valid subset of  `.var_names`.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. It is expected that groupby is
        a categorical. If groupby is not a categorical observation, it would be
        subdivided into `num_categories`.
    log : `bool`, optional (default: `False`)
        Use the log of the values
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value determines
        the number of groups into which the groupby observation should be subdivided.
    figsize : (float, float), optional (default: None)
        Figure size (width, height. If not set, the figure width is set based on the
        number of  `var_names` and the height is set to 10.
    show : bool, optional (default: `None`)
         Show the plot.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the
        default filename. Infer the filetype if ending on {{'.pdf', '.png', '.svg'}}.
    **kwargs : keyword arguments
        Are passed to `matplotlib.pyplot.scatter`.

    Returns
    -------
    A list of `matplotlib.Axes` where the first ax is the groupby categories colorcode, the
    second axis is the heatmap and the third axis is the colorbar.

    """
    from scipy.sparse import issparse
    sanitize_anndata(adata)
    if isinstance(var_names, str):
        var_names = [var_names]
    if groupby is not None:
        if groupby not in adata.obs_keys():
            raise ValueError('groupby has to be a valid observation. Given value: {}, '
                             'valid observations: {}'.format(groupby, adata.obs_keys()))

    if use_raw:
        matrix = adata.raw[:, var_names].X
    else:
        matrix = adata[:, var_names].X

    if issparse(matrix):
        matrix = matrix.toarray()
    if log:
        matrix = np.log1p(matrix)

    obs_tidy = pd.DataFrame(matrix, columns=var_names)
    if groupby is None:
        # assign all data to a single category to allow groupby opeerations
        categorical = pd.Series(np.repeat(0, obs_tidy.shape[0]), dtype='category')
    else:
        if not is_categorical_dtype(adata.obs[groupby]):
            # if the groupby column is not categorical, turn it into one
            # by subdividing into 7 categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]
    obs_tidy.set_index(categorical, groupby, inplace=True)
    categories = obs_tidy.index.categories
    # for if category defined by groupby (if any) compute for each var_name
    # 1. the mean value over the category
    # 2. the fraction of cells in the category having a value > 0

    # 1. compute mean value
    mean_obs = obs_tidy.groupby(level=0).mean()

    # 2. compute fraction of cells having value >0
    # transform obs_tidy into boolean matrix
    obs_bool = obs_tidy.astype(bool)

    # compute the sum per group which in the boolean matrix this is the number
    # of values >0, and divide the result by the total number of values in the group
    # (given by `count()`)
    fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()

    if figsize is None:
        height = len(categories) * 0.3 + 1  # +1 for labels
        # if the number of categories is small (eg 1 or 2) use
        # a larger height
        height = max([1.5, height])
        heatmap_width = len(var_names) * 0.5
        width = heatmap_width + 1.6 + 1  # +1.6 to account for the colorbar and  + 1 to account for labels
    else:
        width, height = figsize
        heatmap_width = width * 0.75

    # colorbar ax width should not change with differences in the width of the image
    # otherwise can become too small
    colorbar_width = 0.4
    colorbar_width_spacer = 0.7
    size_legend_width = 0.5

    # define a layout of 1 row x 4 columns
    # First ax is for the main figure
    # Second ax is for the color bar legend
    # Third ax is for an spacer that avoids the ticks
    # from the color bar to be hidden beneath the next axis
    # Fourth ax is to plot the size legend
    from matplotlib import gridspec
    fig = pl.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=1, ncols=4, left=0.05, right=0.48, wspace=0.05,
                            width_ratios=[heatmap_width, colorbar_width, colorbar_width_spacer, size_legend_width])
    if len(categories) < 4:
        #  hen two few categories are shown, the colorbar and size legend
        # are compressed, thus, the dotplot ax is split into two:
        axs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[0], height_ratios=[len(categories) * 0.3, 1])
        dot_ax = fig.add_subplot(axs2[0])
    else:
        dot_ax = fig.add_subplot(axs[0])

    color_legend = fig.add_subplot(axs[1])

    # to keep the size_legen of about the same height, irrespective
    # of the number of categories, the fourth ax is subdivided into two parts
    size_legend_height = min(1.3, height)
    # wspace is proportional to the width but a constant value is
    # needed such that the spacing is the same for thinner or wider images.
    wspace = 10.5/ width
    axs3 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[3], wspace=wspace,
                                            height_ratios=[size_legend_height / height,
                                                           (height - size_legend_height) / height])
    size_legend = fig.add_subplot(axs3[0])

    # make scatter plot in which
    # x = var_names
    # y = groupby category
    # size = fraction
    # color = mean expression

    y, x = np.indices(mean_obs.shape)
    y = y.flatten()
    x = x.flatten()
    frac = fraction_obs.values.flatten()
    mean_flat = mean_obs.values.flatten()
    cmap = pl.get_cmap('Reds')

    size = (frac * 10) ** 2

    import matplotlib.colors
    normalize = matplotlib.colors.Normalize(vmin=min(mean_flat), vmax=max(mean_flat))
    colors = [cmap(normalize(value)) for value in mean_flat]

    dot_ax.scatter(x, y, color=colors, s=size, cmap=cmap, norm=None, edgecolor='none', **kwargs)
    y_ticks = range(mean_obs.shape[0])
    dot_ax.set_yticks(y_ticks)
    dot_ax.set_yticklabels([mean_obs.index[idx] for idx in y_ticks])

    x_ticks = range(mean_obs.shape[1])
    dot_ax.set_xticks(x_ticks)
    dot_ax.set_xticklabels([mean_obs.columns[idx] for idx in x_ticks], rotation=90)
    dot_ax.grid(False)
    dot_ax.set_xlim(-0.5, len(var_names) + 0.5)
    dot_ax.set_ylabel(groupby)

    # to be consistent with the heatmap plot, is better to
    # invert the order of the y-axis, such that the first group is on
    # top
    ymin, ymax = dot_ax.get_ylim()
    dot_ax.set_ylim(ymax+0.5, ymin - 0.5)

    dot_ax.set_xlim(-1, len(var_names) + 0.5)

    # plot colorbar
    import matplotlib.colorbar
    matplotlib.colorbar.ColorbarBase(color_legend, cmap=cmap, norm=normalize)

    # plot size bar
    fracs_legend = np.array([0.25, 0.50, 0.75, 1])
    size = (fracs_legend * 10) ** 2
    color = [cmap(normalize(value)) for value in np.repeat(max(mean_flat) * 0.7, len(size))]
    size_legend.scatter(np.repeat(0, len(size)), range(len(size)), s=size, color=color)
    size_legend.set_yticks(range(len(size)))
    size_legend.set_yticklabels(["{:.0%}".format(x) for x in fracs_legend])

    size_legend.tick_params(axis='y', left=False, labelleft=False, labelright=True)

    # remove x ticks and labels
    size_legend.tick_params(axis='x', bottom=False, labelbottom=False)

    # remove surrounding lines
    size_legend.spines['right'].set_visible(False)
    size_legend.spines['top'].set_visible(False)
    size_legend.spines['left'].set_visible(False)
    size_legend.spines['bottom'].set_visible(False)
    size_legend.grid(False)

    ymin, ymax = size_legend.get_ylim()
    size_legend.set_ylim(ymin, ymax+0.5)

    utils.savefig_or_show('heatmap', show=show, save=save)
    return axs
