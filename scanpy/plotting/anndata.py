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
from .docs import doc_scatter_bulk, doc_show_save_ax

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
        use_raw=None,
        layers='X',
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
        frameon=None,
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
        Keys for annotations of observations/cells or variables/genes, e.g.,
        `'ann1'` or `['ann1', 'ann2']`.
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    layers : `str` or tuple of strings, optional (default: `X`)
        Use the `layers` attribute of `adata` if present: specify the layer for
        `x`, `y` and `color`. If `layers` is a string, then it is expanded to
        `(layers, layers, layers)`.
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
            layers=layers,
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
            frameon=frameon,
            right_margin=right_margin,
            left_margin=left_margin,
            size=size,
            title=title,
            show=show,
            save=save,
            ax=ax)
    elif x is not None and y is not None:
        if ((x in adata.obs.keys() or x in adata.var.index)
            and (y in adata.obs.keys() or y in adata.var.index)
            and (color is None or color in adata.obs.keys() or color in adata.var.index)):
            axs = _scatter_obs(
                adata=adata,
                x=x,
                y=y,
                color=color,
                use_raw=use_raw,
                layers=layers,
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
                frameon=frameon,
                right_margin=right_margin,
                left_margin=left_margin,
                size=size,
                title=title,
                show=show,
                save=save,
                ax=ax)
        elif ((x in adata.var.keys() or x in adata.obs.index)
                and (y in adata.var.keys() or y in adata.obs.index)
                and (color is None or color in adata.var.keys() or color in adata.obs.index)):
            axs = _scatter_var(
                adata=adata,
                x=x,
                y=y,
                color=color,
                use_raw=use_raw,
                layers=layers,
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
                frameon=frameon,
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
        use_raw=None,
        layers='X',
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
        frameon=None,
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
        layers=layers,
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
        frameon=frameon,
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
        use_raw=None,
        layers='X',
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
        frameon=None,
        right_margin=None,
        left_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """See docstring of scatter."""
    sanitize_anndata(adata)
    from scipy.sparse import issparse
    if use_raw is None and adata.raw is not None: use_raw = True

    # process layers
    if layers is None:
        layers = 'X'
    if isinstance(layers, str) and (layers == 'X' or layers in adata.layers.keys()):
        layers = (layers, layers, layers)
    elif isinstance(layers, (tuple, list)) and len(layers) == 3:
        for layer in layers:
            if layer not in adata.layers.keys() and layer != 'X':
                raise ValueError(
                    '`layers` should have elements that are either \'X\' or in adata.layers.keys().')
    else:
        raise ValueError('`layers` should be a string or a list/tuple of length 3.')
    if use_raw and (layers != ('X', 'X', 'X') or layers != ['X', 'X', 'X']):
        ValueError('`use_raw` must be `False` if layers other than \'X\' are used.')

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
        x_arr = adata._get_obs_array(x, use_raw=use_raw, layer=layers[0])
        y_arr = adata._get_obs_array(y, use_raw=use_raw, layer=layers[1])

        x_arr = x_arr.toarray().flatten() if issparse(x_arr) else x_arr
        y_arr = y_arr.toarray().flatten() if issparse(y_arr) else y_arr

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

    # generate the colors
    color_ids = []
    categoricals = []
    colorbars = []
    for ikey, key in enumerate(keys):
        c = 'white'
        categorical = False  # by default, assume continuous or flat color
        colorbar = None
        # test whether we have categorial or continuous annotation
        if key in adata.obs_keys():
            if is_categorical_dtype(adata.obs[key]):
                categorical = True
            else:
                c = adata.obs[key]
        # coloring according to gene expression
        elif (use_raw
              and adata.raw is not None
              and key in adata.raw.var_names):
            c = adata.raw[:, key].X
        elif key in adata.var_names:
            c = adata[:, key].X if layers[2] == 'X' else adata[:, key].layers[layers[2]]
            c = c.toarray().flatten() if issparse(c) else c
        elif is_color_like(key):  # a flat color
            c = key
            colorbar = False
        else:
            raise ValueError(
                'key \'{}\' is invalid! pass valid observation annotation, '
                'one of {} or a gene name {}'
                .format(key, adata.obs_keys(), adata.var_names))
        if colorbar is None:
            colorbar = not categorical
        colorbars.append(colorbar)
        if categorical: categoricals.append(ikey)
        color_ids.append(c)

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

    # loop over all categorical annotation and plot it
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
            groups = [groups] if isinstance(groups, str) else groups
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

    # draw a frame around the scatter
    frameon = settings._frameon if frameon is None else frameon
    if not frameon:
        for ax in axs:
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_frame_on(False)

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
def violin(adata, keys, groupby=None, log=False, use_raw=None, stripplot=True, jitter=True,
           size=1, scale='width', order=None, multi_panel=None, show=None,
           xlabel='', rotation=None, save=None, ax=None, **kwds):
    """\
    Violin plot.

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
    use_raw : `bool`, optional (default: `None`)
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
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `seaborn.violinplot`.

    Returns
    -------
    A `matplotlib.Axes` object if `ax` is `None` else `None`.
    """
    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None: use_raw = True
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
                      scale=scale, order=keys, **kwds)
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
        if ax is None:
            axs, _, _, _ = setup_axes(
                ax=ax, panels=['x'] if groupby is None else keys, show_ticks=True, right_margin=0.3)
        else:
            axs = [ax]
        for ax, y in zip(axs, ys):
            ax = sns.violinplot(x, y=y, data=obs_tidy, inner=None, order=order,
                                orient='vertical', scale=scale, ax=ax, **kwds)
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
        adata, obs_keys=None, use_raw=None, show=None, save=None, **kwds):
    """\
    Hierarchically-clustered heatmap.

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
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    {show_save_ax}
    **kwds : keyword arguments
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
    if use_raw is None and adata.raw is not None: use_raw = True
    X = adata.raw.X if use_raw else adata.X
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    if obs_keys is not None:
        row_colors = adata.obs[obs_keys]
        utils.add_colors_for_categorical_sample_annotation(adata, obs_keys)
        # do this more efficiently... just a quick solution
        lut = dict(zip(
            row_colors.cat.categories,
            adata.uns[obs_keys + '_colors']))
        row_colors = adata.obs[obs_keys].map(lut)
        g = sns.clustermap(df, row_colors=row_colors, **kwds)
    else:
        g = sns.clustermap(df, **kwds)
    show = settings.autoshow if show is None else show
    if show: pl.show()
    else: return g


@doc_params(show_save_ax=doc_show_save_ax)
def stacked_violin(adata, var_names, groupby=None, log=False, use_raw=None, num_categories=7,
                   stripplot=False, jitter=False, size=1, scale='width', order=None,
                   show=None, save=None, figsize=None, var_group_positions=None,
                   var_group_labels=None, var_group_rotation=None, swap_axes=False, **kwds):
    """\
    Stacked violin plots.

    Makes a compact image composed of individual violin plots (from `seaborn.violinplot`)
    stacked on top of each other. Useful to visualize gene expression per cluster.

    Wraps `seaborn.violinplot` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names : `str` or list of `str`
        `var_names` should be a valid subset of  `adata.var_names`.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider.
    log : `bool`, optional (default: `False`)
        Plot on logarithmic axis.
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
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
    figsize : (float, float), optional (default: None)
        Figure size when multi_panel = True. Otherwise the rcParam['figure.figsize] value is used.
        Format is (width, height)
    var_group_positions :  list of `tuples`.
        Use this parameter to highlight groups of `var_names` (only when swap_axes=False).
        This will draw a 'bracket'
        on top of the plot between the given start and end positions. If the
        parameter `var_group_labels` is set, the corresponding labels is added on
        top of the bracket. E.g. var_group_positions = [(4,10)] will add a bracket
        between the fourth var_name and the tenth var_name. By giving more
        positions, more brackets are drawn.
    var_group_labels : list of `str`
        Labels for each of the var_group_positions that want to be highlighted.
    swap_axes: `bool`, optional (default: `False`)
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby `categories.
         By setting `swap_axes` then y are the `groupby` categories and x the `var_names`. When swapping
         axes var_group_positions are no longer used
    var_group_rotation : `float` (default: `None`)
        Label rotation degrees. By default, labels larger than 4 characters are rotated 90 degrees
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `seaborn.violinplot`.

    Returns
    -------
    A list of `matplotlib.Axes` where each ax corresponds to each row in the image
    """
    if use_raw is None and adata.raw is not None: use_raw = True

    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories)
    from matplotlib import gridspec

    if swap_axes is False:
        # plot image in which x = var_names and y = groupby categories
        if figsize is None:
            height = len(categories) * 0.2 + 3
            width = len(var_names) * 0.2 + 1
        else:
            width, height = figsize

        num_rows = len(categories)
        height_ratios = None
        if var_group_positions is not None and len(var_group_positions) > 0:
            # add some space in case 'brackets' want to be plotted on top of the image
            num_rows += 2  # +2 to add the row for the brackets and a spacer
            height_ratios = [0.2, 0.05] + [float(height) / len(categories)] * len(categories)
            categories = [None, None] + list(categories)
        fig = pl.figure(figsize=(width, height))
        # define a layout of nrows = len(categories) rows x 1 column
        # each row is one violin plot.
        # if var_group_positions is defined, a new row is added
        axs = gridspec.GridSpec(nrows=num_rows, ncols=1, height_ratios=height_ratios)

        ax0 = None

        for idx in range(num_rows)[::-1]:  # iterate in reverse to start on the bottom plot
                                           # this facilitates adding the brackets plot (if
                                           # needed) by sharing the x axis with a previous
                                           # violing plot.
            category = categories[idx]
            if var_group_positions is not None and len(var_group_positions) > 0 and idx <= 1:
                # if var_group_positions is given, axs[0] and axs[1] are the location for the
                # brackets and a spacer (axs[1])
                if idx == 0:
                    brackets_ax = fig.add_subplot(axs[0], sharex=ax0)
                    _plot_gene_groups_brackets(brackets_ax, group_positions=var_group_positions,
                                               group_labels=var_group_labels,
                                               rotation=var_group_rotation)

                continue

            df = pd.melt(obs_tidy[obs_tidy.index == category], value_vars=var_names)

            if ax0 is None:
                ax = fig.add_subplot(axs[idx])
                ax0 = ax

            else:
                ax = fig.add_subplot(axs[idx], sharey=ax0)

            ax = sns.violinplot('variable', y='value', data=df, inner=None, order=order,
                                orient='vertical', scale=scale, ax=ax, **kwds)

            if stripplot:
                ax = sns.stripplot('variable', y='value', data=df, order=order,
                                   jitter=jitter, color='black', size=size, ax=ax)

            # remove the grids because in such a compact plot are unnecessary
            ax.grid(False)

            ax.tick_params(axis='y', left=False, right=True, labelright=True,
                           labelleft=False, labelsize='x-small')
            ax.set_ylabel(category, rotation=0, fontsize='small', labelpad=8, ha='right', va='center')
            ax.set_xlabel('')
            if log:
                ax.set_yscale('log')

            if idx < num_rows - 1:
                # remove the xticks labels except for the last processed plot (first from bottom-up).
                # Because the plots share the x axis it is redundant and less compact to plot the
                # axis ticks and labels for each plot
                ax.set_xticklabels([])
                ax.tick_params(axis='x', bottom=False, top=False, labeltop=False, labelbottom=False)

        ax0.tick_params(axis='x', labelrotation=90, labelsize='small')

    else:
        # plot image in which x = group by and y = var_names
        if figsize is None:
            height = len(var_names) * 0.2 + 3
            width = len(categories) * 0.2 + 1
        else:
            width, height = figsize
        fig, axs = pl.subplots(nrows=len(var_names), ncols=1, sharex=True, sharey=True,
                               figsize=(width, height))
        for idx, y in enumerate(var_names):
            if len(var_names) > 1:
                ax = axs[idx]
            else:
                ax = axs

            ax = sns.violinplot(x=obs_tidy.index, y=y, data=obs_tidy, inner=None, order=order,
                                orient='vertical', scale=scale, ax=ax, **kwds)
            if stripplot:
                ax = sns.stripplot(x=obs_tidy.index, y=y, data=obs_tidy, order=order,
                                   jitter=jitter, color='black', size=size, ax=ax)

            ax.set_ylabel(y, rotation=0, fontsize='small', labelpad=8, ha='right', va='center')
            # remove the grids because in such a compact plot are unnecessary
            ax.grid(False)
            ax.tick_params(axis='y', right=True, labelright=True,
                           labelleft=False, labelrotation=0, labelsize='x-small')
            ax.tick_params(axis='x', labelsize='small')

            # remove the xticks labels except for the last processed plot (first from bottom-up).
            # Because the plots share the x axis it is redundant and less compact to plot the
            # axis for each plot
            if idx < len(var_names) - 1:
                ax.set_xticklabels([])
            if log:
                ax.set_yscale('log')

            if max([len(x) for x in categories]) > 1:
                ax.tick_params(axis='x', labelrotation=90)

    # remove the spacing between subplots
    pl.subplots_adjust(wspace=0, hspace=0)

    utils.savefig_or_show('stacked_violin', show=show, save=save)

    return axs


@doc_params(show_save_ax=doc_show_save_ax)
def heatmap(adata, var_names, groupby=None, use_raw=None, log=False, num_categories=7,
            var_group_positions=None, var_group_labels=None,
            var_group_rotation=None, show=None, save=None, figsize=None, **kwds):
    """\
    Heatmap of the expression values of set of genes..

    If `groupby` is given, the heatmap is ordered by the respective group. For
    example, a list of marker genes can be plotted, ordered by clustering. If
    the `groupby` observation annotation is not categorical the observation
    annotation is turned into a categorical by binning the data into the number
    specified in `num_categories`.

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
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    figsize : (float, float), optional (default: None)
        Figure size (width, height). If not set, the figure width is set based on the
        number of  `var_names` and the height is set to 10.
    var_group_positions :  list of `tuples`.
        Use this parameter to highlight groups of `var_names`. This will draw a 'bracket'
        on top of the plot between the given start and end positions. If the
        parameter `var_group_labels` is set, the corresponding labels is added on
        top of the bracket. E.g. var_group_positions = [(4,10)] will add a bracket
        between the fourth var_name and the tenth var_name. By giving more
        positions, more brackets are drawn.
    var_group_labels : list of `str`
        Labels for each of the var_group_positions that want to be highlighted.
    var_group_rotation : `float` (default: `None`)
        Label rotation degrees. By default, labels larger than 4 characters are rotated 90 degrees
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `seaborn.heatmap`.

    Returns
    -------
    A list of `matplotlib.Axes` where the first ax is the groupby categories
    colorcode, the second axis is the heatmap and the third axis is the
    colorbar.
    """
    if use_raw is None and adata.raw is not None: use_raw = True
    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories)

    if figsize is None:
        height = 6
        heatmap_width = len(var_names) * 0.25
        width = heatmap_width + 3  # +3 to account for the colorbar and labels
    else:
        width, height = figsize
    ax_frac2width = 0.25

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.15, height]
    else:
        height_ratios = [0, height]

    # define a layout of 2 rows x 3 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
    # second row is for main content. This second row is divided into three axes:
    #   first ax is for the categories defined by `groupby`
    #   second ax is for the heatmap
    #   third ax is for colorbar

    from matplotlib import gridspec
    fig = pl.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=2, ncols=3, left=0.05, right=0.48, wspace=0.5 / width,
                            hspace=0.13 / height,
                            width_ratios=[ax_frac2width, width, ax_frac2width],
                            height_ratios=height_ratios)

    groupby_ax = fig.add_subplot(axs[1, 0])
    heatmap_ax = fig.add_subplot(axs[1, 1])
    heatmap_cbar_ax = fig.add_subplot(axs[1, 2])
    heatmap_cbar_ax.tick_params(axis='y', labelsize='small')

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
    groupby_ax.tick_params(axis='y', left=False, labelsize='small')
    # remove x ticks and labels
    groupby_ax.tick_params(axis='x', bottom=False, labelbottom=False)

    # remove surrounding lines
    groupby_ax.spines['right'].set_visible(False)
    groupby_ax.spines['top'].set_visible(False)
    groupby_ax.spines['left'].set_visible(False)
    groupby_ax.spines['bottom'].set_visible(False)

    groupby_ax.set_ylabel(groupby)
    groupby_ax.grid(False)

    sns.heatmap(obs_tidy, yticklabels='none', ax=heatmap_ax, cbar_ax=heatmap_cbar_ax, **kwds)
    heatmap_ax.tick_params(axis='y', left=False, labelleft=False)
    heatmap_ax.tick_params(axis='x', labelsize='small')
    heatmap_ax.set_ylabel('')
    heatmap_ax.set_xticks(np.arange(len(var_names)) + 0.5)
    heatmap_ax.set_xticklabels(var_names)

    # plot group legends on top of heatmap_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 1], sharex=heatmap_ax)
        _plot_gene_groups_brackets(gene_groups_ax, group_positions=var_group_positions,
                                   group_labels=var_group_labels, rotation=var_group_rotation,
                                   left_adjustment=0.2, right_adjustment=0.8)

    utils.savefig_or_show('heatmap', show=show, save=save)

    return axs


@doc_params(show_save_ax=doc_show_save_ax)
def dotplot(adata, var_names, groupby=None, use_raw=None, log=False, num_categories=7,
            color_map='Reds', figsize=None, var_group_positions=None, var_group_labels=None,
            var_group_rotation=None, show=None, save=None, **kwds):
    """\
    Makes a _dot plot_ of the expression values of `var_names`.

    For each var_name and each `groupby` category a dot is plotted. Each dot
    represents two values: mean expression within each category (visualized by
    color) and fraction of cells expressing the var_name in the
    category. (visualized by the size of the dot).  If groupby is not given, the
    dotplot assumes that all data belongs to a single category. A gene is not
    considered expressed if the expression value in the adata (or adata.raw) is
    equal to zero.

    For instance, for each marker gene, the mean value and the percentage of cells
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
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value determines
        the number of groups into which the groupby observation should be subdivided.
    color_map : `str`, optional (default: `Reds`)
        String denoting matplotlib color map.
    figsize : (float, float), optional (default: None)
        Figure size (width, height. If not set, the figure width is set based on the
        number of  `var_names` and the height is set to 10.
    var_group_positions :  list of `tuples`.
        Use this parameter to highlight groups of `var_names`. This will draw a 'bracket'
        on top of the plot between the given start and end positions. If the
        parameter `var_group_labels` is set, the corresponding labels is added on
        top of the bracket. E.g. var_group_positions = [(4,10)] will add a bracket
        between the fourth var_name and the tenth var_name. By giving more
        positions, more brackets are drawn.
    var_group_labels : list of `str`
        Labels for each of the var_group_positions that want to be highlighted.
    var_group_rotation : `float` (default: `None`)
        Label rotation degrees. By default, labels larger than 4 characters are rotated 90 degrees
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `matplotlib.pyplot.scatter`.

    Returns
    -------
    A list of `matplotlib.Axes` where the first ax is the groupby categories colorcode, the
    second axis is the heatmap and the third axis is the colorbar.
    """
    if use_raw is None and adata.raw is not None: use_raw = True
    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories)

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

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.5, 10]
    else:
        height_ratios = [0, 10.5]

    # define a layout of 2 rows x 4 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
    # second row is for main content. This second row
    # is divided into 4 axes:
    #   first ax is for the main figure
    #   second ax is for the color bar legend
    #   third ax is for an spacer that avoids the ticks
    #    from the color bar to be hidden beneath the size lengend axis
    #   fourth ax is to plot the size legend
    from matplotlib import gridspec
    fig = pl.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=2, ncols=4, left=0.05, right=0.48, wspace=0.05, hspace=0.04,
                            width_ratios=[heatmap_width, colorbar_width, colorbar_width_spacer, size_legend_width],
                            height_ratios=height_ratios)
    if len(categories) < 4:
        # when few categories are shown, the colorbar and size legend
        # need to be larger than the main plot, otherwise they would look
        # compressed. For this, the dotplot ax is split into two:
        axs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[1, 0],
                                                height_ratios=[len(categories) * 0.3, 1])
        dot_ax = fig.add_subplot(axs2[0])
    else:
        dot_ax = fig.add_subplot(axs[1, 0])

    color_legend = fig.add_subplot(axs[1, 1])

    # to keep the size_legen of about the same height, irrespective
    # of the number of categories, the fourth ax is subdivided into two parts
    size_legend_height = min(1.3, height)
    # wspace is proportional to the width but a constant value is
    # needed such that the spacing is the same for thinner or wider images.
    wspace = 10.5 / width
    axs3 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[1, 3], wspace=wspace,
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
    cmap = pl.get_cmap(color_map)

    size = (frac * 10) ** 2

    import matplotlib.colors
    normalize = matplotlib.colors.Normalize(vmin=min(mean_flat), vmax=max(mean_flat))
    colors = [cmap(normalize(value)) for value in mean_flat]

    dot_ax.scatter(x, y, color=colors, s=size, cmap=cmap, norm=None, edgecolor='none', **kwds)
    y_ticks = range(mean_obs.shape[0])
    dot_ax.set_yticks(y_ticks)
    dot_ax.set_yticklabels([mean_obs.index[idx] for idx in y_ticks])

    x_ticks = range(mean_obs.shape[1])
    dot_ax.set_xticks(x_ticks)
    dot_ax.set_xticklabels([mean_obs.columns[idx] for idx in x_ticks], rotation=90)
    dot_ax.tick_params(axis='both', labelsize='small')
    dot_ax.grid(False)
    dot_ax.set_xlim(-0.5, len(var_names) + 0.5)
    dot_ax.set_ylabel(groupby)

    # to be consistent with the heatmap plot, is better to
    # invert the order of the y-axis, such that the first group is on
    # top
    ymin, ymax = dot_ax.get_ylim()
    dot_ax.set_ylim(ymax+0.5, ymin - 0.5)

    dot_ax.set_xlim(-1, len(var_names) + 0.5)

    # plot group legends on top of dot_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)
        _plot_gene_groups_brackets(gene_groups_ax, group_positions=var_group_positions,
                                   group_labels=var_group_labels,
                                   rotation=var_group_rotation)

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

    utils.savefig_or_show('dotplot', show=show, save=save)
    return axs


@doc_params(show_save_ax=doc_show_save_ax)
def matrixplot(adata, var_names, groupby=None, use_raw=None, log=False, num_categories=7,
               figsize=None, var_group_positions=None, var_group_labels=None,
               var_group_rotation=None, show=None, save=None, **kwds):
    """\
    Creates a heatmap of the mean expression values per cluster of each var_names
    If groupby is not given, the matrixplot assumes that all data belongs to a single
    category.

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
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value determines
        the number of groups into which the groupby observation should be subdivided.
    figsize : (float, float), optional (default: None)
        Figure size (width, height. If not set, the figure width is set based on the
        number of  `var_names` and the height is set to 10.
    var_group_positions :  list of `tuples`.
        Use this parameter to highlight groups of `var_names`. This will draw a 'bracket'
        on top of the plot between the given start and end positions. If the
        parameter `var_group_labels` is set, the corresponding labels is added on
        top of the bracket. E.g. var_group_positions = [(4,10)] will add a bracket
        between the fourth var_name and the tenth var_name. By giving more
        positions, more brackets are drawn.
    var_group_labels : list of `str`
        Labels for each of the var_group_positions that want to be highlighted.
    var_group_rotation : `float` (default: `None`)
        Label rotation degrees. By default, labels larger than 4 characters are rotated 90 degrees
    {show_save_ax}
    **kwds : keyword arguments
        Are passed to `matplotlib.pyplot.pcolor`.

    Returns
    -------
    A list of `matplotlib.Axes` where the first ax is the groupby categories colorcode, the
    second axis is the heatmap and the third axis is the colorbar.
    """
    if use_raw is None and adata.raw is not None: use_raw = True
    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories)

    mean_obs = obs_tidy.groupby(level=0).mean()
    if figsize is None:
        height = len(categories) * 0.2 + 1  # +1 for labels
        heatmap_width = len(var_names) * 0.6
        width = heatmap_width + 1.6 + 1  # +1.6 to account for the colorbar and  + 1 to account for labels
    else:
        width, height = figsize
        heatmap_width = width * 0.75

    # colorbar ax width should not change with differences in the width of the image
    colorbar_width = 0.4

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.5, 10]
        height += 0.5
    else:
        height_ratios = [0, 10.5]

    # define a layout of 2 rows x 2 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
    # second row is for main content. This second row
    # is divided into two axes:
    #   first ax is for the main matrix figure
    #   second ax is for the color bar legend
    from matplotlib import gridspec
    fig = pl.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=2, ncols=2, left=0.05, right=0.48, wspace=0.05, hspace=0.04,
                            width_ratios=[heatmap_width, colorbar_width],
                            height_ratios=height_ratios)
    matrix_ax = fig.add_subplot(axs[1, 0])

    color_legend = fig.add_subplot(axs[1, 1])

    pc = matrix_ax.pcolor(mean_obs, edgecolor='gray', **kwds)
    y_ticks = np.arange(mean_obs.shape[0]) + 0.5
    matrix_ax.set_yticks(y_ticks)
    matrix_ax.set_yticklabels([mean_obs.index[idx] for idx in range(mean_obs.shape[0])])
    # invert y axis to show categories ordered from top to bottom
    matrix_ax.set_ylim(mean_obs.shape[0], 0)

    x_ticks = np.arange(mean_obs.shape[1]) + 0.5
    matrix_ax.set_xticks(x_ticks)
    matrix_ax.set_xticklabels([mean_obs.columns[idx] for idx in range(mean_obs.shape[1])], rotation=90)
    matrix_ax.tick_params(axis='both', labelsize='small')
    matrix_ax.grid(False)
    matrix_ax.set_xlim(-0.5, len(var_names) + 0.5)
    matrix_ax.set_ylabel(groupby)
    matrix_ax.set_xlim(0, mean_obs.shape[1])
    # plot group legends on top of matrix_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=matrix_ax)
        _plot_gene_groups_brackets(gene_groups_ax, group_positions=var_group_positions,
                                   group_labels=var_group_labels, rotation=var_group_rotation,
                                   left_adjustment=0.2, right_adjustment=0.8)

    # plot colorbar
    pl.colorbar(pc, cax=color_legend)

    utils.savefig_or_show('matrixplot', show=show, save=save)
    return axs


def _prepare_dataframe(adata, var_names, groupby=None, use_raw=None, log=False, num_categories=7):
    """
    Given the anndata object, prepares a data frame in which the row index are the categories
    defined by group by and the columns correspond to var_names.

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
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present.
    num_categories : `int`, optional (default: `7`)
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.

    Returns
    -------
    Tuple of `pandas.DataFrame` and list of categories.
    """
    from scipy.sparse import issparse
    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None: use_raw = True
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
            # by subdividing into  `num_categories` categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]

    obs_tidy.set_index(categorical, groupby, inplace=True)
    categories = obs_tidy.index.categories

    return categories, obs_tidy


def _plot_gene_groups_brackets(gene_groups_ax, group_positions, group_labels,
                               left_adjustment=-0.3, right_adjustment=0.3, rotation=None):
    """
    Draws brackets that represent groups of genes on the give axis.
    For best results, this axis is located on top of an image whose
    x axis contains gene names.

    The gene_groups_ax should share the x axis with the main ax.

    Eg: gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)

    This function is used by dotplot, heatmap etc.

    Parameters
    ----------
    gene_groups_ax : matplotlib axis
        In this axis the gene marks are drawn
    group_positions : list of `tuples`
        Each item in the list, should contain the start and end position that the
        bracket should cover.
        Eg. [(0, 4), (5, 8)] means that there are two brackets, one for the var_names (eg genes)
        in positions 0-4 and other for positions 5-8
    group_labels :  list
        List of group labels
    left_adjustment : `float`
        adjustment to plot the bracket start slightly before or after the first gene position.
        If the value is negative the start is moved before.
    right_adjustment : `float`
        adjustment to plot the bracket end slightly before or after the last gene position
        If the value is negative the start is moved before.
    rotation : `float` (default None)
        rotation degrees for the labels. If not given, small labels (<4 characters) are not
        rotated, otherwise, they are rotated 90 degrees

    Returns
    -------
    None
    """
    import matplotlib.patches as patches
    from matplotlib.path import Path

    # get the 'brackets' coordinates as lists of start and end positions

    left = [x[0] + left_adjustment for x in group_positions]
    right = [x[1] + right_adjustment for x in group_positions]

    # verts and codes are used by PathPatch to make the brackets
    verts = []
    codes = []

    # rotate labels if any of them is longer than 4 characters
    if rotation is None and group_labels is not None and len(group_labels) > 0:
        if max([len(x) for x in group_labels]) > 4:
            rotation = 90
        else:
            rotation = 0
    for idx in range(len(left)):
        verts.append((left[idx], 0))  # lower-left
        verts.append((left[idx], 0.6))  # upper-left
        verts.append((right[idx], 0.6))  # upper-right
        verts.append((right[idx], 0))  # lower-right

        codes.append(Path.MOVETO)
        codes.append(Path.LINETO)
        codes.append(Path.LINETO)
        codes.append(Path.LINETO)

        try:
            group_x_center = left[idx] + float(right[idx] - left[idx]) / 2
            gene_groups_ax.text(group_x_center, 1.1, group_labels[idx], ha='center',
                                va='bottom', rotation=rotation)
        except:
            pass

    path = Path(verts, codes)

    patch = patches.PathPatch(path, facecolor='none', lw=1.5)

    gene_groups_ax.add_patch(patch)
    gene_groups_ax.spines['right'].set_visible(False)
    gene_groups_ax.spines['top'].set_visible(False)
    gene_groups_ax.spines['left'].set_visible(False)
    gene_groups_ax.spines['bottom'].set_visible(False)
    gene_groups_ax.grid(False)

    # remove y ticks
    gene_groups_ax.tick_params(axis='y', left=False, labelleft=False)
    # remove x ticks and labels
    gene_groups_ax.tick_params(axis='x', bottom=False, labelbottom=False, labeltop=False)
