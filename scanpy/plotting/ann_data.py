# Authors: F. Alex Wolf <http://falexwolf.de>
#          P. Angerer
"""Plotting functions for AnnData.
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import is_color_like
import seaborn as sns

from .. import settings as sett
from . import utils
from .utils import scatter_base, scatter_group


def scatter(
        adata,
        x=None,
        y=None,
        color='grey',
        basis=None,
        groups=None,
        components=None,
        projection='2d',
        legend_loc='right margin',
        legend_fontsize=None,
        color_map=None,
        palette=None,
        right_margin=None,
        size=None,
        title=None,
        show=None,
        save=None,
        ax=None):
    """Scatter plot.

    Color with sample annotation (`color in adata.smp_keys()`) or gene
    expression (`color in adata.var_names`).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    x : str or None
        x coordinate.
    y : str or None
        y coordinate.
    color : string or list of strings, optional (default: None)
        Keys for sample/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
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
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    A list of matplotlib.Axis objects.
    """
    if components is None: components = '1,2' if '2d' in projection else '1,2,3'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1
    title = None if title is None else title.split(',') if isinstance(title, str) else title
    color_keys = ['grey'] if color is None else color.split(',') if isinstance(color, str) else color
    groups = None if groups is None else groups.split(',') if isinstance(groups, str) else groups
    highlights = adata.add['highlights'] if 'highlights' in adata.add else []
    if basis is not None:
        try:
            Y = adata.smp['X_' + basis][:, components]
        except KeyError:
            raise KeyError('compute coordinates using visualization tool {} first'
                           .format(basis))
    elif x is not None and y is not None:
        x_arr = adata.get_smp_array(x)
        y_arr = adata.get_smp_array(y)
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
        palettes = [palette for i in range(len(color_keys))]
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
    color_ids = [None if not is_color_like(color_key)
                 else color_key for color_key in color_keys]
    categoricals = []
    colorbars = []
    for icolor_key, color_key in enumerate(color_keys):
        if color_ids[icolor_key] is not None:
            c = color_ids[icolor_key]
            continuous = True
            categorical = False
            colorbars.append(False)
        else:
            c = 'white' if projection == '2d' else 'white'
            categorical = False
            continuous = False
            # test whether we have categorial or continuous annotation
            if color_key in adata.smp_keys():
                if adata.smp[color_key].dtype.char in ['S', 'U']:
                    categorical = True
                else:
                    continuous = True
                    c = adata.smp[color_key]
                # sett.m(0, '... coloring according to', color_key)
            # coloring according to gene expression
            elif color_key in set(adata.var_names):
                c = adata[:, color_key].X
                continuous = True
                # sett.m(0, '... coloring according to expression of gene', color_key)
            else:
                raise ValueError('"' + color_key + '" is invalid!'
                                 + ' specify valid sample annotation, one of '
                                 + str(adata.smp_keys()) + ' or a gene name '
                                 + str(adata.var_names))
            colorbars.append(True if continuous else False)
        if categorical: categoricals.append(icolor_key)
        color_ids[icolor_key] = c

    if right_margin is None:
        if legend_loc == 'right margin':
            right_margin = 0.5
        # else:
        #     right_margin = 0.1
    if title is None and color_keys[0] is not None:
        title = [color_key.replace('_', ' ') if not is_color_like(color_key) else '' for color_key in color_keys]

    axs = scatter_base(Y,
                       title=title,
                       component_name=component_name,
                       axis_labels=axis_labels,
                       component_indexnames=components + 1,
                       projection=projection,
                       colors=color_ids,
                       highlights=highlights,
                       colorbars=colorbars,
                       right_margin=right_margin,
                       sizes=[size for c in color_keys],
                       color_map=color_map,
                       show_ticks=show_ticks,
                       ax=ax)

    def add_centroid(centroids, name, Y, mask):
        masked_values = Y[mask]
        if masked_values.shape[0] == 0: return
        median = np.median(masked_values, axis=0)
        i = np.argmin(np.sum(np.abs(masked_values - median), axis=1))
        centroids[name] = masked_values[i]

    for i, icolor_key in enumerate(categoricals):
        palette = palettes[i]
        color_key = color_keys[icolor_key]
        if (not color_key + '_colors' in adata.add or not palette_was_none
            or len(adata.add[color_key + '_names']) != len(adata.add[color_key + '_colors'])):
            utils.add_colors_for_categorical_sample_annotation(adata, color_key, palette)
        # actually plot the groups
        mask_remaining = np.ones(Y.shape[0], dtype=bool)
        centroids = {}
        if groups is None:
            for iname, name in enumerate(adata.add[color_key + '_names']):
                if name not in sett._ignore_categories:
                    mask = scatter_group(axs[icolor_key], color_key, iname,
                                         adata, Y, projection, size=size)
                    mask_remaining[mask] = False
                    if legend_loc == 'on data': add_centroid(centroids, name, Y, mask)
        else:
            for name in names:
                if name not in set(adata.add[color_key + '_names']):
                    raise ValueError('"' + name + '" is invalid!'
                                     + ' specify valid name, one of '
                                     + str(adata.add[color_key + '_names']))
                else:
                    iname = np.flatnonzero(adata.add[color_key + '_names'] == name)[0]
                    mask = scatter_group(axs[icolor_key], color_key, iname,
                                         adata, Y, projection, size=size)
                    if legend_loc == 'on data': add_centroid(centroids, name, Y, mask)
                    mask_remaining[mask] = False
        if mask_remaining.sum() > 0:
            data = [Y[mask_remaining, 0], Y[mask_remaining, 1]]
            if projection == '3d': data.append(Y[mask_remaining, 2])
            axs[icolor_key].scatter(*data, marker='.', c='grey', s=size,
                                    edgecolors='none', zorder=-1)
        legend = None
        if legend_loc == 'on data':
            for name, pos in centroids.items():
                axs[icolor_key].text(pos[0], pos[1], name,
                                     verticalalignment='center',
                                     horizontalalignment='center',
                                     fontsize=legend_fontsize)
        elif legend_loc == 'right margin':
            legend = axs[icolor_key].legend(frameon=False, loc='center left',
                                            bbox_to_anchor=(1, 0.5),
                                            ncol=(1 if len(adata.add[color_key + '_names']) <= 14
                                                  else 2 if len(adata.add[color_key + '_names']) <= 30 else 3),
                                            fontsize=legend_fontsize)
        elif legend_loc != 'none':
            legend = axs[icolor_key].legend(frameon=False, loc=legend_loc,
                                            fontsize=legend_fontsize)
        if legend is not None:
            for handle in legend.legendHandles: handle.set_sizes([300.0])
    utils.savefig_or_show('scatter' if basis is None else basis, show=show, save=save)
    return axs


def ranking(adata, attr, keys, labels=None, color='black', n_points=30,
            log=False):
    """Plot rankings.

    See, for example, how this is used in pl.pca_ranking.

    Parameters
    ----------
    adata : AnnData
        The data.
    attr : {'var', 'add', 'smp'}
        The attribute of AnnData that contains the score.
    keys : str or list of str
        The scores to look up an array from the attribute of adata.

    Returns
    -------
    Returns matplotlib gridspec with access to the axes.
    """
    scores = getattr(adata, attr)[keys]
    n_panels = len(keys) if isinstance(keys, list) else 1
    if n_panels == 1: scores, keys = scores[:, None], [keys]
    if log: scores = np.log(scores)
    if labels is None:
        labels = adata.var_names if attr == 'var' else np.arange(scores.shape[0]).astype(str)
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


def violin(adata, keys, group_by=None, jitter=True, size=1, scale='width',
           multi_panel=False, show=None, save=None, ax=None):
    """Violin plot.

    Wraps seaborn.violinplot for AnnData.

    Parameters
    ----------
    keys : str
        Keys for accessing fields of adata.smp.
    group_by : str
        Key that denotes grouping (categorical annotation) to index adata.smp.
    multi_panel : bool
        Show fields in multiple panels. Returns a seaborn FacetGrid in that case.
    jitter : float or bool (default: True)
        See sns.stripplot.
    scale : str (default: 'width')
        See sns.violinplot.
    show : bool, optional (default: None)
         Show the plot.
    save : bool or str, optional (default: None)
         If True or a str, save the figure. A string is appended to the
         default filename.
    ax : matplotlib.Axes
         A matplotlib axes object.

    Returns
    -------
    A matplotlib.Axes object.
    """
    if group_by is not None and isinstance(keys, list):
        raise ValueError('Pass a single key as string if using `group_by`.')
    if not isinstance(keys, list): keys = [keys]
    smp_keys = False
    for key in keys:
        if key in adata.smp_keys():
            smp_keys = True
        if smp_keys and key not in set(adata.smp_keys()):
            raise ValueError('Either use sample keys or variable names, but do not mix. '
                             'Did not find {} in adata.smp_keys().'.format(key))
    if smp_keys:
        smp_df = adata.smp.to_df()
    else:
        if group_by is None:
            smp_df = pd.DataFrame()
        else:
            smp_df = adata.smp.to_df()
        for key in keys:
            X_col = adata[:, key].X
            if issparse(X_col): X_col = X_col.toarray().flatten()
            smp_df[key] = X_col
    if group_by is None:
        smp_tidy = pd.melt(smp_df, value_vars=keys)
        x = 'variable'
        y = 'value'
        order = None
    else:
        smp_tidy = smp_df
        x = group_by
        y = keys[0]
        if not group_by + '_names' in adata.add:
            from .. import utils as sc_utils
            sc_utils.check_adata(adata)
        order = adata.add[group_by + '_names']
    if multi_panel:
        sns.set_style('whitegrid')
        g = sns.FacetGrid(smp_tidy, col=x, sharey=False)
        g = g.map(sns.violinplot, y, inner=None, orient='vertical', scale=scale)
        g = g.map(sns.stripplot, y, orient='vertical', jitter=jitter, size=size,
                     color='black').set_titles(
                         col_template='{col_name}').set_xlabels('')
        ax = g
    else:
        ax = sns.violinplot(x=x, y=y, data=smp_tidy, inner=None, order=order,
                            orient='vertical', scale=scale, ax=ax)
        ax = sns.stripplot(x=x, y=y, data=smp_tidy, order=order,
                           jitter=jitter, color='black', size=size, ax=ax)
        ax.set_xlabel('' if group_by is None else group_by.replace('_', ' '))
    utils.savefig_or_show('violin', show=show, save=save)
    return ax
