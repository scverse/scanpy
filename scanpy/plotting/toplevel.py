# Authors: F. Alex Wolf <http://falexwolf.de>
#          P. Angerer
"""Generic plotting functions for AnnData.
"""

import os
import numpy as np
from scipy.sparse import issparse

from ..compat.matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.figure import SubplotParams as sppars
from matplotlib.colors import is_color_like
from .. import settings as sett
from .. import logging as logg
from .. import utils as sc_utils
from .. import readwrite
from . import utils
from .utils import scatter_base, scatter_group


# -------------------------------------------------------------------------------
# Generic Helper Functions
# -------------------------------------------------------------------------------


def savefig(writekey, dpi=None, ext=None):
    """Save current figure to file.

    The filename is generated as follows:
    ```
    if sett.run_name != '': writekey = sett.run_name + '_' + writekey
    filename = sett.figdir + writekey + sett.plotsuffix + '.' + sett.file_format_figs
    ```
    """
    if dpi is None:
        if rcParams['savefig.dpi'] < 300:
            dpi = 300
            if sett._low_resolution_warning:
                logg.m('... you are using a very low resolution for saving figures, adjusting to dpi=300')
                sett._low_resolution_warning = False
        else:
            dpi = rcParams['savefig.dpi']
    if not os.path.exists(sett.figdir): os.makedirs(sett.figdir)
    if sett.run_name != '': writekey = sett.run_name + '_' + writekey
    if sett.figdir[-1] != '/': sett.figdir += '/'
    if ext is None: ext = sett.file_format_figs
    filename = sett.figdir + writekey + sett.plotsuffix + '.' + ext
    logg.info('... saving figure to file', filename)
    pl.savefig(filename, dpi=dpi)


def savefig_or_show(writekey, show=None, dpi=None, ext=None, save=None):
    if isinstance(save, str):
        writekey += save
        save = True
    save = sett.savefigs if save is None else save
    show = (sett.autoshow and not save) if show is None else show
    if save: savefig(writekey, dpi=dpi, ext=ext)
    if show: pl.show()
    if save: pl.close()  # clear figure


# -------------------------------------------------------------------------------
# Generic plotting functions
# -------------------------------------------------------------------------------


def matrix(matrix, xlabels=None, ylabels=None, cshrink=0.5,
           cmap='Greys', show=None, save=None, ax=None):
    """Plot a matrix."""
    if ax is None: ax = pl.gca()
    ax.imshow(matrix, cmap=cmap)
    if xlabels is not None:
        ax.set_xticks(range(len(xlabels)), xlabels, rotation='vertical')
    if ylabels is not None:
        ax.set_yticks(range(len(ylabels)), ylabels)
    ax.colorbar(shrink=cshrink)
    savefig_or_show('matrix', show=show, save=save)


def scatter(adata,
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
            pal=None,
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
    else:
        x_arr = adata.get_smp_array(x)
        y_arr = adata.get_smp_array(y)
        Y = np.c_[x_arr[:, None], y_arr[:, None]]

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
                       color_map='viridis' if color_map is None else color_map,
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
            for iname, name in enumerate(adata.add[color_key + '_groups']):
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
    savefig_or_show('scatter' if basis is None else basis, show=show, save=save)
    return axs


def ranking(adata, attr, keys, labels=None, color='black', n_points=30, log=False):
    """Plot rankings.

    See, for example, how this is used in pl.pca_ranking.

    Parameters
    ----------
    adata : AnnData
        The data.
    attr : {'var', 'add', 'smp'}
        An attribute of AnnData.
    keys : str or list of str
        Used to look up an array from the attribute of adata.

    Returns
    -------
    Returns matplotlib gridspec with access to the axes.
    """
    scores = getattr(adata, attr)[keys]
    n_panels = len(keys) if isinstance(keys, list) else 1
    if n_panels == 1: scores, keys = scores[:, None], [keys]
    if log: scores = np.log10(scores)
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


def ranking_deprecated(adata, toolkey, n_genes=20):
    """Plot ranking of genes

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_genes : int
        Number of genes.
    """

    # one panel for each ranking
    scoreskey = adata.add[toolkey + '_scoreskey']
    n_panels = len(adata.add[toolkey + '_rankings_names'])

    def get_scores(irank):
        allscores = adata.add[toolkey + '_' + scoreskey][irank]
        scores = allscores[adata.add[toolkey + '_rankings_geneidcs'][irank, :n_genes]]
        scores = np.abs(scores)
        return scores

    # number of panels
    if n_panels <= 5:
        n_panels_y = 1
        n_panels_x = n_panels
    else:
        n_panels_y = 2
        n_panels_x = int(n_panels/2+0.5)

    fig = pl.figure(figsize=(n_panels_x * 4, n_panels_y * 4))

    from matplotlib import gridspec
    left = 0.2/n_panels_x
    bottom = 0.13/n_panels_y
    gs = gridspec.GridSpec(nrows=n_panels_y,
                           ncols=n_panels_x,
                           left=left,
                           right=1-(n_panels_x-1)*left-0.01/n_panels_x,
                           bottom=bottom,
                           top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
                           wspace=0.17)

    count = 1
    for irank in range(len(adata.add[toolkey + '_rankings_names'])):
        pl.subplot(gs[count-1])
        scores = get_scores(irank)
        ymin = np.min(scores)
        ymax = np.max(scores)
        ymax += 0.3*(ymax-ymin)
        for ig, g in enumerate(adata.add[toolkey + '_rankings_geneidcs'][irank, :n_genes]):
            marker = (r'\leftarrow' if adata.add[toolkey + '_zscores'][irank, g] < 0
                                    else r'\rightarrow')
            pl.text(ig, scores[ig],
                    adata.var_names[g],
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center',
                    fontsize=8)
        title = adata.add[toolkey + '_rankings_names'][irank]
        pl.title(title)
        if n_panels <= 5 or count > n_panels_x:
            pl.xlabel('ranking')
        if count == 1 or count == n_panels_x+1:
            pl.ylabel(scoreskey)
        pl.ylim([ymin, ymax])
        pl.xlim(-0.9, ig+1-0.1)
        count += 1


def timeseries(X, **kwargs):
    """Plot X. See timeseries_subplot."""
    pl.figure(figsize=(2*rcParams['figure.figsize'][0], rcParams['figure.figsize'][1]),
              subplotpars=sppars(left=0.12, right=0.98, bottom=0.13))
    timeseries_subplot(X, **kwargs)


def timeseries_subplot(X,
                       c=None,
                       varnames=(),
                       highlightsX=(),
                       xlabel='',
                       ylabel='gene expression',
                       yticks=None,
                       xlim=None,
                       legend=True,
                       palette=None,
                       color_map='viridis'):
    """Plot X.

    Call this with:
    X with one column, c categorical
    X with one column, c continuous
    X with n columns, c is of length n
    """

    if c is not None:
        use_color_map = isinstance(c[0], float) or isinstance(c[0], np.float32)
    palette = utils.default_palette(palette)
    x_range = np.arange(X.shape[0])
    if X.shape[1] > 1:
        colors = palette[:X.shape[1]].by_key()['color']
        subsets = [(x_range, X[:, i]) for i in range(X.shape[1])]
    elif use_color_map:
        colors = [c]
        subsets = [(x_range, X[:, 0])]
    else:
        levels, _ = np.unique(c, return_inverse=True)
        colors = np.array(palette[:len(levels)].by_key()['color'])
        subsets = [(x_range[c == l], X[c == l, :]) for l in levels]

    for i, (x, y) in enumerate(subsets):
        pl.scatter(
            x, y,
            marker='.',
            edgecolor='face',
            s=rcParams['lines.markersize'],
            c=colors[i],
            label=varnames[i] if len(varnames) > 0 else '',
            color_map=color_map,
        )
    ylim = pl.ylim()
    for ih, h in enumerate(highlightsX):
        pl.plot([h, h], [ylim[0], ylim[1]], '--', color='black')
    pl.ylim(ylim)
    if xlim is not None:
        pl.xlim(xlim)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if yticks is not None:
        pl.yticks(yticks)
    if len(varnames) > 0 and legend == True:
        pl.legend(frameon=False)


def timeseries_as_heatmap(X, varnames=None, highlightsX=None, color_map='viridis'):
    """Plot timeseries as heatmap.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    varnames : array_like
        Array of strings naming variables stored in columns of X.
    """
    if highlightsX is None:
        highlightsX = []
    if varnames is None:
        varnames = []
    if len(varnames) == 0:
        varnames = np.arange(X.shape[1])
    if varnames.ndim == 2:
        varnames = varnames[:, 0]

    # transpose X
    X = X.T
    minX = np.min(X)

    # insert space into X
    if False:
        # generate new array with highlightsX
        space = 10  # integer
        Xnew = np.zeros((X.shape[0], X.shape[1] + space*len(highlightsX)))
        hold = 0
        _hold = 0
        space_sum = 0
        for ih, h in enumerate(highlightsX):
            _h = h + space_sum
            Xnew[:, _hold:_h] = X[:, hold:h]
            Xnew[:, _h:_h+space] = minX * np.ones((X.shape[0], space))
            # update variables
            space_sum += space
            _hold = _h + space
            hold = h
        Xnew[:, _hold:] = X[:, hold:]

    fig = pl.figure(figsize=(1.5*4, 2*4))
    im = pl.imshow(np.array(X, dtype=np.float_), aspect='auto',
                   interpolation='nearest', color_map=color_map)
    pl.colorbar(shrink=0.5)
    pl.yticks(range(X.shape[0]), varnames)
    for ih, h in enumerate(highlightsX):
        pl.plot([h, h], [0, X.shape[0]], '--', color='black')
    pl.xlim([0, X.shape[1]-1])
    pl.ylim([0, X.shape[0]-1])


def violin(adata, keys, group_by=None, jitter=True, size=1, scale='width',
           multi_panel=False, show=None, save=None, ax=None):
    """Violin plot.

    Wraps seaborn.violinplot.

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
    import pandas as pd
    import seaborn as sns
    utils.init_plotting_params()  # reset fig_params, seaborn overwrites settings
    sett.set_dpi()  # reset resolution
    smp_keys = False
    for key in keys:
        if key in adata.smp_keys():
            smp_keys = True
        if smp_keys and key not in adata.smp_keys():
            raise ValueError('Either use sample keys or variable names, but do not mix.')
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
    else:
        smp_tidy = smp_df
        x = group_by
        y = keys[0]
    if multi_panel:
        sns.set_style('whitegrid')
        g = sns.FacetGrid(smp_tidy, col=x, sharey=False)
        g = g.map(sns.violinplot, y, inner=None, orient='vertical', scale=scale)
        g = g.map(sns.stripplot, y, orient='vertical', jitter=jitter, size=size,
                     color='black').set_titles(
                         col_template='{col_name}').set_xlabels('')
        ax = g
    else:
        ax = sns.violinplot(x=x, y=y, data=smp_tidy, inner=None,
                            orient='vertical', scale=scale, ax=ax)
        ax = sns.stripplot(x=x, y=y, data=smp_tidy,
                           jitter=jitter, color='black', size=size, ax=ax)
        ax.set_xlabel('' if group_by is None else group_by)
    savefig_or_show('violin', show=show, save=save)
    return ax
