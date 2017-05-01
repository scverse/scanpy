# Authors: F. Alex Wolf <http://falexwolf.de>
#          P. Angerer
"""Toplevel plotting functions for AnnData.
"""

import os
import numpy as np

from ..compat.matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib.figure import SubplotParams as sppars
from matplotlib.colors import is_color_like
from .. import settings as sett
from .. import utils as sc_utils
from .. import readwrite
from . import utils
from .utils import scatter_base, scatter_group

# -------------------------------------------------------------------------------
# Toplevel Helper Functions
# -------------------------------------------------------------------------------


def savefig(writekey):
    """Save figure to `sett.figdir + writekey + '.' + sett.file_format_figures`
    """
    filename = sett.figdir + writekey + '.' + sett.file_format_figures
    sett.m(0, 'saving figure to file', filename)
    pl.savefig(filename)


# -------------------------------------------------------------------------------
# Toplevel Plotting Functions
# -------------------------------------------------------------------------------


def violin(adata, smp, show=True):
    import pandas as pd
    import seaborn as sns
    smp_df = adata.smp.to_df()
    smp_tidy = pd.melt(smp_df, value_vars=smp)
    sns.set_style('whitegrid')
    g = sns.FacetGrid(smp_tidy, col='variable', sharey=False).map(
        sns.violinplot, 'value', inner='quartile', orient='vertical').set_titles(
            col_template='{col_name}').set_xlabels('')
    if show:
        pl.show()
    return g


def scatter(adata,
            x=None,
            y=None,
            color='grey',
            basis=None,
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
            show=True):
    """Scatter plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    x : str or None
        x coordinate.
    y : str or None
        y coordinate.
    color : str or list of strings, optional (default: 'grey')
        Sample/Cell annotation key for coloring (either a key for adata.smp or a
        var_name or a uniform matplotlib color). String annotation is plotted assuming categorical annotation,
        float and integer annotation is plotted assuming continuous
        annoation.
    basis : {'pca', 'tsne', 'diffmap'}
        String that denotes a plotting tool.
    names : str, optional (default: all names in color)
        Allows to restrict groups in sample annotation (color) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
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
    # write params to a config file
    params = locals()
    del params['adata']
    if os.path.exists('.scanpy/config_plotting.txt'):
        params = sc_utils.update_params(readwrite.read_params('.scanpy/config_plotting.txt', verbosity=1), params)
        if right_margin != params['right_margin']:
            right_margin = params['right_margin']
            sett.m(1, '... setting right_margin to saved value', right_margin)
    readwrite.write_params('.scanpy/config_plotting.txt', params)
    del params
    # compute components
    if comps is None:
        comps = '1,2' if '2d' in layout else '1,2,3'
    if isinstance(comps, str):
        comps = comps.split(',')
    comps = np.array(comps).astype(int) - 1
    titles = None if titles is None else titles.split(',') if isinstance(titles, str) else titles
    color_keys = [None] if color is None else color.split(',') if isinstance(color, str) else color
    names = None if names is None else names.split(',') if isinstance(names, str) else names
    # highlights
    highlights = adata.add['highlights'] if 'highlights' in adata.smp else []
    if basis is not None:
        try:
            Y = adata.smp['X_' + basis][:, comps]
        except KeyError:
            sett.mi('--> compute the basis using plotting tool', basis, 'first')
            raise
    else:
        x_arr = adata.get_smp_array(x)
        y_arr = adata.get_smp_array(y)
        Y = np.c_[x_arr[:, None], y_arr[:, None]]

    pal = utils.default_pal(pal)

    component_name = ('DC' if basis == 'diffmap'
                      else 'Spring' if basis == 'spring'
                      else 'tSNE' if basis == 'tsne'
                      else 'PC' if basis == 'pca'
                      else None)
    axis_labels = (x, y) if component_name is None else None

    # the actual color ids, e.g. 'grey' or '#109482'
    color_ids = [None if not is_color_like(color_key) else color_key for color_key in color_keys]
    categoricals = []
    colorbars = []
    sizes = []
    for icolor_key, color_key in enumerate(color_keys):
        if color_ids[icolor_key] is not None:
            c = color_ids[icolor_key]
            continuous = True
            categorical = False
            colorbars.append(False)
        else:
            c = 'grey' if layout == '2d' else 'white'
            categorical = False
            continuous = False
            if len(adata.smp_keys()) > 0:
                if color_key is None:
                    for color_key_ in adata.smp_keys():
                        if color_key_ not in color_keys:
                            color_key = color_key_
                            color_keys[icolor_key] = color_key
                            break
                # test whether we have categorial or continuous annotation
                if color_key in adata.smp_keys():
                    if adata.smp[color_key].dtype.char in ['S', 'U']:
                        categorical = True
                        if cont is True:
                            c = adata.smp[color_key]
                    else:
                        continuous = True
                        c = adata.smp[color_key]
                    # sett.m(0, '... coloring according to', color_key)
                # coloring according to gene expression
                elif color_key in set(adata.var_names):
                    c = adata.X[:, color_key]
                    continuous = True
                    # sett.m(0, '... coloring according to expression of gene', color_key)
                else:
                    raise ValueError('"' + color_key + '" is invalid!'
                                     + ' specify valid sample annotation, one of '
                                     + str(adata.smp_keys()) + ' or a gene name '
                                     + str(adata.var_names))
            if cont is not None:
                categorical = not cont
                continuous = cont
            colorbars.append(True if continuous else False)
        sizes.append(size if continuous else 0.66*size)
        if categorical:
            categoricals.append(icolor_key)
        color_ids[icolor_key] = c

    if right_margin is None and legendloc == 'right margin':
        right_margin = 0.3
    if titles is None and color_keys[0] is not None:
        titles = [color_key.replace('_', ' ') for color_key in color_keys]

    axs = scatter_base(Y,
                       titles=titles,
                       component_name=component_name,
                       axis_labels=axis_labels,
                       component_indexnames=comps + 1,
                       layout=layout,
                       colors=color_ids,
                       highlights=highlights,
                       colorbars=colorbars,
                       right_margin=right_margin,
                       sizes=sizes,
                       cmap='viridis' if cmap is None else cmap)

    for icolor_key in categoricals:
        color_key = color_keys[icolor_key]
        if (color_key != 'groups' and 'groups_names' in adata.add
            and len(np.setdiff1d(adata.add['groups_names'], adata.add[color_key + '_names']))
                < len(adata.add['groups_names'])):
            # if there is a correspondence between color_key and the 'groups' defined
            # in adata, that is, if color_key has corresponding categories with those
            # in adata.add['groups_names']
            adata.add[color_key + '_colors'] = pal[:len(adata.add['groups_names'])].by_key()['color']
        elif not color_key + '_colors' in adata.add:
            adata.add[color_key + '_colors'] = pal[:len(adata.add[color_key + '_names'])].by_key()['color']
        if len(adata.add[color_key + '_names']) > len(adata.add[color_key + '_colors']):
            sett.m(0, 'number of categories/names in', color_key, 'so large that color map "jet" is used')
            adata.add[color_key + '_colors'] = pl.cm.get_cmap(cmap)(
                pl.Normalize()(np.arange(len(adata.add[color_key + '_names']), dtype=int)))
        for iname, name in enumerate(adata.add[color_key + '_names']):
            if names is None or (names is not None and name in names):
                scatter_group(axs[icolor_key], color_key, iname, adata, Y, layout, size=size)
        if legendloc == 'right margin':
            legend = axs[icolor_key].legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
        elif legendloc != 'none':
            axs[icolor_key].legend(frameon=False, loc=legendloc)
    return color_keys


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
                       pal=None,
                       cmap='viridis'):
    """Plot X.

    Call this with:
    X with one column, c categorical
    X with one column, c continuous
    X with n columns, c is of length n
    """

    if c is not None:
        use_cmap = isinstance(c[0], float) or isinstance(c[0], np.float32)
    pal = utils.default_pal(pal)
    x_range = np.arange(X.shape[0])
    if X.shape[1] > 1:
        colors = pal[:X.shape[1]].by_key()['color']
        subsets = [(x_range, X[:, i]) for i in range(X.shape[1])]
    elif use_cmap:
        colors = [c]
        subsets = [(x_range, X[:, 0])]
    else:
        levels, _ = np.unique(c, return_inverse=True)
        colors = np.array(pal[:len(levels)].by_key()['color'])
        subsets = [(x_range[c == l], X[c == l, :]) for l in levels]

    for i, (x, y) in enumerate(subsets):
        pl.scatter(
            x, y,
            marker='.',
            edgecolor='face',
            s=rcParams['lines.markersize'],
            c=colors[i],
            label=varnames[i] if len(varnames) > 0 else '',
            cmap=cmap,
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


def timeseries_as_heatmap(X, varnames=None, highlightsX=None, cmap='viridis'):
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
              interpolation='nearest', cmap=cmap)
    pl.colorbar(shrink=0.5)
    pl.yticks(range(X.shape[0]), varnames)
    for ih, h in enumerate(highlightsX):
        pl.plot([h, h], [0, X.shape[0]], '--', color='black')
    pl.xlim([0, X.shape[1]-1])
    pl.ylim([0, X.shape[0]-1])


def ranking(adata, toolkey, n_genes=20):
    """
    Plot ranking of genes

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

    # the limits for the y axis
    ymin = 1e100
    ymax = -1e100
    for irank in range(len(adata.add[toolkey + '_rankings_names'])):
        scores = get_scores(irank)
        ymin = np.min([ymin, np.min(scores)])
        ymax = np.max([ymax, np.max(scores)])
    ymax += 0.3*(ymax-ymin)

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
                           wspace=0)

    count = 1
    for irank in range(len(adata.add[toolkey + '_rankings_names'])):
        pl.subplot(gs[count-1])
        scores = get_scores(irank)
        for ig, g in enumerate(adata.add[toolkey + '_rankings_geneidcs'][irank, :n_genes]):
            marker = (r'\leftarrow' if adata.add[toolkey + '_zscores'][irank, g] < 0
                                    else r'\rightarrow')
            pl.text(ig, scores[ig],
                    r'$ ' + marker + '$ ' + adata.var_names[g],
                    color='red' if adata.add[toolkey + '_zscores'][irank, g] < 0 else 'green',
                    rotation='vertical', verticalalignment='bottom',
                    horizontalalignment='center',
                    fontsize=8)
        title = adata.add[toolkey + '_rankings_names'][irank]
        pl.title(title)
        if n_panels <= 5 or count > n_panels_x:
            pl.xlabel('ranking')
        if count == 1 or count == n_panels_x+1:
            pl.ylabel(scoreskey)
        else:
            pl.yticks([])
        pl.ylim([ymin, ymax])
        pl.xlim(-0.9, ig+1-0.1)
        count += 1
