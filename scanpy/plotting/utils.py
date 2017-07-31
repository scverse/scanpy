# Authors: F. Alex Wolf (http://falexwolf.de)
#          P. Angerer

import os
import numpy as np
import networkx as nx
from matplotlib import pyplot as pl
from matplotlib import rcParams, ticker
from matplotlib.colors import is_color_like
from matplotlib.figure import SubplotParams as sppars
from cycler import Cycler, cycler
from .. import logging as logg
from .. import settings as sett
from . import palettes


# -------------------------------------------------------------------------------
# Simple plotting functions
# -------------------------------------------------------------------------------


def matrix(matrix, xlabels=None, ylabels=None, colorbar_shrink=0.5,
           color_map=None, show=None, save=None, ax=None):
    """Plot a matrix."""
    if ax is None: ax = pl.gca()
    img = ax.imshow(matrix, cmap=color_map)
    if xlabels is not None:
        ax.set_xticks(range(len(xlabels)), xlabels, rotation='vertical')
    if ylabels is not None:
        ax.set_yticks(range(len(ylabels)), ylabels)
    pl.colorbar(img, shrink=colorbar_shrink, ax=ax)  # need a figure instance for colorbar
    savefig_or_show('matrix', show=show, save=save)


def timeseries(X, **kwargs):
    """Plot X. See timeseries_subplot."""
    pl.figure(figsize=(2*rcParams['figure.figsize'][0], rcParams['figure.figsize'][1]),
              subplotpars=sppars(left=0.12, right=0.98, bottom=0.13))
    timeseries_subplot(X, **kwargs)


def timeseries_subplot(X,
                       time=None,
                       color=None,
                       var_names=(),
                       highlightsX=(),
                       xlabel='',
                       ylabel='gene expression',
                       yticks=None,
                       xlim=None,
                       legend=True,
                       palette=None,
                       color_map='viridis'):
    """Plot X.

    Parameters
    ----------
    X : np.ndarray
        Call this with:
        X with one column, color categorical.
        X with one column, color continuous.
        X with n columns, color is of length n.
    """

    if color is not None:
        use_color_map = isinstance(color[0], float) or isinstance(color[0], np.float32)
    palette = default_palette(palette)
    x_range = np.arange(X.shape[0]) if time is None else time
    if X.ndim == 1: X = X[:, None]
    if X.shape[1] > 1:
        colors = palette[:X.shape[1]].by_key()['color']
        subsets = [(x_range, X[:, i]) for i in range(X.shape[1])]
    elif use_color_map:
        colors = [color]
        subsets = [(x_range, X[:, 0])]
    else:
        levels, _ = np.unique(color, return_inverse=True)
        colors = np.array(palette[:len(levels)].by_key()['color'])
        subsets = [(x_range[color == l], X[color == l, :]) for l in levels]

    for i, (x, y) in enumerate(subsets):
        pl.scatter(
            x, y,
            marker='.',
            edgecolor='face',
            s=rcParams['lines.markersize'],
            c=colors[i],
            label=var_names[i] if len(var_names) > 0 else '',
            cmap=color_map)
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
    if len(var_names) > 0 and legend:
        pl.legend(frameon=False)


def timeseries_as_heatmap(X, var_names=None, highlightsX=None, color_map=None):
    """Plot timeseries as heatmap.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    var_names : array_like
        Array of strings naming variables stored in columns of X.
    """
    if highlightsX is None:
        highlightsX = []
    if var_names is None:
        var_names = []
    if len(var_names) == 0:
        var_names = np.arange(X.shape[1])
    if var_names.ndim == 2:
        var_names = var_names[:, 0]

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
                   interpolation='nearest', cmap=color_map)
    pl.colorbar(shrink=0.5)
    pl.yticks(range(X.shape[0]), var_names)
    for ih, h in enumerate(highlightsX):
        pl.plot([h, h], [0, X.shape[0]], '--', color='black')
    pl.xlim([0, X.shape[1]-1])
    pl.ylim([0, X.shape[0]-1])


# -------------------------------------------------------------------------------
# Colors in additional to matplotlib's colors
# -------------------------------------------------------------------------------
    

additional_colors = {'gold2': '#eec900', 'firebrick3': '#cd2626', 'khaki2':
            '#eee685', 'slategray3': '#9fb6cd', 'palegreen3': '#7ccd7c',
            'tomato2': '#ee5c42', 'grey80': '#cccccc', 'grey90': '#e5e5e5',
            'wheat4': '#8b7e66', 'grey65': '#a6a6a6', 'grey10': '#1a1a1a',
            'grey20': '#333333', 'grey50': '#7f7f7f', 'grey30': '#4d4d4d',
            'grey40': '#666666', 'antiquewhite2': '#eedfcc', 'grey77':
            '#c4c4c4', 'snow4': '#8b8989', 'chartreuse3': '#66cd00', 'yellow4':
            '#8b8b00', 'darkolivegreen2': '#bcee68', 'olivedrab3': '#9acd32',
            'azure3': '#c1cdcd', 'violetred': '#d02090', 'mediumpurple3':
            '#8968cd', 'purple4': '#551a8b', 'seagreen4': '#2e8b57'}

    
# -------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------


def savefig(writekey, dpi=None, ext=None):
    """Save current figure to file.

    The filename is generated as follows:
    ```
    if sett.run_name != '': writekey = sett.run_name + '_' + writekey
    filename = sett.figdir + writekey + sett.plot_suffix + '.' + sett.file_format_figs
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
    filename = sett.figdir + writekey + sett.plot_suffix + '.' + ext
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


def default_palette(palette=None):
    if palette is None: return rcParams['axes.prop_cycle']
    elif not isinstance(palette, Cycler): return cycler(color=palette)
    else: return palette


def adjust_palette(palette, length):
    if len(palette.by_key()['color']) < length:
        if length <= 28:
            palette = palettes.default_26
        else:
            palette = palettes.default_64
        logg.m('... updating the color palette to provide enough colors')
        return cycler(color=palette)
    elif not isinstance(palette, Cycler):
        return cycler(color=palette)
    else:
        return palette


def add_colors_for_categorical_sample_annotation(adata, key, palette=None):
    if (key + '_colors' in adata.add
        and len(adata.add[key + '_names']) > len(adata.add[key + '_colors'])):
        logg.info('    number of defined colors does not match number of categories,'
                  ' using palette')
    else:
        logg.m('generating colors for {} using palette'.format(key), v=4)
    palette = default_palette(palette)
    palette_adjusted = adjust_palette(palette, length=len(adata.add[key + '_names']))
    adata.add[key + '_colors'] = palette_adjusted[:len(adata.add[key + '_names'])].by_key()['color']
    if len(adata.add[key + '_names']) > len(adata.add[key + '_colors']):
        raise ValueError('Cannot plot more than {} categories, which is not enough for {}.'
                         .format(len(adata.add[key + '_colors']), key))


def scatter_group(ax, name, imask, adata, Y, projection='2d', size=3):
    """Scatter of group using representation of data Y.
    """
    if name + '_masks' in adata.add:
        mask = adata.add[name + '_masks'][imask]
    else:
        if adata.add[name + '_names'][imask] in adata.smp[name]:
            mask = adata.add[name + '_names'][imask] == adata.smp[name]
        else:
            mask = str(imask) == adata.smp[name]
    color = adata.add[name + '_colors'][imask]
    if not isinstance(color[0], str):
        from matplotlib.colors import rgb2hex
        color = rgb2hex(adata.add[name + '_colors'][imask])
    if not is_color_like(color):
        raise ValueError('"{}" is not a valid matplotlib color.'.format(color))
    data = [Y[mask, 0], Y[mask, 1]]
    if projection == '3d': data.append(Y[mask, 2])
    ax.scatter(*data,
               marker='.',
               c=color,
               edgecolors='none',
               s=size,
               label=adata.add[name + '_names'][imask])
    return mask


def scatter_base(Y,
                 colors='blue',
                 highlights=[],
                 right_margin=None,
                 projection='2d',
                 title=None,
                 component_name='DC',
                 component_indexnames=[1, 2, 3],
                 axis_labels=None,
                 colorbars=[False],
                 sizes=[1],
                 color_map='viridis',
                 show_ticks=True,
                 ax=None):
    """Plot scatter plot of data.

    Parameters
    ----------
    Y : np.ndarray
        Data array.
    projection : {'2d', '3d'}

    Returns
    -------
    axs : matplotlib.axis or list of matplotlib.axis
        Depending on whether supplying a single array or a list of arrays,
        return a single axis or a list of axes.
    """
    if '3d' in projection: from mpl_toolkits.mplot3d import Axes3D
    if isinstance(highlights, dict):
        highlights_indices = sorted(highlights)
        highlights_labels = [highlights[i] for i in highlights_indices]
    else:
        highlights_indices = highlights
        highlights_labels = []
    # if we have a single array, transform it into a list with a single array
    avail_projections = {'2d', '3d'}
    if projection not in avail_projections:
        raise ValueError('choose projection from', avail_projections)
    if type(colors) == str: colors = [colors]
    if len(sizes) != len(colors):
        if len(sizes) == 1:
            sizes = [sizes[0] for i in range(len(colors))]
    # grid of axes for plotting and legends/colorbars
    if np.any(colorbars) and right_margin is None: right_margin = 0.25
    elif right_margin is None: right_margin = 0.10
    # make a list of right margins for each panel
    if not isinstance(right_margin, list):
        right_margin_list = [right_margin for i in range(len(colors))]
    else:
        right_margin_list = right_margin
    # make a figure with len(colors) panels in a row side by side
    top_offset = 1 - rcParams['figure.subplot.top']
    bottom_offset = 0.15 if show_ticks else 0.08
    left_offset = 1 if show_ticks else 0.3  # in units of base_height
    base_height = rcParams['figure.figsize'][1]
    height = base_height
    base_width = rcParams['figure.figsize'][0]
    if show_ticks: base_width *= 1.1
    draw_region_width = base_width - left_offset - top_offset - 0.5  # this is kept constant throughout

    right_margin_factor = sum([1 + right_margin for right_margin in right_margin_list])
    width_without_offsets = right_margin_factor * draw_region_width  # this is the total width that keeps draw_region_width

    right_offset = (len(colors) - 1) * left_offset
    figure_width = width_without_offsets + left_offset + right_offset
    draw_region_width_frac = draw_region_width / figure_width
    left_offset_frac = left_offset / figure_width
    right_offset_frac = 1 - (len(colors) - 1) * left_offset_frac

    if ax is None:
        fig = pl.figure(figsize=(figure_width, height),
                        subplotpars=sppars(left=0, right=1, bottom=bottom_offset))
    left_positions = [left_offset_frac, left_offset_frac + draw_region_width_frac]
    for i in range(1, len(colors)):
        right_margin = right_margin_list[i-1]
        left_positions.append(left_positions[-1] + right_margin * draw_region_width_frac)
        left_positions.append(left_positions[-1] + draw_region_width_frac)
    panel_pos = [[bottom_offset], [1-top_offset], left_positions]
    axs_passed = []
    if ax is not None: axs_passed = ax if isinstance(ax, list) else [ax]
    axs = []
    for icolor, color in enumerate(colors):
        left = panel_pos[2][2*icolor]
        bottom = panel_pos[0][0]
        width = draw_region_width / figure_width
        height = panel_pos[1][0] - bottom
        if projection == '2d':
            if axs_passed: ax = axs_passed[icolor]
            else: ax = pl.axes([left, bottom, width, height])
            data = Y[:, 0], Y[:, 1]
        elif projection == '3d':
            ax = pl.axes([left, bottom, width, height], projection='3d')
            data = Y[:, 0], Y[:, 1], Y[:, 2]
        if not isinstance(color, str) or color != 'white':
            sct = ax.scatter(*data,
                             marker='.',
                             c=color,
                             edgecolors='none',  # 'face',
                             s=sizes[icolor],
                             cmap=color_map)
        if colorbars[icolor]:
            width = 0.006 * draw_region_width
            left = panel_pos[2][2*icolor+1] + (1.2 if projection == '3d' else 0.2) * width
            rectangle = [left, bottom, width, height]
            ax_cb = fig.add_axes(rectangle)
            cb = pl.colorbar(sct, format=ticker.FuncFormatter(ticks_formatter),
                             cax=ax_cb)
        # set the title
        if title is not None: ax.set_title(title[icolor])
        # output highlighted data points
        for iihighlight, ihighlight in enumerate(highlights_indices):
            data = [Y[ihighlight, 0]], [Y[ihighlight, 1]]
            if '3d' in projection:
                data = [Y[ihighlight, 0]], [Y[ihighlight, 1]], [Y[ihighlight, 2]]
            ax.scatter(*data, c='black',
                       facecolors='black', edgecolors='black',
                       marker='x', s=10, zorder=20)
            highlight_text = (highlights_labels[iihighlight] if
                              len(highlights_labels) > 0
                              else str(ihighlight))
            # the following is a Python 2 compatibility hack
            ax.text(*([d[0] for d in data] + [highlight_text]),
                    zorder=20,
                    fontsize=10,
                    color='black')
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
            if '3d' in projection: ax.set_zticks([])
        axs.append(ax)
    # set default axis_labels
    if axis_labels is None:
        axis_labels = [[component_name + str(i) for i in idcs]
                       for idcs in
                       [component_indexnames for iax in range(len(axs))]]
    else:
        axis_labels = [[axis_labels[0], axis_labels[1]] for i in range(len(axs))]
    for iax, ax in enumerate(axs):
        ax.set_xlabel(axis_labels[iax][0])
        ax.set_ylabel(axis_labels[iax][1])
        if '3d' in projection:
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[iax][2], labelpad=-7)
    for ax in axs:
        # scale limits to match data
        ax.autoscale_view()
    return axs


def scatter_single(ax, Y, *args, **kwargs):
    """Plot scatter plot of data.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis to plot on.
    Y : np.array
        Data array, data to be plotted needs to be in the first two columns.
    """
    if 's' not in kwargs:
        kwargs['s'] = 2 if Y.shape[0] > 500 else 10
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'face'
    ax.scatter(Y[:, 0], Y[:, 1], **kwargs)
    ax.set_xticks([])
    ax.set_yticks([])


def arrows_transitions(ax, X, indices, weight=None):
    """
    Plot arrows of transitions in data matrix.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    X : np.array
        Data array, any representation wished (X, psi, phi, etc).
    indices : array_like
        Indices storing the transitions.
    """
    step = 1
    width = axis_to_data(ax, 0.001)
    if X.shape[0] > 300:
        step = 5
        width = axis_to_data(ax, 0.0005)
    if X.shape[0] > 500:
        step = 30
        width = axis_to_data(ax, 0.0001)
    head_width = 10*width
    for ix, x in enumerate(X):
        if ix % step == 0:
            X_step = X[indices[ix]] - x
            # don't plot arrow of length 0
            for itrans in range(X_step.shape[0]):
                alphai = 1
                widthi = width
                head_widthi = head_width
                if weight is not None:
                    alphai *= weight[ix, itrans]
                    widthi *= weight[ix, itrans]
                if np.any(X_step[itrans, :1]):
                    ax.arrow(x[0], x[1],
                             X_step[itrans, 0], X_step[itrans, 1],
                             length_includes_head=True,
                             width=widthi,
                             head_width=head_widthi,
                             alpha=alphai,
                             color='grey')


def ticks_formatter(x, pos):
    # pretty scientific notation
    if False:
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
    else:
        return ('%.3f' % (x)).rstrip('0').rstrip('.')


def pimp_axis(x_or_y_ax):
    """Remove trailing zeros.
    """
    x_or_y_ax.set_major_formatter(ticker.FuncFormatter(ticks_formatter))


def scale_to_zero_one(x):
    """Take some 1d data and scale it so that min matches 0 and max 1.
    """
    xscaled = x - np.min(x)
    xscaled /= np.max(xscaled)
    return xscaled


def hierarchy_pos(G, root, levels=None, width=1., height=1.):
    """Tree layout for networkx graph.

       See https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
       answer by burubum.

       If there is a cycle that is reachable from root, then this will see
       infinite recursion.

       Parameters
       ----------
       G: the graph
       root: the root node
       levels: a dictionary
               key: level number (starting from 0)
               value: number of nodes in this level
       width: horizontal space allocated for drawing
       height: vertical space allocated for drawing
    """
    TOTAL = "total"
    CURRENT = "current"

    def make_levels(levels, node=root, currentLevel=0, parent=None):
        """Compute the number of nodes for each level
        """
        if currentLevel not in levels:
            levels[currentLevel] = {TOTAL: 0, CURRENT: 0}
        levels[currentLevel][TOTAL] += 1
        neighbors = G.neighbors(node)
        if parent is not None:
            neighbors.remove(parent)
        for neighbor in neighbors:
            levels = make_levels(levels, neighbor, currentLevel + 1, node)
        return levels

    def make_pos(pos, node=root, currentLevel=0, parent=None, vert_loc=0):
        dx = 1/levels[currentLevel][TOTAL]
        left = dx/2
        pos[node] = ((left + dx*levels[currentLevel][CURRENT])*width,
                     vert_loc)
        levels[currentLevel][CURRENT] += 1
        neighbors = G.neighbors(node)
        if parent is not None:
            neighbors.remove(parent)
        for neighbor in neighbors:
            pos = make_pos(pos, neighbor, currentLevel + 1, node, vert_loc-vert_gap)
        return pos

    if levels is None:
        levels = make_levels({})
    else:
        levels = {l: {TOTAL: levels[l], CURRENT: 0} for l in levels}
    vert_gap = height / (max([l for l in levels])+1)
    return make_pos({})


def hierarchy_sc(G, root, node_sets):
    def make_sc_tree(sc_G, node=root, parent=None):
        sc_G.add_node(node)
        neighbors = G.neighbors(node)
        if parent is not None:
            sc_G.add_edge(parent, node)
            neighbors.remove(parent)
        old_node = node
        for n in node_sets[int(node)]:
            new_node = str(node) + '_' + str(n)
            sc_G.add_node(new_node)
            sc_G.add_edge(old_node, new_node)
            old_node = new_node
        for neighbor in neighbors:
            sc_G = make_sc_tree(sc_G, neighbor, node)
        return sc_G
    return make_sc_tree(nx.Graph())


def zoom(ax, xy='x', factor=1):
    """Zoom into axis.

    Parameters
    ----------
    """
    limits = ax.get_xlim() if xy == 'x' else ax.get_ylim()
    new_limits = (0.5*(limits[0] + limits[1])
                  + 1./factor * np.array((-0.5, 0.5)) * (limits[1] - limits[0]))
    if xy == 'x':
        ax.set_xlim(new_limits)
    else:
        ax.set_ylim(new_limits)


def get_ax_size(ax, fig):
    """Get axis size

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    fig : matplotlib.Figure
        Figure.
    """
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi


def axis_to_data(ax, width):
    """For a width in axis coordinates, return the corresponding in data
    coordinates.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    width : float
        Width in xaxis coordinates.
    """
    xlim = ax.get_xlim()
    widthx = width*(xlim[1] - xlim[0])
    ylim = ax.get_ylim()
    widthy = width*(ylim[1] - ylim[0])
    return 0.5*(widthx + widthy)


def axis_to_data_points(ax, points_axis):
    """Map points in axis coordinates to data coordinates.

    Uses matplotlib.transform.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    points_axis : np.array
        Points in axis coordinates.
    """
    axis_to_data = ax.transAxes + ax.transData.inverted()
    return axis_to_data.transform(points_axis)


def data_to_axis_points(ax, points_data):
    """Map points in data coordinates to axis coordinates.

    Uses matplotlib.transform.

    Parameters
    ----------
    ax : matplotlib.axis
        Axis object from matplotlib.
    points_axis : np.array
        Points in data coordinates.
    """
    data_to_axis = axis_to_data.inverted()
    return data_to_axis(points_data)
