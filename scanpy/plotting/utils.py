# Authors: F. Alex Wolf <http://falexwolf.de>
#          P. Angerer

import os
import numpy as np
from ..compat.matplotlib import pyplot as pl
from matplotlib import rcParams, ticker, is_interactive
from matplotlib.figure import SubplotParams as sppars
from cycler import Cycler, cycler

from .. import sett

def init_plotting_params():
    """Init default plotting parameters.

    Is called when importing sc.plotting.
    """
    # figure
    rcParams['figure.figsize'] = (4, 4)
    rcParams['figure.subplot.left'] = 0.18
    rcParams['figure.subplot.right'] = 0.96
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.top'] = 0.91

    rcParams['lines.linewidth'] = 1.5
    rcParams['lines.markersize'] = 6
    rcParams['lines.markeredgewidth'] = 1
    # font
    rcParams['font.sans-serif'] = ['Arial',
                                   'Helvetica',
                                   'DejaVu Sans',
                                   'Bitstream Vera Sans',
                                   'sans-serif']
    fontsize = 14
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = 0.92 * fontsize
    rcParams['axes.titlesize'] = fontsize
    # legend
    rcParams['legend.numpoints'] = 1
    rcParams['legend.scatterpoints'] = 1
    rcParams['legend.handlelength'] = 0.5
    rcParams['legend.handletextpad'] = 0.4
    # resolution of png output
    rcParams['savefig.dpi'] = 400
    # color palette (is default in matplotlib 2.0 anyway)
    # see 'category20' on https://github.com/vega/vega/wiki/Scales#scale-range-literals
    rcParams['axes.prop_cycle'] = cycler(
        color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
               '#9467bd', '#8c564b', '#e377c2',  # '#7f7f7f', remove grey
               '#bcbd22', '#17becf',
               '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
               '#c5b0d5', '#c49c94', '#f7b6d2',  # '#c7c7c7', remove grey
               '#dbdb8d', '#9edae5'])

    if 'DISPLAY' not in os.environ:
        sett.m(0, '    setting `sett.savefigs = True`')
        sett.savefigs = True

    # we show plots instead of saving when running interactively, including in a
    # jupyter notebook/terminal
    sett.savefigs = not is_interactive()


def default_pal(pal=None):
    if pal is None:
        return rcParams['axes.prop_cycle']
    elif not isinstance(pal, Cycler):
        return cycler(color=pal)


def scatter_group(ax, name, imask, adata, Y, layout='2d', size=3):
    """
    Plot group using representation of data Y.
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
    data = [Y[mask, 0], Y[mask, 1]]
    if layout == '3d':
        data.append(Y[mask, 2])
    ax.scatter(*data,
               c=color,
               edgecolors='face',
               s=size,
               label=adata.add[name + '_names'][imask])

def scatter_base(Y,
                 colors='blue',
                 highlights=[],
                 highlights_labels=[],
                 title='',
                 right_margin=None,
                 layout='2d',
                 titles=None,
                 component_name='DC',
                 component_indexnames=[1, 2, 3],
                 axis_labels=None,
                 colorbars=[False],
                 sizes=[1],
                 cmap='viridis',
                 show_ticks=True):
    """Plot scatter plot of data.

    Parameters
    ----------
    Y : np.ndarray
        Data array.
    layout : str
        Either '2d' or '3d'.
    comps : iterable
        Iterable that stores the component indices.

    Returns
    -------
    axs : matplotlib.axis or list of matplotlib.axis
        Depending on whether supplying a single array or a list of arrays,
        return a single axis or a list of axes.
    """
    from matplotlib import gridspec
    # if we have a single array, transform it into a list with a single array
    avail_layouts = {'2d', '3d'}
    if layout not in avail_layouts:
        raise ValueError('choose layout from', avail_layouts)
    if type(colors) == str:
        colors = [colors]
    if len(sizes) != len(colors):
        if len(sizes) == 1:
            sizes = [sizes[0] for i in range(len(colors))]
    # grid of axes for plotting and legends/colorbars
    if np.any(colorbars) and right_margin is None: right_margin = 0.25
    elif right_margin is None: right_margin = 0.01
    # make a figure with panels len(colors) x 1
    top_offset = 1 - rcParams['figure.subplot.top']
    bottom_offset = 0.08
    left_offset = 0.3  # in units of base_height
    base_height = rcParams['figure.figsize'][1]
    height = base_height
    base_width = rcParams['figure.figsize'][0]
    draw_region_width = base_width - left_offset - top_offset - 0.5  # this is kept constant throughout

    right_margin_factor = (1 + right_margin/(1-right_margin))
    width_without_offsets = right_margin_factor * len(colors) * draw_region_width  # this is the total width that keeps draw_region_width

    right_offset = (len(colors) - 1) * left_offset
    width = width_without_offsets + left_offset + right_offset
    left_offset_frac = left_offset / width
    right_offset_frac = 1 - (len(colors) - 1) * left_offset_frac

    figsize = (width, height)
    fig = pl.figure(figsize=figsize,
                    subplotpars=sppars(left=0, right=1, bottom=bottom_offset))
    gs = gridspec.GridSpec(nrows=1,
                           ncols=2*len(colors),
                           width_ratios=[r for i in range(len(colors))
                                         for r in [1-right_margin, right_margin]],
                           left=left_offset_frac,
                           right=right_offset_frac,
                           wspace=0)
    pos = gs.get_grid_positions(fig)
    fig.suptitle(title)
    count = 1
    bool3d = True if layout == '3d' else False
    axs = []
    for icolor, color in enumerate(colors):
        # set up panel
        if layout == 'unfolded 3d' and count != 3:
            ax = fig.add_subplot(2, 2, count)
            bool3d = False
        elif layout == 'unfolded 3d' and count == 3:
            ax = fig.add_subplot(2, 2, count,
                                 projection='3d')
            bool3d = True
        elif layout == '2d':
            ax = pl.subplot(gs[2*(count-1)])
        elif layout == '3d':
            ax = pl.subplot(gs[2*(count-1)], projection='3d')
        if not bool3d:
            data = Y[:, 0], Y[:, 1]
        else:
            data = Y[:, 0], Y[:, 1], Y[:, 2]
        # do the plotting
        if type(color) != str or color != 'white':
            sct = ax.scatter(*data,
                             c=color,
                             edgecolors='face',
                             s=sizes[icolor],
                             cmap=cmap)
        if colorbars[icolor]:
            pos = gs.get_grid_positions(fig)
            left = pos[2][2*(count-1)+1]
            bottom = pos[0][0]
            width = 0.006 * draw_region_width
            # print(0.2*(pos[3][2*(count-1)+1] - left))
            height = pos[1][0] - bottom
            # again shift to left
            left = pos[3][2*(count-1)] + (1 if layout == '3d' else 0.2) * width
            rectangle = [left, bottom, width, height]
            ax_cb = fig.add_axes(rectangle)
            cb = pl.colorbar(sct, format=ticker.FuncFormatter(ticks_formatter),
                             cax=ax_cb)
        # set the titles
        if titles is not None:
            ax.set_title(titles[icolor])
        # output highlighted data points
        for iihighlight, ihighlight in enumerate(highlights):
            data = [Y[ihighlight, 0]], [Y[ihighlight, 1]]
            if bool3d:
                data = [Y[ihighlight, 0]], [Y[ihighlight, 1]], [Y[ihighlight, 2]]
            ax.scatter(*data, c='black',
                       facecolors='black', edgecolors='black',
                       marker='x', s=40, zorder=20)
            highlight = (highlights_labels[iihighlight] if
                         len(highlights_labels) > 0
                         else str(ihighlight))
            # the following is a Python 2 compatibility hack
            ax.text(*([d[0] for d in data] + [highlight]), zorder=20)
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
            if bool3d: ax.set_zticks([])
        # scale limits to match data
        ax.autoscale_view()
        axs.append(ax)
        count += 1
    # set default axis_labels
    if axis_labels is None:
        if layout == '2d':
            axis_labels = [[component_name + str(i) for i in idcs]
                         for idcs in
                         [component_indexnames for iax in range(len(axs))]]
        elif layout == '3d':
            axis_labels = [[component_name + str(i) for i in idcs]
                         for idcs in
                         [component_indexnames for iax in range(len(axs))]]
        elif layout == 'unfolded 3d':
            axis_labels = [[component_name
                         + str(component_indexnames[i-1]) for i in idcs]
                         for idcs in [[2, 3], [1, 2], [1, 2, 3], [1, 3]]]
    else:
        axis_labels = [[axis_labels[0], axis_labels[1]] for i in range(len(axs))]
    # set axis_labels
    bool3d = True if layout == '3d' else False
    for iax, ax in enumerate(axs):
        if layout == 'unfolded 3d' and iax != 2:
            bool3d = False
        elif layout == 'unfolded 3d' and iax == 2:
            bool3d = True
        if axis_labels is not None:
            ax.set_xlabel(axis_labels[iax][0])
            ax.set_ylabel(axis_labels[iax][1])
            if bool3d:
                # shift the label closer to the axis
                ax.set_zlabel(axis_labels[iax][2], labelpad=-7)
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
        return ('%.3f'%(x)).rstrip('0').rstrip('.')


def pimp_axis(ax):
    """Remove trailing zeros.
    """
    ax.set_major_formatter(ticker.FuncFormatter(ticks_formatter))


def scale_to_zero_one(x):
    """Take some 1d data and scale it so that min matches 0 and max 1.
    """
    xscaled = x - np.min(x)
    xscaled /= np.max(xscaled)
    return xscaled


def zoom(ax,xy='x',factor=1):
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


def get_ax_size(ax,fig):
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


def axis_to_data(ax,width):
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


def axis_to_data_points(ax,points_axis):
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
