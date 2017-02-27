# Authors: F. Alex Wolf <http://falexwolf.de>
#          P. Angerer        
"""
Plotting
"""

import numpy as np
from cycler import Cycler, cycler

from .compat.matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib import ticker
from matplotlib.figure import SubplotParams as sppars
from . import settings as sett

#--------------------------------------------------------------------------------
# Scanpy Plotting Functions
#--------------------------------------------------------------------------------

def savefig(writekey):
    if sett.savefigs:
        filename = sett.figdir + writekey + '.' + sett.extf
        sett.m(0, 'saving figure to file', filename)
        pl.savefig(filename)

def default_pal(pal=None):
    if pal is None:
        return rcParams['axes.prop_cycle']
    elif not isinstance(pal, Cycler):
        return cycler(color=pal)

def plot_tool(adata,
              toolkey,
              basis='pca',
              smp=None,
              names=None,
              comps=None,
              cont=None,
              layout='2d',
              legendloc='right margin',
              cmap=None,
              pal=None,
              right_margin=None,
              size=3,
              subtitles=None):
    """
    Scatter plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    toolkey : str
        Toolkey to use for generating filename.
    basis : {'pca', 'tsne', 'diffmap'}
        String that denotes a plotting tool.
    smp : str, optional (default: first annotation)
        Sample/Cell annotation for coloring in the form "ann1,ann2,...". String
        annotation is plotted assuming categorical annotation, float and integer
        annotation is plotted assuming continuous annoation. Option 'cont'
        allows to switch between these default choices.
    names : str, optional (default: all names in smp)
        Allows to restrict groups in sample annotation (smp) to a few.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
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
    """
    # compute components
    from numpy import array
    if comps is None:
        comps = '1,2' if '2d' in layout else '1,2,3'
    comps = array(comps.split(',')).astype(int) - 1
    smps = [None] if smp is None else smp.split(',') if isinstance(smp, str) else smp
    names = None if names is None else names.split(',') if isinstance(names, str) else names
    # highlights
    highlights = adata['highlights'] if 'highlights' in adata else []
    try:
        Y = adata['X_' + basis][:, comps]
    except KeyError:
        sett.mi('--> compute the basis using plotting tool', basis, 'first')
        raise

    pal = default_pal(pal)

    component_name = ('DC' if basis == 'diffmap'
                      else 'Spring' if basis == 'spring'
                      else 'tSNE' if basis == 'tsne'
                      else 'PC')

    colors = []
    categoricals = []
    colorbars = []
    for ismp, smp in enumerate(smps):
        c = 'grey'
        categorical = False
        continuous = False
        if len(adata.smp_keys()) > 0:
            if smp is None:
                for smp_ in adata.smp_keys():
                    if smp_ not in smps:
                        smp = smp_
                        smps[ismp] = smp
            # test whether we have categorial or continuous annotation
            if smp in adata.smp_keys():
                if adata.smp[smp].dtype.char in ['S', 'U']:
                    categorical = True
                    if cont is True:
                        c = adata.smp[smp]
                else:
                    continuous = True
                    c = adata.smp[smp]
                sett.m(0, '... coloring according to', smp)
            # coloring according to gene expression
            elif smp in adata.var_names:
                c = adata.X[:, np.where(smp==adata.var_names)[0][0]]
                continuous = True
                sett.m(0, '... coloring according to expression of gene', smp)
            else:
                raise ValueError('"' + smp + '" is invalid!'
                                 + ' specify valid sample annotation, one of '
                                 + str(adata.smp_keys()) + ' or a gene name '
                                 + str(adata.var_names))
        if cont is not None:
            categorical = not cont
            continuous = cont
        colors.append(c)
        colorbars.append(True if continuous else False)
        if categorical:
            categoricals.append(ismp)

    if right_margin is None and legendloc == 'right margin':
        right_margin = 0.24
    if subtitles is None and smps[0] is not None:
        subtitles = [smp.replace('_', ' ') for smp in smps]

    axs = scatter(Y,
                  subtitles=subtitles,
                  component_name=component_name,
                  component_indexnames=comps + 1,
                  layout=layout,
                  c=colors,
                  highlights=highlights,
                  colorbars=colorbars,
                  right_margin=right_margin,
                  cmap='viridis' if cmap is None else cmap,
                  s=size-2)

    for ismp in categoricals:
        smp = smps[ismp]
        if (smp != 'groups' and 'groups_names' in adata
            and len(np.setdiff1d(adata['groups_names'], adata[smp + '_names']))
                < len(adata['groups_names'])):
            # if there is a correspondence between smp and the 'groups' defined
            # in adata, that is, if smp has corresponding categories with those
            # in adata['groups_names']
            adata[smp + '_colors'] = pal[:len(adata['groups_names'])].by_key()['color']
        elif not smp + '_colors' in adata:
            adata[smp + '_colors'] = pal[:len(adata[smp + '_names'])].by_key()['color']
        for iname, name in enumerate(adata[smp + '_names']):
            if (names is None or (names != None and name in names)):
                group(axs[ismp], smp, iname, adata, Y, layout, size)
        # for the last panel, let's put the legend outside panel
        if legendloc == 'right margin':
            axs[ismp].legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
        elif legendloc != 'none':
            axs[ismp].legend(frameon=False, loc=legendloc)

    writekey = sett.basekey + '_' + toolkey
    writekey += '_' + ('-'.join(smps) if smps[0] is not None else '') + sett.plotsuffix
    savefig(writekey)
    if not sett.savefigs and sett.autoshow:
        pl.show()

def group(ax, name, imask, adata, Y, layout='2d', size=3):
    """
    Plot group using representation of data Y.
    """
    if name + '_masks' in adata:
        mask = adata[name + '_masks'][imask]
    else:
        if adata[name + '_names'][imask] in adata.smp[name]:
            mask = adata[name + '_names'][imask] == adata.smp[name]
        else:
            mask = str(imask) == adata.smp[name]
    color = adata[name + '_colors'][imask]
    if not isinstance(color[0], str):
        from matplotlib.colors import rgb2hex
        color = rgb2hex(adata[name + '_colors'][imask])
    data = [Y[mask, 0], Y[mask, 1]]
    if layout == '3d':
        data.append(Y[mask, 2])
    ax.scatter(*data,
               c=color,
               edgecolors='face',
               s=size,
               alpha=1,
               label=adata[name + '_names'][imask])

def timeseries(X, **kwargs):
    """
    Plot X. See timeseries_subplot.
    """
    pl.figure(figsize=(2*4,4),
              subplotpars=sppars(left=0.12,right=0.98,bottom=0.13))
    timeseries_subplot(X, **kwargs)

def timeseries_subplot(X,
                       c=None,
                       varnames=(),
                       highlightsX=(),
                       xlabel='segments / pseudotime order',
                       ylabel='gene expression',
                       yticks=None,
                       xlim=None,
                       legend=True,
                       pal=None,
                       cmap='viridis'):
    """
    Plot X. Call this with:
    X with one column, c categorical
    X with one column, c continuous
    X with n columns, c is of length n
    """
    
    if c is not None:
        use_cmap = isinstance(c[0], float)
    pal = default_pal(pal)
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
    for ih,h in enumerate(highlightsX):
        pl.plot([h,h],[ylim[0],ylim[1]],
                '--',color='black')
    pl.ylim(ylim)
    if xlim is not None:
        pl.xlim(xlim)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if yticks is not None:
        pl.yticks(yticks)
    if len(varnames) > 0 and legend==True:
        pl.legend(frameon=False)

def timeseries_as_heatmap(X, varnames=None, highlightsX = None, cmap='viridis'):
    """
    Plot timeseries as heatmap.

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
        varnames = varnames[:,0]

    # transpose X
    X = X.T
    minX = np.min(X)

    # insert space into X
    if False:
        # generate new array with highlightsX
        space = 10 # integer
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

    fig = pl.figure(figsize=(1.5*4,2*4))
    im = pl.imshow(np.array(X,dtype=np.float_), aspect='auto',
              interpolation='nearest', cmap=cmap)
    pl.colorbar(shrink=0.5)
    pl.yticks(range(X.shape[0]), varnames)
    for ih,h in enumerate(highlightsX):
        pl.plot([h,h], [0,X.shape[0]],
                '--', color='black')
    pl.xlim([0, X.shape[1]-1])
    pl.ylim([0, X.shape[0]-1])
    pl.xlabel('segments / pseudotime order')
    # pl.tight_layout()

def scatter(Ys,
            c='blue',
            highlights=[],
            highlights_labels=[],
            title='',
            right_margin=None,
            layout='2d',
            subtitles=None,
            component_name='DC',
            component_indexnames=[1, 2, 3],
            axlabels=None,
            colorbars=[False],
            **kwargs):
    """
    Plot scatter plot of data.

    Parameters
    ----------
    Ys : np.ndarray or list of np.ndarray
        Single data array or list of data arrays. Rows store observations,
        columns store variables. For example X, or phi or [phi,psi]. Arrays must
        be of dimension ndim=2.
    layout : str
        Choose from '2d', '3d' and 'unfolded 3d', default '2d'.
    comps : iterable
        Iterable that stores the component indices.

    Returns
    -------
    axs : matplotlib.axis or list of matplotlib.axis
        Depending on whether supplying a single array or a list of arrays,
        return a single axis or a list of axes.
    """
    # if we have a single array, transform it into a list with a single array
    avail_layouts = ['2d', '3d', 'unfolded 3d']
    if layout not in avail_layouts:
        raise ValueError('choose layout from',avail_layouts)
    colors = c
    if type(Ys) == np.ndarray:
        Ys = [Ys]
    if len(colors) == len(Ys[0]) or type(colors) == str:
        colors = [colors]
    # make a figure with panels len(colors) x len(Ys)
    figsize = (4*len(colors), 4*len(Ys))
    # checks
    if layout == 'unfolded 3d':
        if len(Ys) != 1:
            raise ValueError('use single 3d array')
        if len(colors) > 1:
            raise ValueError('choose a single color')
        figsize = (4*2, 4*2)
        Y = Ys[0]
        Ys = [Y[:,[1,2]], Y[:,[0,1]], Y, Y[:,[0,2]]]
    # try importing Axes3D
    if '3d' in layout:
        from mpl_toolkits.mplot3d import Axes3D

    fig = pl.figure(figsize=figsize,
                    subplotpars=sppars(left=0, right=1, bottom=0.08))

    from matplotlib import gridspec
    # grid of axes for plotting and legends/colorbars
    if np.any(colorbars) and right_margin is None:
        right_margin = 0.25
    elif right_margin is None:
        right_margin = 0.01

    from inspect import getcallargs
    gs = gridspec.GridSpec(nrows=len(Ys),
                           ncols=2*len(colors),
                           width_ratios=[r for i in range(len(colors))
                                         for r in [1-right_margin, right_margin]],
                           left=0.08/len(colors),
                           right=1-(len(colors)-1)*0.08/len(colors),
                           wspace=0)
    fig.suptitle(title)
    count = 1
    bool3d = True if layout == '3d' else False
    axs = []
    for Y in Ys:
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
                data = Y[:,0], Y[:,1]
            else:
                data = Y[:,0], Y[:,1], Y[:,2]
            # do the plotting
            if type(color) != str or color != 'white':
                sct = ax.scatter(*data,
                                 c=color,
                                 edgecolors='face',
                                 **kwargs)
            if colorbars[icolor]:
                pos = gs.get_grid_positions(fig)
                left = pos[2][2*(count-1)+1]
                bottom = pos[0][0]
                width = 0.2*(pos[3][2*(count-1)+1] - left)
                height = pos[1][0] - bottom
                # again shift to left
                left = pos[3][2*(count-1)] + 0.1*width
                rectangle = [left, bottom, width, height]
                ax_cb = fig.add_axes(rectangle)
                cb = pl.colorbar(sct, format=ticker.FuncFormatter(ticks_formatter),
                                 cax=ax_cb)
            # set the subtitles
            if subtitles is not None:
                ax.set_title(subtitles[icolor])
            # output highlighted data points
            for iihighlight,ihighlight in enumerate(highlights):
                data = [Y[ihighlight,0]], [Y[ihighlight,1]]
                if bool3d:
                    data = [Y[ihighlight,0]], [Y[ihighlight,1]], [Y[ihighlight,2]]
                ax.scatter(*data, c='black',
                           facecolors='black', edgecolors='black',
                           marker='x', s=40, zorder=20)
                highlight = (highlights_labels[iihighlight] if
                             len(highlights_labels) > 0
                             else str(ihighlight))
                # the following is a Python 2 compatibility hack
                ax.text(*([d[0] for d in data]+[highlight]),zorder=20)
            ax.set_xticks([]); ax.set_yticks([])
            if bool3d:
                ax.set_zticks([])
            axs.append(ax)
            count += 1
    # set default axlabels
    if axlabels is None:
        if layout == '2d':
            axlabels = [[component_name + str(i) for i in idcs]
                         for idcs in
                         [component_indexnames for iax in range(len(axs))]]
        elif layout == '3d':
            axlabels = [[component_name + str(i) for i in idcs]
                         for idcs in
                         [component_indexnames for iax in range(len(axs))]]
        elif layout == 'unfolded 3d':
            axlabels = [[component_name
                         + str(component_indexnames[i-1]) for i in idcs]
                         for idcs in [[2, 3], [1, 2], [1, 2, 3], [1, 3]]]
    # set axlabels
    bool3d = True if layout == '3d' else False
    for iax,ax in enumerate(axs):
        if layout == 'unfolded 3d' and iax != 2:
            bool3d = False
        elif layout == 'unfolded 3d' and iax == 2:
            bool3d = True
        if axlabels is not None:
            ax.set_xlabel(axlabels[iax][0])
            ax.set_ylabel(axlabels[iax][1])
            if bool3d:
                # shift the label closer to the axis
                ax.set_zlabel(axlabels[iax][2],labelpad=-7)
    return axs

# the following doesn't use a fixed grid, therefore panels have
# different sizes depending on whether there is a suplot or not
def _scatter_auto_grid(Ys,
             layout='2d',
             subtitles=['pseudotime', 'segments', 'experimental labels'],
             c='blue',
             highlights=[],
             highlights_labels=[],
             title='',
             colorbars=[False],
             **kwargs):
    # if we have a single array, transform it into a list with a single array
    avail_layouts = ['2d', '3d', 'unfolded 3d']
    if layout not in avail_layouts:
        raise ValueError('choose layout from',avail_layouts)
    colors = c
    if type(Ys) == np.ndarray:
        Ys = [Ys]
    if len(colors) == len(Ys[0]) or type(colors) == str:
        colors = [colors]
    # make a figure with panels len(colors) x len(Ys)
    figsize = (4*len(colors), 4*len(Ys))
    # checks
    if layout == 'unfolded 3d':
        if len(Ys) != 1:
            raise ValueError('use single 3d array')
        if len(colors) > 1:
            raise ValueError('choose a single color')
        figsize = (4*2, 4*2)
        Y = Ys[0]
        Ys = [Y[:,[1,2]], Y[:,[0,1]], Y, Y[:,[0,2]]]
    # try importing Axes3D
    if '3d' in layout:
        from mpl_toolkits.mplot3d import Axes3D
    fig = pl.figure(figsize=figsize,
                    subplotpars=sppars(left=0.07,right=0.98,bottom=0.08))
    fig.suptitle(title)
    count = 1
    bool3d = True if layout == '3d' else False
    axs = []
    for Y in Ys:
        markersize = (2 if Y.shape[0] > 500 else 10)
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
                ax = fig.add_subplot(len(Ys), len(colors), count)
            elif layout == '3d':
                ax = fig.add_subplot(len(Ys), len(colors), count,
                                     projection='3d')
            if not bool3d:
                data = Y[:,0], Y[:,1]
            else:
                data = Y[:,0], Y[:,1], Y[:,2]
            # do the plotting
            if type(color) != str or color != 'white':
                sct = ax.scatter(*data,
                                 c=color,
                                 edgecolors='face',
                                 s=markersize,
                                 **kwargs)
            if colorbars[icolor]:
                cb = pl.colorbar(sct, format=ticker.FuncFormatter(ticks_formatter))
            # set the subtitles
            ax.set_title(subtitles[icolor])
            # output highlighted data points
            for iihighlight,ihighlight in enumerate(highlights):
                data = [Y[ihighlight,0]], [Y[ihighlight,1]]
                if bool3d:
                    data = [Y[ihighlight,0]], [Y[ihighlight,1]], [Y[ihighlight,2]]
                ax.scatter(*data, c='black',
                           facecolors='black', edgecolors='black',
                           marker='x', s=40, zorder=20)
                highlight = (highlights_labels[iihighlight] if
                             len(highlights_labels) > 0
                             else str(ihighlight))
                # the following is a Python 2 compatibility hack
                ax.text(*([d[0] for d in data]+[highlight]),zorder=20)
            ax.set_xticks([]); ax.set_yticks([])
            if bool3d:
                ax.set_zticks([])
            axs.append(ax)
            count += 1
    # scatter.set_edgecolors = scatter.set_facecolors = lambda *args:None
    return axs

def _scatter_single(ax,Y,*args,**kwargs):
    """
    Plot scatter plot of data. Just some wrapper of matplotlib.Axis.scatter.

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
    ax.scatter(Y[:,0],Y[:,1],**kwargs)
    ax.set_xticks([]); ax.set_yticks([])

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
    scoreskey = adata[toolkey + '_scoreskey']
    n_panels = len(adata[toolkey + '_rankings_names'])

    def get_scores(irank):
        allscores = adata[toolkey + '_' + scoreskey][irank]
        scores = allscores[adata[toolkey + '_rankings_geneidcs'][irank, :n_genes]]
        scores = np.abs(scores)
        return scores

    # the limits for the y axis
    ymin = 1e100
    ymax = -1e100
    for irank in range(len(adata[toolkey + '_rankings_names'])):
        scores = get_scores(irank)
        ymin = np.min([ymin,np.min(scores)])
        ymax = np.max([ymax,np.max(scores)])
    ymax += 0.3*(ymax-ymin)

    # number of panels
    if n_panels <= 5:
        n_panels_y = 1
        n_panels_x = n_panels
    else:
        n_panels_y = 2
        n_panels_x = int(n_panels/2+0.5)

    fig = pl.figure(figsize=(n_panels_x*4,n_panels_y*4))
    pl.subplots_adjust(left=0.15,top=0.9,right=0.98,bottom=0.13)

    count = 1
    for irank in range(len(adata[toolkey + '_rankings_names'])):
        fig.add_subplot(n_panels_y,n_panels_x,count)
        scores = get_scores(irank)
        for ig,g in enumerate(adata[toolkey + '_rankings_geneidcs'][irank, :n_genes]):
            marker = (r'\leftarrow' if adata[toolkey + '_zscores'][irank,g] < 0
                                    else r'\rightarrow')
            pl.text(ig,scores[ig],
                    r'$ ' + marker + '$ ' + adata.var_names[g],
                    color = 'red' if adata[toolkey + '_zscores'][irank,g] < 0 else 'green',
                    rotation='vertical',verticalalignment='bottom',
                    horizontalalignment='center',
                    fontsize=8)
        title = adata[toolkey + '_rankings_names'][irank]
        pl.title(title)
        if n_panels <= 5 or count > n_panels_x:
            pl.xlabel('ranking')
        if count == 1 or count == n_panels_x+1:
            pl.ylabel(scoreskey)
        else:
            pl.yticks([])
        pl.ylim([ymin,ymax])
        pl.xlim(-0.9,ig+1-0.1)
        count += 1

def arrows_transitions(ax,X,indices,weight=None):
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
    width = axis_to_data(ax,0.001)
    if X.shape[0] > 300:
        step = 5
        width = axis_to_data(ax,0.0005)
    if X.shape[0] > 500:
        step = 30
        width = axis_to_data(ax,0.0001)
    head_width = 10*width
    for ix,x in enumerate(X):
        if ix%step == 0:
            X_step = X[indices[ix]] - x
            # don't plot arrow of length 0
            for itrans in range(X_step.shape[0]):
                alphai = 1
                widthi = width
                head_widthi = head_width
                if weight is not None:
                    alphai *= weight[ix,itrans]
                    widthi *= weight[ix,itrans]
                if np.any(X_step[itrans,:1]):
                    ax.arrow(x[0], x[1],
                             X_step[itrans,0], X_step[itrans,1],
                             length_includes_head=True,
                             width=widthi,
                             head_width=head_widthi,
                             alpha=alphai,
                             color='grey')
#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

def ticks_formatter(x, pos):
    # pretty scientific notation
    if False:
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
    else:
        return ('%.3f'%(x)).rstrip('0').rstrip('.')

def pimp_axis(ax):
    """
    Remove trailing zeros.
    """
    ax.set_major_formatter(ticker.FuncFormatter(ticks_formatter))

def scale_to_zero_one(x):
    """
    Take some 1d data and scale it so that min matches 0 and max 1.
    """
    xscaled = x - np.min(x)
    xscaled /= np.max(xscaled)
    return xscaled

def zoom(ax,xy='x',factor=1):
    """
    Zoom into axis.

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
    """
    Get axis size

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
    """
    For a width in axis coordinates, return the corresponding in data
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
    """
    Map points in axis coordinates to data coordinates.

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

def data_to_axis_points(ax,points_data):
    """
    Map points in data coordinates to axis coordinates.

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

#--------------------------------------------------------------------------------
# Global Plotting Variables
#--------------------------------------------------------------------------------

# standard linewidth
lw0 = rcParams['lines.linewidth']
# list of markers
ml = ['o', 's', '^', 'd']

# init default values for rcParams
def init_fig_params():
    # args figure
    rcParams['figure.figsize'] = (5,4)
    rcParams['figure.subplot.left'] = 0.18
    rcParams['figure.subplot.right'] = 0.96
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.top'] = 0.91

    rcParams['lines.linewidth'] = 1.5
    rcParams['lines.markersize'] = 6
    rcParams['lines.markeredgewidth'] = 1

    # args font
    rcParams['font.family'] = ['Arial','Helvetica']
    fontsize = 14
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = 0.92*fontsize
    rcParams['axes.titlesize'] = fontsize

    # legend
    rcParams['legend.numpoints'] = 1
    rcParams['legend.scatterpoints'] = 1
    rcParams['legend.handlelength'] = 0.5
    rcParams['legend.handletextpad'] = 0.4

    # resolution of png output
    rcParams['savefig.dpi'] = 200

    # color palette
    # see 'category20' on 
    # https://github.com/vega/vega/wiki/Scales#scale-range-literals
    rcParams['axes.prop_cycle'] = cycler(color=
                                         ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                                          '#9467bd', '#8c564b', '#e377c2', #'#7f7f7f', remove grey
                                          '#bcbd22', '#17becf',
                                          '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
                                          '#c5b0d5', '#c49c94', '#f7b6d2', #'#c7c7c7', remove grey
                                          '#dbdb8d', '#9edae5'])

init_fig_params()

