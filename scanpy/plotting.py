# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Plotting
"""   

import numpy as np
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
        pl.savefig(sett.figdir + writekey + '.' + sett.extf)

def plot_tool(dplot, ddata,
              rowcat='',
              comps='1,2',
              layout='2d',
              legendloc='lower right',
              cmap='jet',
              adjust_right=0.75,
              subtitles=['one title'],
              component_name='comp'): 
    """
    Plot the results of a DPT analysis.

    Parameters
    ----------
    dplot : dict
        Dict returned by plotting tool.
    ddata : dict
        Data dictionary.
    comps : str
         String in the form "comp1,comp2,comp3".
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: jet)
         String denoting matplotlib color map. 
    """
    params = locals(); del params['ddata']; del params['dplot']
    from numpy import array
    comps = array(params['comps'].split(',')).astype(int) - 1
    # highlights
    highlights = []
    if False:
        if 'highlights' in ddata:
            highlights = ddata['highlights']
    # base figure
    try:
        Y = dplot['Y'][:, comps]
    except IndexError:
        sett.mi('IndexError: Only computed', dplot['Y'].shape[1], ' components')
        sett.mi('--> recompute using scanpy exkey diffmap -p nr_comps YOUR_NR')
        from sys import exit
        exit(0)
    axs = scatter(Y,
                  subtitles=subtitles,
                  component_name=component_name,
                  component_indexnames=comps + 1,
                  layout=params['layout'],
                  c='grey',
                  highlights=highlights,
                  cmap=params['cmap'])
    # rowcategories in ddata
    if 'rowcat' in ddata:
        if rowcat == '':
            # simply take a random key
            rowcat = list(ddata['rowcat'].keys())[0]
        elif rowcat not in ddata:
            print('specify valid row category class')
        # colors for categories
        if not rowcat + '_colors' in ddata:
            ddata[rowcat + '_colors'] = pl.cm.get_cmap(params['cmap'])(
                                                  pl.Normalize()(ddata[k + '_ids']))
        for icat in ddata[rowcat + '_ids']:
            group(axs[0], rowcat, icat, ddata, dplot['Y'][:, comps], params['layout'])
        axs[0].legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
        # right margin
        pl.subplots_adjust(right=params['adjust_right'])

    if sett.savefigs:
        pl.savefig(sett.figdir+dplot['writekey']+'.'+sett.extf)
    elif sett.autoshow:
        pl.show()

def group(ax, name, imask, ddata, Y, layout='2d'):
    """
    Plot group using representation of data Y.
    """
    mask = ddata[name + '_masks'][imask]
    color = ddata[name + '_colors'][imask]
    if not isinstance(color[0], str):
        from matplotlib.colors import rgb2hex
        color = rgb2hex(ddata[name + '_colors'][imask])
    data = [Y[mask,0], Y[mask,1]]
    if layout == '3d':
        data.append(Y[mask,2])
    markersize = 3
    ax.scatter(*data,
               c=color,
               edgecolors='face',
               s=markersize,
               alpha=1,
               label=ddata[name + '_names'][imask])

def timeseries(*args,**kwargs):
    """ 
    Plot X. See timeseries_subplot.
    """
    pl.figure(figsize=(2*4,4),
              subplotpars=sppars(left=0.12,right=0.98,bottom=0.13))
    timeseries_subplot(*args,**kwargs)

def timeseries_subplot(X,
                       varnames=[],
                       highlightsX=[],
                       c = None,
                       xlabel='segments / pseudotime order',
                       ylabel='gene expression',
                       yticks=None,
                       xlim=None,
                       legend=True,
                       cmap='jet'): # consider changing to 'viridis'
    """ 
    Plot X.
    """
    for i in range(X.shape[1]):
        pl.scatter(
            np.arange(X.shape[0]),X[:,i],
            marker='.',
            edgecolor='face',
            s=rcParams['lines.markersize'],
            c=cl[i] if c is None else c,
            label = (varnames[i] if len(varnames) > 0 else ''),
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

def timeseries_as_heatmap(X, varnames=None, highlightsX = None):
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
        Xnew = np.zeros((X.shape[0], X.shape[1]+space*len(highlightsX)))
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
              interpolation='nearest', cmap='viridis')
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
            layout='2d',
            subtitles=['pseudotime', 'segments', 'experimental labels'],
            component_name='DC',
            component_indexnames=[1, 2, 3],
            axlabels=None,
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
    axs = _scatter(Ys,layout=layout,subtitles=subtitles,**kwargs)
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

def _scatter(Ys,
             layout='2d',
             subtitles=['pseudotime', 'segments', 'experimental labels'],
             c='blue',
             highlights=[],
             highlights_labels=[],
             title='', 
             cmap='jet'):
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
        figsize = (4*2,4*2)
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
        for icolor,color in enumerate(colors):
            # set up panel
            if layout == 'unfolded 3d' and count != 3:
                ax = fig.add_subplot(2,2,count)
                bool3d = False
            elif layout == 'unfolded 3d' and count == 3:
                ax = fig.add_subplot(2,2,count,
                                     projection='3d')
                bool3d = True
            elif layout == '2d':
                ax = fig.add_subplot(len(Ys),len(colors),count)
            elif layout == '3d':
                ax = fig.add_subplot(len(Ys),len(colors),count,
                                     projection='3d')
            if not bool3d:
                data = Y[:,0],Y[:,1]
            else:
                data = Y[:,0],Y[:,1],Y[:,2]
            # do the plotting
            if type(color) != str or 'white' != color:
                ax.scatter(*data,
                           c=color,
                           edgecolors='face',
                           s=markersize,
                           cmap=cmap)
            # set the subsubtitles
            if icolor == 0:
                ax.set_title(subtitles[0])
            if icolor == 1:
                ax.set_title(subtitles[1])
            if icolor == 2:
                ax.set_title(subtitles[2])
            # output highlighted data points
            for iihighlight,ihighlight in enumerate(highlights):
                data = [Y[ihighlight,0]],[Y[ihighlight,1]]
                if bool3d:
                    data = [Y[ihighlight,0]],[Y[ihighlight,1]],[Y[ihighlight,2]]
                ax.scatter(*data,c='black',
                           facecolors='black',edgecolors='black', 
                           marker='x',s=40,zorder=20)
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

def ranking(drankings, ddata, nr=20):
    """ 
    Plot ranking of genes

    Parameters
    ----------
    drankings : dict containing
        scoreskey : str
            Key to identify scores.
        scores : np.ndarray
            Array of scores for genes according to which to rank them.
    ddata : dict
        Data dictionary.
    nr : int
        Number of genes.
    """

    scoreskey = drankings['scoreskey']

    # one panel for each ranking
    nr_panels = len(drankings['testnames'])
    # maximal number of genes that is shown
    nr_genes = nr

    def get_scores(irank):
        allscores = drankings[scoreskey][irank]
        allscores = np.ma.masked_invalid(allscores)
        scores = allscores[drankings['genes_sorted'][irank,:nr_genes]]
        scores = np.abs(scores)
        return scores

    # the limits for the y axis
    ymin = 1e100
    ymax = -1e100
    for irank in range(len(drankings['testnames'])):
        scores = get_scores(irank)
        ymin = np.min([ymin,np.min(scores)])
        ymax = np.max([ymax,np.max(scores)])
    ymax += 0.3*(ymax-ymin)

    if nr_panels <= 5:
        nr_panels_y = 1
        nr_panels_x = nr_panels
    else:
        nr_panels_y = 2
        nr_panels_x = int(nr_panels/2+0.5)

    fig = pl.figure(figsize=(nr_panels_x*4,nr_panels_y*4))
    pl.subplots_adjust(left=0.15,top=0.9,right=0.98,bottom=0.13)

    count = 1
    for irank in range(len(drankings['testnames'])):
        fig.add_subplot(nr_panels_y,nr_panels_x,count)
        scores = get_scores(irank)
        for ig,g in enumerate(drankings['genes_sorted'][irank,:nr_genes]):
            marker = (r'\leftarrow' if drankings['zscores'][irank,g] < 0 
                                    else r'\rightarrow')
            pl.text(ig,scores[ig],
                    r'$ ' + marker + '$ '+ddata['colnames'][g],
                    color = 'red' if drankings['zscores'][irank,g] < 0 else 'green',
                    rotation='vertical',verticalalignment='bottom',
                    horizontalalignment='center',
                    fontsize=8)
        title = drankings['testnames'][irank]
        pl.title(title)
        if nr_panels <= 5 or count > nr_panels_x:
            pl.xlabel('ranking')
        if count == 1 or count == nr_panels_x+1: 
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
    
# color dict
cDict = {
    '0':'black',
    # red
    '1':'#8B0000',  # darkred
    '11':'#EE0000', # red
    # orange
     '2':'#DD2F00',   # darkorange1                                        
     '21':'#DD8C00',  # darkorange                                         
     '22':'#FF7F50',  # darkorange                            
     '23':'#FFA500',  # orange                                             
     '24':'#FFBB00',  # orange                                             
    # green
    '3':'#006400',   # darkgreen
    '31':'#556B2F',  # DarkOliveGreen
    '32':'#228B22',  # forestgreen
    '33':'#66CD00',  # chartreuse3
    '34':'#42CD42',  # limegreen
    '35':'#7CFC00',  # LawnGreen
    # blue
    '4':'#00008B',   # darkblue
    '41':'#104E8B',  # DodgerBlue4
    '42':'#1874CD',  # dodgerblue3
    '43':'#1E90FF',  # DodgerBlue
    '44':'#1E90FF',  # DodgerBlue
    '45':'#00BFFF',  # DeepSkyBlue
    # violett
    '5':'#68228B',  # DarkOrchid4
    '51':'#9932CC', # darkorchid
    '52':'#8B008B', # darkmagenta
    '53':'#FF34B3', # maroon1
    # yellow
    '6':'#FFA500',   # orange
    '61':'#FFB90F',  # DarkGoldenrod1
    '62':'#FFD700',  # gold
    '63':'#FFFF00',  # yellow
    # other
    '7':'turquoise', # turquoise
    '8':'#212121',   # darkgrey
    '81':'#424242',
    '82':'grey',
    '99':'white',
}
# list of colors
cl = [cDict[k] for k in ['4','3','1','21','42','33','11',
                          '53','23','35','7','81','0']]
# list of colors (rainbow)
clrb = [cDict[k] for k in ['1','11',                 # red (dark to light)
                            '2','21','23',            # orange (dark to light)
                            '61','62',                # yellow (dark to light)
                            '35','34','33','3',       # green (light to dark)
                            '45','43','42','41','4',  # blue (light to dark)
                            '5',                      # violet
                            '8'                       # dark grey 
                            ]]
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

    # resolution of png output
    rcParams['savefig.dpi'] = 200

init_fig_params()

