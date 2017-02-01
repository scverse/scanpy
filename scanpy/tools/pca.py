# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Principal Component Analysis

Notes
-----
There are two PCA versions, which are automatically chosen
- sklearn.decomposition.PCA
- function _pca_fallback
"""

from ..compat.matplotlib import pyplot as pl
from .. import settings as sett
from .. import plotting as plott
from .. import utils
from ..preprocess import preprocess as pp

def pca(ddata_or_X, nr_comps=10):
    """
    Embed data using PCA.

    Parameters
    ----------
    ddata_or_X : dict containing
        X : np.ndarray
            Data array, rows store observations, columns variables.
    nr_comps : int, optional (default: 2)
        Number of PCs.

    Returns
    -------
    dtsne : dict containing
        Y : np.ndarray
            PCA representation of the data.
    """
    if isinstance(ddata_or_X, dict):      
        X = ddata_or_X['X']
        isddata = True
    else:
        X = ddata_or_X
        isddata = False
    Y = pp(X, 'pca', nr_comps)
    if isddata:
        return {'type': 'pca', 'Y': Y}
    else:
        return Y

def plot(dpca, ddata,
         comps='1,2,3',
         layout='2d',
         legendloc='lower right',
         cmap='jet',
         adjust_right=0.75): # consider changing to 'viridis'
    """
    Plot the results of a DPT analysis.

    Parameters
    ----------
    dpca : dict
        Dict returned by PCA tool.
    ddata : dict
        Data dictionary.
    comps : str
         String in the form "comp1,comp2,comp3".
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right') 
         Options for keyword argument 'loc'.
    cmap : str, optional (default: jet)
         String denoting matplotlib color map. 
    """
    params = locals(); del params['ddata']; del params['dpca']
    # highlights
    highlights = []
    if False:
        if 'highlights' in ddata:
            highlights = ddata['highlights']
    # base figure
    from numpy import array
    comps = array(params['comps'].split(',')).astype(int) - 1
    try:
        Y = dpca['Y'][:, comps]
    except IndexError:
        sett.mi('IndexError: Only computed', dpca['Y'].shape[1], ' components')
        sett.mi('--> recompute using scanpy exkey pca -p nr_comps YOUR_NR')
        from sys import exit
        exit(0)
    axs = plott.scatter(Y,
                        subtitles=['PCA'],
                        component_name='PC',
                        component_indexnames=comps + 1,
                        layout=params['layout'],
                        c='grey',
                        highlights=highlights,
                        cmap=params['cmap'])
    # annotated groups
    if 'groupmasks' in ddata:
        for igroup, group in enumerate(ddata['groupmasks']):
            plott.group(axs[0], igroup, ddata, dpca['Y'][:, comps], params['layout'])
        axs[0].legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
        # right margin
        pl.subplots_adjust(right=params['adjust_right'])

    if sett.savefigs:
        pl.savefig(sett.figdir+dpca['writekey']+'.'+sett.extf)
    elif sett.autoshow:
        pl.show()

