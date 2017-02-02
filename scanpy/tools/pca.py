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

def plot(dplot, ddata,
         comps='1,2,3',
         layout='2d',
         legendloc='lower right',
         cmap='jet',
         adjust_right=0.75): # consider changing to 'viridis'
    """
    Plot the results of a DPT analysis.

    Parameters
    ----------
    dplot : dict
        Dict returned by diffmap tool.
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
    from .. import plotting as plott
    plott.plot_tool(dplot, ddata,
                    comps,
                    layout,
                    legendloc,
                    cmap,
                    adjust_right,
                    # defined in plotting
                    subtitles=['PCA'],
                    component_name='PC')
