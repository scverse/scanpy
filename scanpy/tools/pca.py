# coding: utf-8
"""
Principal Component Analysis
============================

From package Scanpy (https://github.com/falexwolf/scanpy).
Written in Python 3 (compatible with 2).
Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).

Notes
-----
There are two PCA versions, which are automatically chosen
- sklearn.decomposition.PCA
- function _pca_fallback
"""

from collections import OrderedDict as odict
import numpy as np
from ..compat.matplotlib import pyplot as pl
from .. import settings as sett
from .. import plotting as plott
from .. import utils
from ..tools import preprocess

def pca(ddata, n_components=10):
    """
    Embed data using PCA.

    Parameters
    ----------
    ddata : dict containing
        X : np.ndarray
            Data array, rows store observations, columns variables.
    n_components : int, optional (default: 2)
        Number of PCs.

    Returns
    -------
    dtsne : dict containing
        Y : np.ndarray
            PCA representation of the data.
    """
    X = ddata['X']
    Y = preprocess.pca(X, n_components=n_components)
    return {'type': 'pca', 'Y': Y}

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
    from numpy import array
    comps = array(params['comps'].split(',')).astype(int) - 1
    # highlights
    highlights = []
    if False:
        if 'highlights' in ddata:
            highlights = ddata['highlights']
    # base figure
    axs = plott.scatter(dpca['Y'][:, comps],
                        subtitles=['PCA'],
                        component_name='PC',
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

