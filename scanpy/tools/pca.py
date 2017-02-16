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
from .. import preprocess as pp
from ..ann_data import AnnData

def pca(adata_or_X, nr_comps=10):
    """
    Embed data using PCA.

    Parameters
    ----------
    adata_or_X : dict containing
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
    isadata = isinstance(adata_or_X, AnnData)
    if isadata:
        X = adata_or_X.X
        adata = adata_or_X
    else:
        X = adata_or_X
    if isadata and 'Xpca' in adata and adata['Xpca'].shape[1] > nr_comps:
        Y = adata['Xpca']
    else:
        Y = pp.pca(X, nr_comps)
    if isadata:
        return {'type': 'pca', 'Y': Y}
    else:
        return Y

def plot(dplot, adata,
         smp=None,
         names=None,
         comps='1,2',
         cont=None,
         layout='2d',
         legendloc='right margin',
         cmap=None,
         right_margin=None,
         size=3):
    """
    Scatter plots.

    Parameters
    ----------
    dplot : dict
        Dict returned by plotting tool.
    adata : AnnData
        Annotated data matrix.
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
    legendloc : {'right margin', see matplotlib.legend}, optional (default: 'right margin')
         Options for keyword argument 'loc'.
    cmap : str (default: continuous: viridis/ categorical: finite palette)
         String denoting matplotlib color map.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    """
    from .. import plotting as plott
    plott.plot_tool(dplot, adata,
                    smp,
                    names,
                    comps,
                    cont,
                    layout,
                    legendloc,
                    cmap,
                    right_margin,
                    size=size,
                    # defined in plotting
                    subtitles=['PCA'],
                    component_name='PC')

