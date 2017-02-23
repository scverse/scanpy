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
from ..classes.ann_data import AnnData

def pca(adata_or_X, n_comps=10):
    """
    Embed data using PCA.

    Parameters
    ----------
    adata_or_X : AnnData object or np.ndarray
        X : np.ndarray
            Data matrix of shape n_samples x n_variables.
    n_comps : int, optional (default: 10)
        Number of principal components to compute.

    Returns
    -------
    X_pca : np.ndarray
         PCA representation of the data with shape n_variables x n_comps.
         Depending on whether an AnnData or a data matrix has been
         provided, the array is written to AnnData or returned directly.
    """
    isadata = isinstance(adata_or_X, AnnData)
    if isadata:
        X = adata_or_X.X
        adata = adata_or_X
    else:
        X = adata_or_X
    if isadata and 'X_pca' in adata and adata['X_pca'].shape[1] >= n_comps:
        return adata
    else:
        X_pca = pp.pca(X, n_comps)
    if isadata:
        adata['X_pca'] = X_pca
        return adata
    else:
        return X_pca

def plot(adata,
         smp=None,
         names=None,
         comps=None,
         cont=None,
         layout='2d',
         legendloc='right margin',
         cmap=None,
         pal=None,
         right_margin=None,
         size=3):
    """
    Scatter plots.

    Parameters
    ----------
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
    cmap : str (default: 'viridis')
         String denoting matplotlib color map.
    pal : list of str (default: matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
         Colors cycle to use for categorical groups.
    right_margin : float (default: None)
         Adjust how far the plotting panel extends to the right.
    size : float (default: 3)
         Point size.
    """
    from .. import plotting as plott
    plott.plot_tool(adata,
                    basis='pca',
                    toolkey='pca',
                    smp=smp,
                    names=names,
                    comps=comps,
                    cont=cont,
                    layout=layout,
                    legendloc=legendloc,
                    cmap=cmap,
                    pal=pal,
                    right_margin=right_margin,
                    size=size)

