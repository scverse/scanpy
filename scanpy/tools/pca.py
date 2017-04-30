# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Principal Component Analysis

Notes
-----
There are two PCA versions, which are automatically chosen
- sklearn.decomposition.PCA
- function _pca_fallback
"""
from .. import settings as sett
from .. import preprocess as pp

pca = pp.pca

def plot_pca(adata,
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
             titles=None):
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
    titles : str, optional (default: None)
         Provide titles for panels as "my title1,another title,...".
    """
    from ..examples import check_adata
    adata = check_adata(adata)
    from .. import plotting as plott
    smps = plott.scatter(adata,
                         basis='pca',
                         smp=smp,
                         names=names,
                         comps=comps,
                         cont=cont,
                         layout=layout,
                         legendloc=legendloc,
                         cmap=cmap,
                         pal=pal,
                         right_margin=right_margin,
                         size=size,
                         titles=titles)
    writekey = sett.basekey + '_pca' + sett.plotsuffix
    plott.savefig(writekey)
    if not sett.savefigs and sett.autoshow:
        from ..compat.matplotlib import pyplot as pl
        pl.show()
