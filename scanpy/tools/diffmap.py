# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Diffusion Maps

Diffusion Maps for analysis of single-cell data.

Reference
---------
- Diffusion Maps: Coifman et al., PNAS 102, 7426 (2005).

See also
--------
- Diffusion Maps applied to single-cell data: Haghverdi et al., Bioinformatics
  31, 2989 (2015).
- Diffusion Maps as a flavour of spectral clustering: von Luxburg,
  arXiv:0711.0189 (2007).
"""

# standard modules
from collections import OrderedDict as odict
# scientific modules
import matplotlib
from ..compat.matplotlib import pyplot as pl
from ..tools import dpt
from .. import utils
from .. import settings as sett
from .. import plotting as plott

def diffmap(adata, nr_comps=10, k=5, knn=False, sigma=0):
    """
    Compute diffusion map embedding as of Coifman et al. (2005).

    Also implements the modifications to diffusion map introduced by Haghverdi
    et al. (2016).

    Return dictionary that stores the new data representation 'Y', which 
    consists of the first few eigenvectors of a kernel matrix of the data, and
    the eigenvalues 'evals'. 

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    nr_comps : int, optional (default: 3)
        The number of dimensions of the representation.
    k : int, optional (default: 5)
        Specify the number of nearest neighbors in the knn graph. If knn ==
        False, set the Gaussian kernel width to the distance of the kth
        neighbor (method 'local').
    knn : bool, optional (default: False)
        If True, use a hard threshold to restrict the number of neighbors to
        k, that is, consider a knn graph. Otherwise, use a Gaussian Kernel
        to assign low weights to neighbors more distant than the kth nearest
        neighbor.
    sigma : float, optional (default: 0)
        If greater 0, ignore parameter 'k', but directly set a global width
        of the Kernel Gaussian (method 'global').

    Returns
    -------
    ddmap : dict containing
        Y : np.ndarray
            Array of shape (number of samples) x (number of eigen
            vectors). DiffMap representation of data, which is the right eigen
            basis of transition matrix with eigenvectors as columns.
        evals : np.ndarray
            Array of size (number of cells). Eigenvalues of transition matrix.
    """
    params = locals(); del params['adata']
    dmap = dpt.DPT(adata, params)
    ddmap = dmap.diffmap()
    ddmap['type'] = 'diffmap'
    # restrict number of components
    ddmap['Y'] = ddmap['Y'][:, :params['nr_comps']]
    return ddmap

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
                    subtitles=['diffusion map'],
                    component_name='DC')

