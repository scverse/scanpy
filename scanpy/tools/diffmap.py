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
         comps='1,2',
         cont=None,
         layout='2d',
         legendloc='lower right',
         cmap='jet',
         adjust_right=0.75):
    """
    Scatter plots.

    Parameters
    ----------
    dplot : dict
        Dict returned by plotting tool.
    adata : AnnData
        Annotated data matrix.
    smp : str, optional (default: first annotation)
        Sample annotation to choose for coloring. String annotation is plotted
        assuming categorical annotation, float and integer annotation is plotted
        assuming continuous annoation. Option 'cont' allows to switch between
        these default choices.
    comps : str, optional (default: '1,2')
         String in the form '1,2,3'.
    cont : bool, None (default: None)
        Switch on continuous layout, switch off categorical layout.
    layout : {'2d', '3d', 'unfolded 3d'}, optional (default: '2d')
         Layout of plot.
    legendloc : see matplotlib.legend, optional (default: 'lower right')
         Options for keyword argument 'loc'.
    cmap : str (default: continuous: inferno/ categorical: finite palette)
         String denoting matplotlib color map.
    adjust_right : float (default: 0.75)
         Adjust how far the plotting panel extends to the right.
    """
    from .. import plotting as plott
    plott.plot_tool(dplot, adata,
                    smp,
                    comps,
                    cont,
                    layout,
                    legendloc,
                    cmap,
                    adjust_right,
                    # defined in plotting
                    subtitles=['diffusion map'],
                    component_name='DC')

