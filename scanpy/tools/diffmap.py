# Author: F. Alex Wolf (http://falexwolf.de)
"""Diffusion Maps
"""

from ..tools import dpt
from .. import logging as logg

def diffmap(adata, n_comps=15, n_neighbors=30, knn=True, n_pcs=50, sigma=0, n_jobs=None,
            flavor='haghverdi16', copy=False):
    """Diffusion Maps

    Visualize data using Diffusion Maps.

    Implements the crucial modifications to diffusion map introduced by
    Haghverdi et al., Nature Methods (2016).

    References
    ----------
    - Diffusion Maps: Coifman et al., PNAS 102, 7426 (2005).
    - Diffusion Maps applied to single-cell data: Haghverdi et al., Bioinformatics
      31, 2989 (2015).
    - Diffusion Pseudotime: Haghverdi et al., Nature Methods 13, 3971 (2016).

    See also
    --------
    - Diffusion Maps as a flavour of spectral clustering: von Luxburg,
      arXiv:0711.0189 (2007).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_comps : int, optional (default: 10)
        The number of dimensions of the representation.
    n_neighbors : int, optional (default: 30)
        Specify the number of nearest neighbors in the knn graph. If knn ==
        False, set the Gaussian kernel width to the distance of the kth
        neighbor (method 'local').
    knn : bool, optional (default: True)
        If True, use a hard threshold to restrict the number of neighbors to
        k, that is, consider a knn graph. Otherwise, use a Gaussian Kernel
        to assign low weights to neighbors more distant than the kth nearest
        neighbor.
    n_pcs : int, optional (default: 50)
        Use n_pcs PCs to compute the Euclidian distance matrix, which is the
        basis for generating the graph. Set to 0 if you don't want preprocessing
        with PCA.
    n_jobs : int or None
        Number of CPUs to use (default: sett.n_cpus).
    copy : bool (default: False)
        Return a copy instead of writing to adata.

    Notes
    -----
    The following is added to adata.smp
        X_diffmap : np.ndarray
            Array of shape n_samples x n_comps. DiffMap representation of data,
            which is the right eigen basis of transition matrix with eigenvectors as
            columns.
    The following is added to adata.add
        diffmap_evals : np.ndarray
            Eigenvalues of the transition matrix.
    """
    adata = adata.copy() if copy else adata
    dmap = dpt.DPT(adata, n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs,
                   n_dcs=n_comps, n_jobs=n_jobs, recompute_graph=True,
                   flavor=flavor)
    dmap.update_diffmap()
    adata.add['data_graph_distance_local'] = dmap.Dsq
    adata.add['data_graph_norm_weights'] = dmap.Ktilde
    adata.smp['X_diffmap'] = dmap.rbasis[:, 1:]
    adata.smp['X_diffmap0'] = dmap.rbasis[:, 0]
    adata.add['diffmap_evals'] = dmap.evals[1:]
    logg.m('    finished', t=True, end=' ')
    logg.m('and added\n'
           '    "X_diffmap", the diffmap coordinates (adata.smp),\n'
           '    "diffmap_evals", the eigenvalues of the transition matrix (adata.add)')
    return adata if copy else None
