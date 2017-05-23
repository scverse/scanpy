# Author: F. Alex Wolf (http://falexwolf.de)
"""Diffusion Map Embedding

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

from ..tools import dpt


def diffmap(adata, n_comps=10, k=30, knn=True, n_pcs=50, sigma=0, n_jobs=None,
            flavor='haghverdi16', copy=False):
    """Diffusion Map Embedding

    Also implements the modifications to diffusion map introduced by Haghverdi
    et al. (2016).

    Return dictionary that stores the new data representation 'Y', which
    consists of the first few eigenvectors of a kernel matrix of the data, and
    the eigenvalues 'evals'.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_comps : int, optional (default: 3)
        The number of dimensions of the representation.
    k : int, optional (default: 30)
        Specify the number of nearest neighbors in the knn graph. If knn ==
        False, set the Gaussian kernel width to the distance of the kth
        neighbor (method 'local').
    knn : bool, optional (default: True)
        If True, use a hard threshold to restrict the number of neighbors to
        k, that is, consider a knn graph. Otherwise, use a Gaussian Kernel
        to assign low weights to neighbors more distant than the kth nearest
        neighbor.
    sigma : float, optional (default: 0)
        If greater 0, ignore parameter 'k', but directly set a global width
        of the Kernel Gaussian (method 'global').
    n_cpus : int or None
        Number of CPUs to use (default: sett.n_cpus).
    copy : bool (default: False)
        Return a copy instead of writing to adata.

    Notes
    -------
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
    dmap = dpt.DPT(adata, k=k, knn=knn, n_pcs=n_pcs,
                   n_jobs=n_jobs, recompute_diffmap=True, flavor=flavor)
    ddmap = dmap.diffmap()
    adata.smp['X_diffmap'] = ddmap['X_diffmap']
    adata.smp['X_diffmap0'] = dmap.rbasis[:, 0]
    adata.add['diffmap_evals'] = ddmap['evals']
    if knn: adata.add['distance'] = dmap.Dsq
    return adata if copy else None
