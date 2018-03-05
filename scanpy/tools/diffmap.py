from ..tools import dpt
from .. import settings
from .. import logging as logg


def diffmap(adata, n_comps=15, n_neighbors=None, knn=True, n_pcs=50, sigma=0,
            n_jobs=None, flavor='haghverdi16', copy=False):
    """Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf17]_.

    Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell
    data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested
    by [Haghverdi16]_ in the implementation of [Wolf17]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    n_comps : `int`, optional (default: 15)
        The number of dimensions of the representation.
    n_jobs : `int` or `None`
        Number of CPUs to use (default: `sc.settings.n_jobs`).
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_diffmap : `adata.obsm`
        Array of shape (#samples) × (#eigen vectors). Diffusion map
        representation of data, which is the right eigen basis of the transition
        matrix with eigenvectors as columns.
    diffmap_evals : `np.ndarray` (`adata.uns`)
        Array of size (number of eigen vectors). Eigenvalues of transition matrix.
    """
    logg.info('computing Diffusion Maps', r=True)
    if n_comps <= 2:
        raise ValueError('Provide any value greater than 2. '
                         'Mind that the 0th component of diffusion maps '
                         'is not used for visualization.')
    adata = adata.copy() if copy else adata
    dmap = dpt.DPT(adata, n_jobs=n_jobs)
    # check whether we need to recompute
    if dmap.evals is None or n_dcs is None or dmap.evals.size < n_dcs:
        dmap.compute_eigen(n_comps=n_comps)
    adata.obsm['X_diffmap'] = dmap.rbasis
    adata.uns['diffmap_evals'] = dmap.evals
    adata.uns['neighbors_distances'] = dmap.distances
    adata.uns['neighbors_similarities'] = dmap.similarities
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'X_diffmap\', tSNE coordinates (adata.obsm)\n'
              '    \'diffmap_evals\', eigenvalues of transition matrix (adata.uns)')
    return adata if copy else None
