from anndata import AnnData
from typing import Optional

from ._dpt import _diffmap
from .._utils import AnyRandom


def diffmap(
    adata: AnnData,
    n_comps: int = 15,
    neighbors_key: Optional[str] = None,
    random_state: AnyRandom = 0,
    copy: bool = False,
):
    """\
    Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf18]_.

    Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell
    data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested
    by [Haghverdi16]_ in the implementation of [Wolf18]_.

    The width ("sigma") of the connectivity kernel is implicitly determined by
    the number of neighbors used to compute the single-cell graph in
    :func:`~scanpy.pp.neighbors`. To reproduce the original implementation
    using a Gaussian kernel, use `method=='gauss'` in
    :func:`~scanpy.pp.neighbors`. To use an exponential kernel, use the default
    `method=='umap'`. Differences between these options shouldn't usually be
    dramatic.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_comps
        The number of dimensions of the representation.
    neighbors_key
        If not specified, diffmap looks .uns['neighbors'] for neighbors settings
        and .obsp['connectivities'], .obsp['distances'] for connectivities and
        distances respectively (default storage places for pp.neighbors).
        If specified, diffmap looks .uns[neighbors_key] for neighbors settings and
        .obsp[.uns[neighbors_key]['connectivities_key']],
        .obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances
        respectively.
    random_state
        A numpy random seed
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    `X_diffmap` : :class:`numpy.ndarray` (`adata.obsm`)
        Diffusion map representation of data, which is the right eigen basis of
        the transition matrix with eigenvectors as columns.
    `diffmap_evals` : :class:`numpy.ndarray` (`adata.uns`)
        Array of size (number of eigen vectors).
        Eigenvalues of transition matrix.

    Notes
    -----
    The 0-th column in `adata.obsm["X_diffmap"]` is the steady-state solution,
    which is non-informative in diffusion maps.
    Therefore, the first diffusion component is at index 1,
    e.g. `adata.obsm["X_diffmap"][:,1]`
    """
    if neighbors_key is None:
        neighbors_key = 'neighbors'

    if neighbors_key not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.'
        )
    if n_comps <= 2:
        raise ValueError('Provide any value greater than 2 for `n_comps`. ')
    adata = adata.copy() if copy else adata
    _diffmap(
        adata, n_comps=n_comps, neighbors_key=neighbors_key, random_state=random_state
    )
    return adata if copy else None
