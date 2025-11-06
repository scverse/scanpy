from __future__ import annotations

from typing import TYPE_CHECKING

from .._compat import old_positionals
from ._dpt import _diffmap

if TYPE_CHECKING:
    from anndata import AnnData

    from .._utils.random import _LegacyRandom


@old_positionals("neighbors_key", "random_state", "copy")
def diffmap(
    adata: AnnData,
    n_comps: int = 15,
    *,
    neighbors_key: str | None = None,
    random_state: _LegacyRandom = 0,
    copy: bool = False,
) -> AnnData | None:
    """Diffusion Maps :cite:p:`Coifman2005,Haghverdi2015,Wolf2018`.

    Diffusion maps :cite:p:`Coifman2005` have been proposed for visualizing single-cell
    data by :cite:t:`Haghverdi2015`. This tool uses the adapted Gaussian kernel suggested
    by :cite:t:`Haghverdi2016` with the implementation of :cite:t:`Wolf2018`.

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
        If not specified, diffmap looks in .uns['neighbors'] for neighbors settings
        and .obsp['connectivities'] and .obsp['distances'] for connectivities and
        distances, respectively (default storage places for pp.neighbors).
        If specified, diffmap looks in .uns[neighbors_key] for neighbors settings and
        .obsp[.uns[neighbors_key]['connectivities_key']] and
        .obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances,
        respectively.
    random_state
        A numpy random seed
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsm['X_diffmap']` : :class:`numpy.ndarray` (dtype `float`)
        Diffusion map representation of data, which is the right eigen basis of
        the transition matrix with eigenvectors as columns.

    `adata.uns['diffmap_evals']` : :class:`numpy.ndarray` (dtype `float`)
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
        neighbors_key = "neighbors"

    if neighbors_key not in adata.uns:
        msg = "You need to run `pp.neighbors` first to compute a neighborhood graph."
        raise ValueError(msg)
    if n_comps <= 2:
        msg = "Provide any value greater than 2 for `n_comps`. "
        raise ValueError(msg)
    adata = adata.copy() if copy else adata
    _diffmap(
        adata, n_comps=n_comps, neighbors_key=neighbors_key, random_state=random_state
    )
    return adata if copy else None
