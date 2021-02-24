from typing import Optional, Union

import numpy as np
import openTSNE
from anndata import AnnData

from .. import logging as logg
from .._settings import settings
from .._utils import AnyRandom, NeighborsView
from ..tools._utils import get_init_pos_from_paga


class UniformAffinities(openTSNE.affinity.Affinities):
    def __init__(self, neighbors, symmetrize=True, verbose=False):
        self.verbose = verbose

        # Ignore the weights and just assign every neighbor equal weight
        P = neighbors > 0
        # Symmetrize the probability matrix
        if symmetrize:
            P = (P + P.T) / 2
        # Convert weights to probabilities
        P /= np.sum(P)

        self.P = P


def tsne(
    adata: AnnData,
    learning_rate: Union[str, float] = "auto",
    early_exaggeration_iter: int = 250,
    early_exaggeration: Union[float, int] = 12,
    n_iter: int = 500,
    exaggeration: Optional[float] = None,
    dof: float = 1,
    init_pos: Union[str, np.ndarray] = "X_pca",
    initial_momentum: float = 0.5,
    final_momentum: float = 0.8,
    max_step_norm: float = 5,
    random_state: AnyRandom = 0,
    n_jobs: Optional[int] = None,
    neighbors_key: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    t-SNE [Maaten08]_ [Amir13]_ [Policar19]_.

    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been
    proposed for visualizating single-cell data by [Amir13]_. We use the
    implementation of `openTSNE <https://github.com/pavlin-policar/openTSNE>`__
    [Policar19]_.

    Parameters
    ----------
    adata
        Annotated data matrix.
    early_exaggeration
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state
        Change this to use different intial states for the optimization.
        If `None`, the initial state is not reproducible.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_tsne** : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    """
    adata = adata.copy() if copy else adata

    if neighbors_key is None:
        neighbors_key = "neighbors"

    if neighbors_key not in adata.uns:
        raise ValueError(
            f"Did not find .uns['{neighbors_key}']. Run `sc.pp.neighbors` first."
        )

    start = logg.info("computing tSNE")

    # We can be smart about which approximation algorithm to use. BH is faster
    # on tiny data sets
    if adata.shape[0] < 1000:
        negative_gradient_method = "bh"
    else:
        negative_gradient_method = "fft"

    tsne_params = dict(
        learning_rate=learning_rate,
        early_exaggeration_iter=early_exaggeration_iter,
        early_exaggeration=early_exaggeration,
        n_iter=n_iter,
        exaggeration=exaggeration,
        dof=dof,
        initial_momentum=initial_momentum,
        final_momentum=final_momentum,
        max_step_norm=max_step_norm,
        negative_gradient_method=negative_gradient_method,
    )
    adata.uns["tsne"] = {"params": dict(tsne_params)}

    verbose = settings.verbosity > 3

    # Construct affinity matrix using precomputed distance matrix
    neighbors = NeighborsView(adata, neighbors_key)
    affinities = UniformAffinities(neighbors["connectivities"])

    # Get initialization
    if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif isinstance(init_pos, str) and init_pos == "paga":
        init_coords = get_init_pos_from_paga(
            adata, random_state=random_state, neighbors_key=neighbors_key
        )
    elif isinstance(init_pos, str) and init_pos == "spectral":
        init_coords = openTSNE.initialization.spectral(
            affinities.P,
            n_components=2,
            random_state=random_state,
            verbose=verbose,
        )
    else:
        raise ValueError(f"`{init_pos}` is not a valid initialization scheme!")

    # t-SNE requires the initialization to be appropriately rescaled
    init_pos = openTSNE.initialization.rescale(init_coords[:, :2])

    if random_state != 0:
        adata.uns["tsne"]["params"]["random_state"] = random_state

    n_jobs = settings.n_jobs if n_jobs is None else n_jobs

    tsne = openTSNE.TSNE(
        affinities=affinities,
        initialization=init_pos,
        random_state=random_state,
        n_jobs=n_jobs,
        verbose=verbose,
        **tsne_params,
    )
    # We have to call fit with a matrix, but it's arbitrary, since we already
    # have all the information we need in affinitis and the initalization,
    # so we'll arbitrarily call it with the initialization
    embedding = tsne.fit(init_pos)

    # update AnnData instance
    adata.obsm["X_tsne"] = embedding.view(np.ndarray)
    logg.info(
        "    finished",
        time=start,
        deep="added\n    'X_tsne', tSNE coordinates (adata.obsm)",
    )
    return adata if copy else None
