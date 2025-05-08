from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np

from .. import logging as logg
from .._settings import settings
from .._utils import _choose_graph

if TYPE_CHECKING:
    from anndata import AnnData

    from .._compat import CSRBase, SpBase


def _choose_representation(
    adata: AnnData,
    *,
    use_rep: str | None = None,
    n_pcs: int | None = None,
    silent: bool = False,
) -> np.ndarray | CSRBase:  # TODO: what else?
    verbosity = settings.verbosity
    if silent and settings.verbosity > 1:
        settings.verbosity = 1
    if use_rep is None and n_pcs == 0:  # backwards compat for specifying `.X`
        use_rep = "X"
    if use_rep is None:
        X = _get_pca_or_small_x(adata, n_pcs)
    elif use_rep in adata.obsm and n_pcs is not None:
        if n_pcs > adata.obsm[use_rep].shape[1]:
            msg = (
                f"{use_rep} does not have enough Dimensions. Provide a "
                "Representation with equal or more dimensions than"
                "`n_pcs` or lower `n_pcs` "
            )
            raise ValueError(msg)
        X = adata.obsm[use_rep][:, :n_pcs]
    elif use_rep in adata.obsm and n_pcs is None:
        X = adata.obsm[use_rep]
    elif use_rep == "X":
        X = adata.X
    else:
        msg = f"Did not find {use_rep} in `.obsm.keys()`. You need to compute it first."
        raise ValueError(msg)
    settings.verbosity = verbosity  # resetting verbosity
    return X


def _get_pca_or_small_x(adata: AnnData, n_pcs: int | None) -> np.ndarray | CSRBase:
    if adata.n_vars <= settings.N_PCS:
        logg.info("    using data matrix X directly")
        return adata.X

    if "X_pca" in adata.obsm:
        if n_pcs is not None and n_pcs > adata.obsm["X_pca"].shape[1]:
            msg = "`X_pca` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`."
            raise ValueError(msg)
        X = adata.obsm["X_pca"][:, :n_pcs]
        logg.info(f"    using 'X_pca' with n_pcs = {X.shape[1]}")
        return X

    from ..preprocessing import pca

    warnings.warn(
        f"You’re trying to run this on {adata.n_vars} dimensions of `.X`, "
        "if you really want this, set `use_rep='X'`.\n         "
        "Falling back to preprocessing with `sc.pp.pca` and default params.",
        stacklevel=3,
    )
    n_pcs_pca = n_pcs if n_pcs is not None else settings.N_PCS
    pca(adata, n_comps=n_pcs_pca)
    return adata.obsm["X_pca"]


def get_init_pos_from_paga(
    adata: AnnData,
    adjacency: SpBase | None = None,
    random_state=0,
    neighbors_key: str | None = None,
    obsp: str | None = None,
):
    np.random.seed(random_state)
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    if "pos" not in adata.uns.get("paga", {}):
        msg = "Plot PAGA first, so that `adata.uns['paga']['pos']` exists."
        raise ValueError(msg)

    groups = adata.obs[adata.uns["paga"]["groups"]]
    pos = adata.uns["paga"]["pos"]
    connectivities_coarse = adata.uns["paga"]["connectivities"]
    init_pos = np.ones((adjacency.shape[0], 2))
    for i, group_pos in enumerate(pos):
        subset = (groups == groups.cat.categories[i]).values
        neighbors = connectivities_coarse[i].nonzero()
        if len(neighbors[1]) > 0:
            connectivities = connectivities_coarse[i][neighbors]
            nearest_neighbor = neighbors[1][np.argmax(connectivities)]
            noise = np.random.random((len(subset[subset]), 2))
            dist = group_pos - pos[nearest_neighbor]
            noise = noise * dist
            init_pos[subset] = group_pos - 0.5 * dist + noise
        else:
            init_pos[subset] = group_pos
    return init_pos
