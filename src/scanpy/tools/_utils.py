from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .. import logging as logg
from .._compat import warn
from .._keys import _existing_preset_keys
from .._settings import Preset, settings
from .._utils import _choose_graph
from ..get.get import LayerAcc, MultiAcc, _get_arr, _resolve_rep

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import CSBase, CSRBase
    from ..get.get import RepAcc


def _choose_representation(
    adata: AnnData,
    *,
    use_rep: RepAcc | str | None,
    n_pcs: int | None,
    silent: bool = False,
) -> np.ndarray | CSRBase:  # TODO: what else?
    """Get the representation to compute on, resolving strings using `anndata.acc`."""
    return _choose_representation_compat(
        adata,
        use_rep=None if use_rep is None else _resolve_rep(use_rep),
        n_pcs=n_pcs,
        silent=silent,
    )


def _choose_representation_compat(
    adata: AnnData,
    *,
    use_rep: RepAcc | str | None,
    n_pcs: int | None,
    silent: bool = False,
) -> np.ndarray | CSRBase:  # TODO: what else?
    """Get the representation to compute on.

    Treats strings as `.obsm` keys (or `'X'`) instead of `anndata.acc` specs when the preset is v1.
    """
    verbosity = settings.verbosity
    if silent and settings.verbosity > 1:
        settings.verbosity = 1
    if use_rep is not None and (
        not isinstance(use_rep, str) or settings.preset is Preset.ScanpyV2Preview
    ):
        use_rep = _resolve_rep(use_rep)
    if use_rep is None and n_pcs == 0:  # backwards compat for specifying `.X`
        use_rep = "X"
    match use_rep:
        case None:
            x = _get_pca_or_small_x(adata, n_pcs)
        case LayerAcc():
            x = _get_arr(adata, use_rep)
        case MultiAcc():
            x = _slice_n_pcs(_get_arr(adata, use_rep), n_pcs, use_rep)
        case str() if use_rep in adata.obsm:
            x = _slice_n_pcs(adata.obsm[use_rep], n_pcs, use_rep)
        case "X":
            x = adata.X
        case _:
            msg = f"Did not find {use_rep} in `.obsm.keys()`. You need to compute it first."
            raise ValueError(msg)
    settings.verbosity = verbosity  # resetting verbosity
    return x


def _slice_n_pcs[A: np.ndarray | CSBase](
    x: A, n_pcs: int | None, use_rep: RepAcc | str
) -> A:
    if n_pcs is None:
        return x
    if n_pcs > x.shape[1]:
        msg = (
            f"{use_rep} does not have enough Dimensions. Provide a "
            "Representation with equal or more dimensions than"
            "`n_pcs` or lower `n_pcs` "
        )
        raise ValueError(msg)
    return x[:, :n_pcs]


def _get_pca_or_small_x(adata: AnnData, n_pcs: int | None) -> np.ndarray | CSRBase:
    from .._keys import _embedding_keys
    from ..preprocessing._pca import pca

    if adata.n_vars <= settings.N_PCS:
        logg.info("    using data matrix X directly")
        return adata.X

    if keys := _existing_preset_keys(adata, "pca"):
        if n_pcs is not None and n_pcs > adata.obsm[keys.obsm].shape[1]:
            msg = f"`adata.obsm[{keys.obsm!r}]` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`."
            raise ValueError(msg)
        x = adata.obsm[keys.obsm][:, :n_pcs]
        logg.info(f"    using {keys.obsm!r} with n_pcs = {x.shape[1]}")
        return x

    msg = (
        f"You’re trying to run this on {adata.n_vars} dimensions of `.X`, "
        "if you really want this, set `use_rep=’X’`.\n         "
        "Falling back to preprocessing with `sc.pp.pca` and default params."
    )
    warn(msg, UserWarning)
    n_pcs_pca = n_pcs if n_pcs is not None else settings.N_PCS
    pca(adata, n_comps=n_pcs_pca)
    return adata.obsm[_embedding_keys("pca").obsm]


def get_init_pos_from_paga(
    adata: AnnData,
    *,
    rng: np.random.Generator,
    adjacency: CSBase | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
) -> NDArray[np.float64]:
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    if "pos" not in adata.uns.get("paga", {}):
        msg = "Plot PAGA first, so that `adata.uns['paga']['pos']` exists."
        raise ValueError(msg)

    groups = adata.obs[adata.uns["paga"]["groups"]]
    pos = adata.uns["paga"]["pos"]
    connectivities_coarse = adata.uns["paga"]["connectivities"]
    init_pos = np.ones((adjacency.shape[0], 2))
    # spawn sub-rngs so this can maybe be parallelized without changing random number generation
    for i, (group_pos, sub_rng) in enumerate(
        zip(pos, rng.spawn(len(pos)), strict=True)
    ):
        subset = (groups == groups.cat.categories[i]).values
        neighbors = connectivities_coarse[i].nonzero()
        if len(neighbors[1]) > 0:
            connectivities = connectivities_coarse[i][neighbors]
            nearest_neighbor = neighbors[1][np.argmax(connectivities)]
            noise = sub_rng.random((len(subset[subset]), 2))
            dist = group_pos - pos[nearest_neighbor]
            noise = noise * dist
            init_pos[subset] = group_pos - 0.5 * dist + noise
        else:
            init_pos[subset] = group_pos
    return init_pos
