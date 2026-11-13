from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scverse_misc import deprecated_arg

from .. import logging
from .._docs import DEPR_NEIGHBORS_KEY, doc_out_key, doc_rng
from .._keys import _embedding_keys
from .._settings import Default, settings
from .._utils import _doc_params
from .._utils.random import _accepts_legacy_random_state
from ..get._out import _out_key, deprecated_copy
from ..neighbors._common import _neighbors_view
from ._docs import doc_use_neighbors

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from anndata.acc import MultiAcc

    from .._utils.random import RNGLike, SeedLike
    from ..get._types import GraphAcc


__all__ = ["diffmap"]


@_doc_params(
    rng=doc_rng,
    use=doc_use_neighbors(meta=True),
    out=doc_out_key(
        "the embedding",
        example="A.obsm['diffmap']",
        default=(
            "`.obsm['X_diffmap']` and the eigenvalues to `.uns['diffmap_evals']` "
            "(`.obsm['diffmap']` and `.uns['diffmap']['evals']` if :attr:`scanpy.settings.preset` is "
            ":attr:`~scanpy.Preset.ScanpyV2Preview`)"
        ),
    ),
)
@_accepts_legacy_random_state(0)
@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_copy
def diffmap(
    adata: AnnData,
    n_comps: int = 15,
    *,
    rng: SeedLike | RNGLike | None = None,
    use: GraphAcc | str | None = None,
    out: MultiAcc | str | bool = False,
    # deprecated
    neighbors_key: str | None = None,
    copy: bool = False,
) -> AnnData | np.ndarray | None:
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

    .. array-support:: tl.diffmap

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_comps
        The number of dimensions of the representation.
    {rng}
    {use}
    {out}
    neighbors_key
        If not specified, diffmap looks in .uns['neighbors'] for neighbors settings
        and .obsp['connectivities'] and .obsp['distances'] for connectivities and
        distances, respectively (default storage places for pp.neighbors).
        If specified, diffmap looks in .uns[neighbors_key] for neighbors settings and
        .obsp[.uns[neighbors_key]['connectivities_key']] and
        .obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances,
        respectively.
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns the embedding if `out=True`, else `None`.
    By default sets the following fields (see `out` for others):

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
    rng = np.random.default_rng(rng)
    if _neighbors_view(adata, use, neighbors_key) is None:
        msg = "You need to run `pp.neighbors` first to compute a neighborhood graph."
        raise ValueError(msg)
    if n_comps <= 2:
        msg = "Provide any value greater than 2 for `n_comps`. "
        raise ValueError(msg)
    key = _out_key(out, None, kind="obsm", default=Default())
    emb = _diffmap(
        adata,
        n_comps=n_comps,
        use=use,
        neighbors_key=neighbors_key,
        key_added=key,
        rng=rng,
    )
    return emb if key is True else None


def _diffmap(
    adata: AnnData,
    n_comps: int = 15,
    *,
    use: GraphAcc | str | None = None,
    neighbors_key: str | None = None,
    key_added: str | Default | Literal[True] | None,
    rng: np.random.Generator,
) -> np.ndarray | None:
    """Compute diffusion maps, returning them for `key_added=True` instead of writing them."""
    from ._dpt import DPT

    start = logging.info(f"computing Diffusion Maps using {n_comps=}(=n_dcs)")
    dpt = DPT(adata, use=use, neighbors_key=neighbors_key)
    dpt.compute_transitions()
    dpt.compute_eigen(n_comps=n_comps, rng=rng)
    if key_added is True:
        logging.info("    finished", time=start)
        return dpt.eigen_basis
    if isinstance(key_added, Default):
        key_added = settings.preset.diffmap.key_added
    keys = _embedding_keys("diffmap", key_added)
    adata.obsm[keys.obsm] = dpt.eigen_basis
    adata.uns[keys.uns], acc = (
        (dpt.eigen_values, "")
        if key_added is None
        else (dict(evals=dpt.eigen_values), "['evals']")
    )
    logging.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {keys.obsm!r}, diffmap coordinates (adata.obsm)\n"
            f"    {keys.uns!r}{acc}, eigenvalues of transition matrix (adata.uns)"
        ),
    )
    return None
