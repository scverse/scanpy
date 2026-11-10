from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scverse_misc import deprecated_arg
from sklearn.utils import check_array

from .. import logging as logg
from .._compat import pkg_version
from .._docs import DEPR_KEY_ADDED, DEPR_NEIGHBORS_KEY, doc_out_key, doc_rng
from .._keys import _embedding_keys
from .._settings import Default, settings
from .._utils import _doc_params
from .._utils.random import (
    _accepts_legacy_random_state,
    _legacy_random_state,
    _LegacyRng,
)
from ..get._out import _out_key, deprecated_copy
from ..get._types import MultiAcc
from ..get._uns import _use_from_uns
from ..get._use import _resolve_obs, _resolve_ref
from ..neighbors._common import _neighbors_view
from ._docs import doc_use_neighbors
from ._utils import _choose_representation_compat, get_init_pos_from_paga

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, Literal

    from anndata import AnnData

    from .._compat import CSBase
    from .._utils.random import RNGLike, SeedLike
    from ..get._types import GraphAcc, RepAcc


type _InitPos = Literal["auto", "paga", "spectral", "random", "recursive"]


@_accepts_legacy_random_state(0)
@_doc_params(
    rng=doc_rng,
    use=doc_use_neighbors(
        meta=True,
        extra=(
            "\n    UMAP lays out more than `2 * n_components` connected components using data,"
            " by default the representation :func:`~scanpy.pp.neighbors` used."
            "\n    If there is none, pass `rep`."
        ),
    ),
    out=doc_out_key(
        "the embedding",
        example="A.obsm['umap']",
        key_added="`key_added='k'` is `out=A.obsm['k']`.",
        default=(
            "`.obsm['X_umap']` and `.uns['umap']` "
            "(`.obsm['umap']` if :attr:`scanpy.settings.preset` is "
            ":attr:`~scanpy.Preset.ScanpyV2Preview`)"
        ),
    ),
)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_copy
def umap(  # noqa: PLR0913
    adata: AnnData,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: int | None = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: _InitPos | MultiAcc | str | np.ndarray | None = "auto",
    rng: SeedLike | RNGLike | None = None,
    a: float | None = None,
    b: float | None = None,
    method: Literal["umap"] = "umap",
    compatibility_layout: bool | Default = Default(
        preset=("umap", "compatibility_layout")
    ),
    use: GraphAcc | str | None = None,
    rep: RepAcc | str | None = None,
    out: MultiAcc | str | bool = False,
    # deprecated
    neighbors_key: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | np.ndarray | None:
    r"""Embed the neighborhood graph using UMAP :cite:p:`McInnes2018`.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.
    We use the implementation of umap-learn_ :cite:p:`McInnes2018`.
    For a few comparisons of UMAP with tSNE, see :cite:t:`Becht2018`.

    .. _umap-learn: https://github.com/lmcinnes/umap

    .. array-support:: tl.umap

    Parameters
    ----------
    adata
        Annotated data matrix.
    min_dist
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components
        The number of dimensions of the embedding.
    maxiter
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha
        The initial learning rate for the embedding optimization.
    gamma
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Options are:

        * An `.obsm` entry, e.g. `A.obsm['pca']`.
          A :class:`str` is :meth:`~anndata.acc.AdAcc.resolve`\ d
          if :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`
          (e.g. `'obsm.pca'`), otherwise it is a key of `.obsm`.
        * 'auto': 'spectral' if `compatibility_layout`, else 'recursive'
          (umap-learn’s choice).
        * 'paga': positions from :func:`~scanpy.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * 'recursive': label propagation-based hierarchical initialization
          (requires umap-learn ≥0.6).
        * A numpy array of initial embedding positions.
    {rng}
    a
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    method
        Chosen implementation.

        ``'umap'``
            Umap’s simplical set embedding.
    compatibility_layout
        Use the layout pipeline of umap-learn 0.5 (`True`)
        or the faster one introduced in umap-learn 0.6 (`False`, requires umap-learn ≥0.6).
        Stored in `.uns['umap']['params']`, so :func:`~scanpy.tl.ingest` can match it.
    {use}
    rep
        The data to lay out more than `2 * n_components` connected components with,
        e.g. `A.obsm['pca']`, if not the representation :func:`~scanpy.pp.neighbors` used.
        :class:`str`\ s are :meth:`~anndata.acc.AdAcc.resolve`\ d, e.g. `'obsm.pca'`.
    {out}
    neighbors_key
        Umap looks in
        :attr:`~anndata.AnnData.uns`\ ``[neighbors_key]`` for neighbors settings and
        :attr:`~anndata.AnnData.obsp`\ ``[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns the embedding if `out=True`, else `None`.
    By default sets the following fields (see `out` for others):

    `adata.obsm['X_umap']` : :class:`numpy.ndarray` (dtype `float`)
        UMAP coordinates of data.
    `adata.uns['umap']` : :class:`dict`
        UMAP parameters.

    """
    non_deterministic = rng is None
    rng = np.random.default_rng(rng)
    key = _out_key(out, key_added, kind="obsm", default=Default())

    neighbors = _neighbors_view(adata, use, neighbors_key)
    if neighbors is None or "connectivities" not in neighbors:
        msg = "Did not find a neighborhood graph. Run `sc.pp.neighbors` first."
        raise ValueError(msg)
    conns = neighbors["connectivities"]
    neigh_params = neighbors["params"] if "params" in neighbors else {}  # noqa: SIM401 (no `.get`)

    start = logg.info("computing UMAP")

    if neigh_params and neigh_params["method"] != "umap":
        logg.warning(
            f'.obsp["{neighbors["connectivities_key"]}"] have not been computed using umap'
        )

    from umap.umap_ import find_ab_params, simplicial_set_embedding

    compatibility_layout = _resolve_compatibility_layout(compatibility_layout)
    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    meta_random_state = (
        dict(random_state=rng.arg) if isinstance(rng, _LegacyRng) else {}
    )
    params = dict(
        a=a, b=b, compatibility_layout=compatibility_layout, **meta_random_state
    )
    init_coords = _init_coords(adata, init_pos, rng=rng, conns=conns)
    x = _umap_data(adata, conns, rep, neigh_params, n_components=n_components)
    if method == "umap":
        default_epochs = 500 if conns.shape[0] <= 10000 else 200
        n_epochs = default_epochs if maxiter is None else maxiter
        x_umap, _ = simplicial_set_embedding(
            data=x,
            graph=conns.tocoo(copy=True),
            n_components=n_components,
            initial_alpha=alpha,
            a=a,
            b=b,
            gamma=gamma,
            negative_sample_rate=negative_sample_rate,
            n_epochs=n_epochs,
            init=init_coords,
            random_state=_legacy_random_state(rng, always_state=True),
            parallel=non_deterministic,  # if True, random_state is ignored
            metric=neigh_params.get("metric", "euclidean"),
            metric_kwds=neigh_params.get("metric_kwds", {}),
            densmap=False,
            densmap_kwds={},
            output_dens=False,
            verbose=settings.verbosity > 3,
            **_umap_layout_kw(compatibility_layout=compatibility_layout),
        )
    else:
        msg = f"Unknown method {method}"
        raise ValueError(msg)
    if key is True:
        logg.info("    finished", time=start)
        return x_umap
    keys = _embedding_keys("umap", key)
    adata.uns[keys.uns] = dict(params=params)
    adata.obsm[keys.obsm] = x_umap  # annotate samples with UMAP coordinates
    logg.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {keys.obsm!r}, UMAP coordinates (adata.obsm)\n"
            f"    {keys.uns!r}, UMAP parameters (adata.uns)"
        ),
    )
    return None


def _init_coords(
    adata: AnnData,
    init_pos: _InitPos | MultiAcc | str | np.ndarray | None,
    *,
    rng: np.random.Generator,
    conns: CSBase,
) -> np.ndarray | str | None:
    """Resolve `init_pos` into coordinates or a method name umap-learn understands."""
    if (
        isinstance(init_pos, str)
        and init_pos == "auto"
        and not _umap_layout_kw(compatibility_layout=False)
    ):
        return "spectral"  # umap-learn <0.6 has no 'auto'
    # method names have no dots
    ref = (
        _resolve_ref(init_pos, vec=False)
        if isinstance(init_pos, str) and "." in init_pos
        else init_pos
    )
    match ref:
        case MultiAcc():
            coords = adata[ref]
        case str() if ref in adata.obsm:
            coords = adata.obsm[ref]
        case "paga":
            coords = get_init_pos_from_paga(adata, rng=rng, adjacency=conns)
        case str() | None:
            return ref  # let umap handle it
        case _:
            coords = ref
    return check_array(coords, dtype=np.float32, accept_sparse=False)


def _umap_data(
    adata: AnnData,
    conns: CSBase,
    rep: RepAcc | str | None,
    neigh_params: Mapping[str, Any],
    *,
    n_components: int,
) -> Any:
    """Get the data UMAP uses to lay out the graph’s connected components relative to each other.

    UMAP only reads it if there are more than `2 * n_components` components,
    and guesses the layout without a seed if it’s missing, so we require it then.
    """
    if rep is not None:
        return adata[_resolve_obs(rep, rep=True)]
    if neigh_params:  # the representation `pp.neighbors` used
        return _choose_representation_compat(
            adata,
            use_rep=_use_from_uns(neigh_params),
            n_pcs=neigh_params.get("n_pcs", None),
            silent=True,
        )
    from scipy.sparse.csgraph import connected_components

    if (n := connected_components(conns)[0]) > 2 * n_components:
        msg = (
            f"The graph has {n} connected components, so UMAP needs data to lay them out. "
            "Pass e.g. `rep=A.obsm['pca']`."
        )
        raise ValueError(msg)
    return None


def _umap_layout_kw(*, compatibility_layout: bool) -> dict[str, bool]:
    """Get `compatibility_layout` kwarg, empty for umap-learn <0.6, which lacks the choice."""
    if pkg_version("umap-learn").release < (0, 6):
        return {}
    return dict(compatibility_layout=compatibility_layout)


def _resolve_compatibility_layout(compatibility_layout: bool | Default) -> bool:  # noqa: FBT001
    if not isinstance(compatibility_layout, Default):
        if not (compatibility_layout or _umap_layout_kw(compatibility_layout=False)):
            msg = "`compatibility_layout=False` requires umap-learn ≥0.6."
            raise ImportError(msg)
        return compatibility_layout
    # umap-learn <0.6 only has the compatibility layout
    return settings.preset.umap.compatibility_layout or not _umap_layout_kw(
        compatibility_layout=False
    )
