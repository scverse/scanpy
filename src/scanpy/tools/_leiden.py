from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, cast

import numpy as np
import pandas as pd
from natsort import natsorted
from scverse_misc import deprecated_arg

from .. import _utils
from .. import logging as logg
from .._compat import warn
from .._docs import (
    DEPR_ADJACENCY,
    DEPR_KEY_ADDED,
    DEPR_NEIGHBORS_KEY,
    DEPR_OBSP_GRAPH,
    doc_out_key,
    doc_rng,
)
from .._settings import Default
from .._settings.presets import LeidenPreset
from .._utils import _doc_params
from .._utils.random import _accepts_legacy_random_state, _LegacyRng, _set_igraph_rng
from ..get._out import _out_key, deprecated_copy
from ..neighbors._common import _use_graph
from ._docs import (
    doc_adata,
    doc_adjacency,
    doc_neighbors_key,
    doc_obsp,
    doc_restrict_to,
    doc_use_neighbors,
)
from ._utils_clustering import rename_groups, restrict_adjacency

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData
    from igraph import Graph
    from numpy.typing import NDArray

    from .._compat import CSBase
    from .._settings.presets import LeidenFlavor
    from .._utils.random import RNGLike, SeedLike
    from ..get._types import GraphAcc, VecRef

    try:  # sphinx-autodoc-typehints + optional dependency
        from leidenalg.VertexPartition import MutableVertexPartition
    except ImportError:
        if not TYPE_CHECKING:
            MutableVertexPartition = type(
                "MutableVertexPartition",
                (),
                dict(__module__="leidenalg.VertexPartition"),
            )


@_doc_params(
    doc_adata=doc_adata,
    restrict_to=doc_restrict_to,
    adjacency=doc_adjacency,
    neighbors_key=doc_neighbors_key.format(method="leiden"),
    obsp=doc_obsp,
    use=doc_use_neighbors(),
    rng=doc_rng,
    out=doc_out_key(
        "the cluster labels",
        example="A.obs['clusters']",
        key_added="`key_added='k'` is `out=A.obs['k']`.",
        default=(
            "`.obs['leiden']` and `.uns['leiden']` "
            "(`'leiden_R'` if `restrict_to` is passed)"
        ),
    ),
)
@_accepts_legacy_random_state(0)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("adjacency", DEPR_ADJACENCY)
@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_arg("obsp", DEPR_OBSP_GRAPH)
@deprecated_copy
def leiden(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    resolution: float = 1,
    *,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    rng: SeedLike | RNGLike | None = None,
    directed: bool | Default | None = Default(preset=("leiden", "directed")),
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: type[MutableVertexPartition] | None = None,
    flavor: LeidenFlavor | Default = Default(preset=("leiden", "flavor")),
    use: GraphAcc | str | None = None,
    out: VecRef | bool = False,
    # deprecated
    adjacency: CSBase | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
    **clustering_args,
) -> AnnData | pd.Categorical | None:
    r"""Cluster cells into subgroups :cite:p:`Traag2019`.

    Cluster cells using the Leiden algorithm :cite:p:`Traag2019`,
    an improved version of the Louvain algorithm :cite:p:`Blondel2008`.
    It was proposed for single-cell analysis by :cite:t:`Levine2015`.

    This requires having run :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    .. array-support:: tl.leiden

    Parameters
    ----------
    {doc_adata}
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.
        Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    {rng}
    {restrict_to}
    directed
        Whether to treat the graph as directed or undirected.
        `None` means directed for `flavor='leidenalg'` and undirected for `flavor='igraph'`.
        `flavor='graspologic'` only supports undirected graphs.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
        2 is faster and the default for underlying packages.
        `flavor='graspologic'` can’t run until convergence, so negative values mean 2.
    partition_type
        Type of partition to use.
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    flavor
        Which package's implementation to use.
    {use}
    {out}
    {adjacency}
    {neighbors_key}
    {obsp}
    copy
        Whether to copy `adata` or modify it inplace.
    **clustering_args
        Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
        or :meth:`igraph.Graph.community_leiden` from `igraph`,
        or ``graspologic_native.leiden_csr``.

    Returns
    -------
    Returns the cluster labels if `out=True`, else `None`.
    By default sets the following fields (see `out` for others):

    `adata.obs['leiden']` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['leiden']['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `n_iterations`,
        and `random_state` (if applicable).

    `adata.uns['leiden']['modularity']` : :class:`float`
        The modularity score of the final clustering,
        as calculated by the `flavor`.
        Use :func:`scanpy.metrics.modularity`\ `(adata, mode='calculate' | 'update')`
        to calculate a score independent of `flavor`.

    """
    flavor, directed = _validate_flavor(
        flavor, partition_type=partition_type, directed=directed
    )
    if flavor != "graspologic":
        _utils.ensure_igraph()
    clustering_args = dict(clustering_args)
    rng = np.random.default_rng(rng)
    meta_random_state = (
        dict(random_state=rng.arg) if isinstance(rng, _LegacyRng) else {}
    )

    key = _out_key(out, key_added, kind="obs", default="leiden")
    # like `key_added`, the default key gets a suffix, keys from `out` are literal
    if restrict_to is not None and key == "leiden" and out is False:
        key += "_R"

    start = logg.info("running Leiden clustering")
    # are we clustering a user-provided graph or the default AnnData one?
    adjacency = _use_graph(
        adata, use, adjacency=adjacency, obsp=obsp, neighbors_key=neighbors_key
    )
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories=restrict_categories,
            adjacency=adjacency,
        )
    # Prepare find_partition arguments as a dictionary,
    # appending to whatever the user provided. It needs to be this way
    # as this allows for the accounting of a None resolution
    # (in the case of a partition variant that doesn't take it on input)
    clustering_args["n_iterations"] = n_iterations
    if flavor == "leidenalg":
        import leidenalg

        if resolution is not None:
            clustering_args["resolution_parameter"] = resolution
        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if partition_type is None:
            partition_type = leidenalg.RBConfigurationVertexPartition
        if use_weights:
            clustering_args["weights"] = np.array(g.es["weight"]).astype(np.float64)
        seed = (
            rng.arg
            if isinstance(rng, _LegacyRng) and isinstance(rng.arg, int | np.integer)
            # for some reason leidenalg only accepts int32 (signed) seeds …
            else rng.integers((i := np.iinfo(np.int32)).min, i.max, dtype=np.int32)
        )
        part = cast(
            "MutableVertexPartition",
            leidenalg.find_partition(g, partition_type, seed=seed, **clustering_args),
        )
    elif flavor == "graspologic":
        part = _leiden_graspologic(
            adjacency,
            resolution=resolution,
            use_weights=use_weights,
            rng=rng,
            **clustering_args,
        )
    else:
        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if use_weights:
            clustering_args["weights"] = "weight"
        if resolution is not None:
            clustering_args["resolution"] = resolution
        clustering_args.setdefault("objective_function", "modularity")
        if (
            directed
            and clustering_args["objective_function"] == "modularity"
            and clustering_args.keys().isdisjoint({"node_weights", "node_in_weights"})
        ):
            clustering_args |= _directed_node_weights(
                g, weights=clustering_args.get("weights")
            )
        with _set_igraph_rng(rng):
            part = g.community_leiden(**clustering_args)
    # store output into adata.obs
    groups = np.array(part.membership)
    if restrict_to is not None:
        groups = rename_groups(
            adata,
            restrict_key=restrict_key,
            restrict_categories=restrict_categories,
            restrict_indices=restrict_indices,
            groups=groups,
        )
    labels = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    if key is True:
        logg.info("    finished", time=start)
        return labels
    adata.obs[key] = labels
    # store information on the clustering parameters
    adata.uns[key] = {}
    adata.uns[key]["params"] = dict(
        resolution=resolution, n_iterations=n_iterations, **meta_random_state
    )
    adata.uns[key]["modularity"] = part.modularity
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return None


def _directed_node_weights(
    g: Graph, *, weights: str | Sequence[float] | None
) -> dict[str, list[float]]:
    """Get in- and out-degrees as node weights; python-igraph otherwise only counts out-degrees."""
    return dict(
        node_weights=g.strength(mode="out", weights=weights),
        node_in_weights=g.strength(mode="in", weights=weights),
    )


class _Partition(NamedTuple):
    membership: NDArray[np.int64]
    modularity: float


def _leiden_graspologic(
    adjacency: CSBase | np.ndarray,
    *,
    resolution: float,
    n_iterations: int,
    use_weights: bool,
    rng: np.random.Generator,
    **kwargs,
) -> _Partition:
    from graspologic_native import leiden_csr
    from scipy import sparse

    a = sparse.csr_array(adjacency)  # noqa: TID251
    # a sorted copy without explicit zeros: graspologic’s result depends on the storage order
    a = a.maximum(a.T).tocoo().tocsr()
    a.eliminate_zeros()
    n = a.shape[0]
    if a.nnz == 0:  # graspologic panics without edges
        return _Partition(np.arange(n), 0.0)
    modularity, membership = leiden_csr(
        a.indptr.astype(np.int64),
        a.indices.astype(np.int32),
        a.data.astype(np.float64) if use_weights else np.ones(a.nnz),
        n,
        resolution=resolution,
        iterations=2 if n_iterations < 0 else n_iterations,
        seed=int(rng.integers(2**63)),
        **kwargs,
    )
    return _Partition(
        np.fromiter(map(membership.__getitem__, range(n)), np.int64, n), modularity
    )


def _validate_flavor(  # noqa: PLR0912
    flavor: str | None,
    *,
    partition_type: object | None,
    directed: bool | Default | None,
) -> LeidenPreset:
    from scanpy import settings

    if was_default := (flavor is None or isinstance(flavor, Default)):
        flavor = settings.preset.leiden.flavor
    if directed_was_default := isinstance(directed, Default):
        directed = settings.preset.leiden.directed
    match flavor:
        case "igraph":
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using igraph."
                raise ValueError(msg)
        case "graspologic":
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using graspologic."
                raise ValueError(msg)
            if directed and not directed_was_default:
                msg = "graspologic’s leiden implementation only supports undirected graphs."
                raise ValueError(msg)
            directed = False
            try:
                import graspologic_native  # noqa: F401
            except ImportError as e:
                e.add_note(
                    "Please install `scanpy[graspologic]` (or `graspologic-native` directly) and try again."
                )
                raise
        case "leidenalg":
            msg = (
                "The `igraph` implementation of leiden clustering is *orders of magnitude faster*. "
                "Set the flavor argument to (and install if needed) 'igraph' to use it."
            )
            if was_default:
                msg += (
                    "\nIn the future, the default backend for leiden will be igraph instead of leidenalg. "
                    "To achieve the future defaults please pass: `flavor='igraph'` and `n_iterations=2`."
                )
            warn(msg, FutureWarning if was_default else UserWarning)
            try:
                import leidenalg  # noqa: F401
            except ImportError as e:
                e.add_note(
                    "Please install `scanpy[leiden]` (or `leidenalg` directly) and try again."
                )
                raise
        case _:
            msg = f"flavor must be 'igraph', 'leidenalg' or 'graspologic', but {flavor!r} was passed."
            raise ValueError(msg)
    if directed is None:
        directed = flavor == "leidenalg"
    return LeidenPreset(flavor, directed)
