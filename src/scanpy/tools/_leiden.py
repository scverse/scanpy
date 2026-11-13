from __future__ import annotations

import importlib.util
import os
from collections.abc import Sequence
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

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
from .._settings import Default, settings
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
    from anndata import AnnData
    from igraph import Graph, VertexClustering

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
def leiden(  # noqa: PLR0913
    adata: AnnData,
    resolution: float = 1,
    *,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    rng: SeedLike | RNGLike | None = None,
    directed: bool | Default | None = Default(
        repr="None (in 2.0: True for 'igraph' and 'leidenalg')"
    ),
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: type[MutableVertexPartition] | None = None,
    flavor: LeidenFlavor | Default = Default(preset=("leiden", "flavor")),
    n_jobs: int | None = None,
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
        `None` means directed for `flavor='leidenalg'` and undirected otherwise.
        `flavor='networkit'` only supports undirected graphs.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
        2 is faster and the default for underlying packages.
    partition_type
        Type of partition to use.
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    flavor
        Which package's implementation to use.
        `'networkit'` runs in parallel (see `n_jobs`).
    n_jobs
        Number of threads for `flavor='networkit'` (`-1` means all cores).
        With multiple threads, results for the same `rng` are similar but not identical,
        so `None` means :attr:`scanpy.settings.n_jobs` if `rng=None` is passed
        and `1` otherwise, for exact reproducibility.
    {use}
    {out}
    {adjacency}
    {neighbors_key}
    {obsp}
    copy
        Whether to copy `adata` or modify it inplace.
    **clustering_args
        Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
        or :meth:`igraph.Graph.community_leiden` from `igraph`.

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
    _utils.ensure_igraph()
    clustering_args = dict(clustering_args)
    if n_jobs is None:
        # unseeded runs needn’t be reproducible, so they can be parallel
        n_jobs = settings.n_jobs if rng is None else 1
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
    elif flavor == "networkit":
        part = _leiden_networkit(
            adjacency,
            rng=rng,
            n_jobs=n_jobs,
            use_weights=use_weights,
            n_iterations=n_iterations,
            resolution=resolution,
        )
    else:
        part = _leiden_igraph(
            adjacency,
            rng=rng,
            directed=directed,
            use_weights=use_weights,
            resolution=resolution,
            clustering_args=clustering_args,
        )
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


def _leiden_networkit(
    adjacency: CSBase,
    *,
    rng: np.random.Generator,
    n_jobs: int,
    use_weights: bool,
    n_iterations: int,
    resolution: float | None,
) -> SimpleNamespace:
    _utils.ensure_networkit()
    import networkit

    # Seeding controls which nodes are tried and in what order, but not
    # which thread wins when two move neighboring nodes simultaneously.
    # So runs are statistically reproducible, but only bit-for-bit identical
    # when single-threaded (a fixed thread count > 1 is not enough).
    seed = int(rng.integers(np.iinfo(np.int64).max))
    networkit.setSeed(seed, useThreadId=True)

    # only undirected for Parallel Leiden
    g = _utils.get_networkit_from_adjacency(adjacency, weighted=use_weights)
    # NetworKit's `iterations` caps Leiden passes (stops early on convergence; default is 3).
    # We map scanpy's `n_iterations=-1` ("run until done") to this default,
    # while any positive value explicitly overrides it.
    iterations = n_iterations if n_iterations > 0 else 3
    gamma = 1.0 if resolution is None else resolution
    # randomization was removed as an option, so it is randomize = True
    algorithm = networkit.community.ParallelLeiden(
        g, iterations=iterations, gamma=gamma
    )
    prev_n_threads = networkit.getCurrentNumberOfThreads()
    networkit.setNumberOfThreads(n_jobs if n_jobs > 0 else os.process_cpu_count() or 1)
    try:
        algorithm.run()
    finally:
        networkit.setNumberOfThreads(prev_n_threads)
    nk_part = algorithm.getPartition()
    # NetworKit's Partition exposes getVector() for the labels and a
    # separate Modularity measure, rather than .membership / .modularity.
    return SimpleNamespace(
        membership=np.asarray(nk_part.getVector()),
        modularity=networkit.community.Modularity().getQuality(nk_part, g),
    )


def _leiden_igraph(
    adjacency: CSBase,
    *,
    rng: np.random.Generator,
    directed: bool,
    use_weights: bool,
    resolution: float | None,
    clustering_args: dict[str, object],
) -> VertexClustering:
    g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
    _maybe_suggest_networkit(g)
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
        weights = clustering_args.get("weights")
        if not (weights is None or isinstance(weights, str | Sequence)):
            msg = f"`weights` must be a str or sequence, got {type(weights).__name__}."
            raise TypeError(msg)
        clustering_args |= _directed_node_weights(g, weights=weights)
    with _set_igraph_rng(rng):
        return g.community_leiden(**clustering_args)


def _directed_node_weights(
    g: Graph, *, weights: str | Sequence[float] | None
) -> dict[str, list[float]]:
    """Get in- and out-degrees as node weights; python-igraph otherwise only counts out-degrees."""
    return dict(
        node_weights=g.strength(mode="out", weights=weights),
        node_in_weights=g.strength(mode="in", weights=weights),
    )


def _validate_flavor(
    flavor: str | None,
    *,
    partition_type: object | None,
    directed: bool | Default | None,
) -> tuple[LeidenFlavor, bool]:
    from scanpy import settings

    if was_default := (flavor is None or isinstance(flavor, Default)):
        flavor = settings.preset.leiden.flavor
    is_directed = _resolve_directed(directed, flavor=flavor)
    match flavor:
        case "igraph":
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using igraph."
                raise ValueError(msg)
        case "networkit":
            if is_directed:
                msg = "Cannot use NetworKit's leiden implementation with a directed graph."
                raise ValueError(msg)
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using networkit."
                raise ValueError(msg)
        case "leidenalg":
            msg = (
                "The `igraph` and `networkit` implementations of leiden clustering are *orders of magnitude faster*. "
                "Set the flavor argument to (and install if needed) 'igraph' or 'networkit' to use them."
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
            msg = f"flavor must be either 'igraph', 'leidenalg', or 'networkit', but {flavor!r} was passed."
            raise ValueError(msg)
    return flavor, is_directed


def _resolve_directed(directed: bool | Default | None, *, flavor: str) -> bool:  # noqa: FBT001
    from scanpy import settings

    resolved = (
        settings.preset.leiden.directed if isinstance(directed, Default) else directed
    )
    if isinstance(resolved, frozenset):
        return flavor in resolved
    return flavor == "leidenalg" if resolved is None else resolved


def _maybe_suggest_networkit(g: Graph) -> None:
    """Encourage users toward the parallel NetworKit backend on large graphs."""
    _networkit_edge_heuristic = 500000
    if g.ecount() < _networkit_edge_heuristic:
        return
    if importlib.util.find_spec("networkit") is None:
        return
    logg.hint(
        f"Graph has {g.ecount():,} edges; NetworKit's parallel Leiden "
        f"(`flavor='networkit'`) is typically several times faster than "
        f"igraph's serial implementation."
    )
