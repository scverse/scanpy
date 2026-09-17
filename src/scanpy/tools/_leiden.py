from __future__ import annotations

import importlib.util
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from natsort import natsorted

from .. import _utils
from .. import logging as logg
from .._compat import warn
from .._docs import doc_rng
from .._settings import Default
from .._utils import _doc_params
from .._utils.random import _accepts_legacy_random_state, _LegacyRng, _set_igraph_rng
from ._docs import (
    doc_adata,
    doc_adjacency,
    doc_neighbors_key,
    doc_obsp,
    doc_restrict_to,
)
from ._utils_clustering import rename_groups, restrict_adjacency

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData
    from igraph import Graph

    from .._compat import CSBase
    from .._settings.presets import LeidenFlavor
    from .._utils.random import RNGLike, SeedLike

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
    rng=doc_rng,
)
@_accepts_legacy_random_state(0)
def leiden(  # noqa: PLR0913, PLR0915
    adata: AnnData,
    resolution: float = 1,
    *,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    rng: SeedLike | RNGLike | None = None,
    key_added: str = "leiden",
    adjacency: CSBase | None = None,
    directed: bool | None = None,
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: type[MutableVertexPartition] | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    copy: bool = False,
    flavor: LeidenFlavor | Default = Default(preset=("leiden", "flavor")),
    **clustering_args,
) -> AnnData | None:
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
    key_added
        `adata.obs` key under which to add the cluster labels.
    {adjacency}
    directed
        Whether to treat the graph as directed or undirected.
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
    {neighbors_key}
    {obsp}
    copy
        Whether to copy `adata` or modify it inplace.
    flavor
        Which package's implementation to use.
    **clustering_args
        Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
        or :meth:`igraph.Graph.community_leiden` from `igraph`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['leiden' | key_added]['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `n_iterations`,
        and `random_state` (if applicable).

    `adata.uns['leiden' | key_added]['modularity']` : :class:`float`
        The modularity score of the final clustering,
        as calculated by the `flavor`.
        Use :func:`scanpy.metrics.modularity`\ `(adata, mode='calculate' | 'update')`
        to calculate a score independent of `flavor`.

    """
    flavor = _validate_flavor(flavor, partition_type=partition_type, directed=directed)
    _utils.ensure_igraph()
    clustering_args = dict(clustering_args)
    rng = np.random.default_rng(rng)
    meta_random_state = (
        dict(random_state=rng.arg) if isinstance(rng, _LegacyRng) else {}
    )

    start = logg.info("running Leiden clustering")
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        adjacency = _utils._choose_graph(adata, obsp, neighbors_key)
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
        directed = True if directed is None else directed
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
        from types import SimpleNamespace

        _utils.ensure_networkit()
        import networkit

        # Seeding controls which nodes are tried and in what order, but not
        # which thread wins when two move neighboring nodes simultaneously.
        # So runs are statistically reproducible for a fixed thread count,
        # but not bit-for-bit identical.
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
        # applying algorithm to the graph
        algorithm.run()
        nk_part = algorithm.getPartition()
        # NetworKit's Partition exposes getVector() for the labels and a
        # separate Modularity measure, rather than .membership / .modularity.

        part = SimpleNamespace(
            # get the actual vector representing the partition data structure
            membership=np.asarray(nk_part.getVector()),
            modularity=networkit.community.Modularity().getQuality(nk_part, g),
        )

    else:
        g = _utils.get_igraph_from_adjacency(adjacency, directed=False)
        _maybe_suggest_networkit(g)
        if use_weights:
            clustering_args["weights"] = "weight"
        if resolution is not None:
            clustering_args["resolution"] = resolution
        clustering_args.setdefault("objective_function", "modularity")
        with _set_igraph_rng(rng):
            part = g.community_leiden(**clustering_args)
    # store output into adata.obs
    groups = np.array(part.membership)
    if restrict_to is not None:
        if key_added == "leiden":
            key_added += "_R"
        groups = rename_groups(
            adata,
            key_added=key_added,
            restrict_key=restrict_key,
            restrict_categories=restrict_categories,
            restrict_indices=restrict_indices,
            groups=groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    # store information on the clustering parameters
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = dict(
        resolution=resolution, n_iterations=n_iterations, **meta_random_state
    )
    adata.uns[key_added]["modularity"] = part.modularity
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key_added!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return adata if copy else None


def _validate_flavor(
    flavor: str | None, *, partition_type: object | None, directed: bool | None
) -> Literal["igraph", "leidenalg", "networkit"]:
    if was_default := (flavor is None or isinstance(flavor, Default)):
        from scanpy import settings

        flavor = settings.preset.leiden.flavor
    match flavor:
        case "igraph":
            if directed:
                msg = "Cannot use igraph’s leiden implementation with a directed graph."
                raise ValueError(msg)
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using igraph."
                raise ValueError(msg)
        case "networkit":
            if directed:
                msg = "Cannot use NetworKit's leiden implementation with a directed graph."
                raise ValueError(msg)
            if partition_type is not None:
                msg = "Do not pass in partition_type argument when using networkit."
                raise ValueError(msg)
        case "leidenalg":
            msg = (
                "The `igraph` implementation of leiden clustering is *orders of magnitude faster*. "
                "Set the flavor argument to (and install if needed) 'igraph' to use it."
            )
            if was_default:
                msg += (
                    "\nIn the future, the default backend for leiden will be igraph instead of leidenalg. "
                    "To achieve the future defaults please pass: `flavor='igraph'` and `n_iterations=2`. "
                    "`directed` must also be `False` to work with igraph’s implementation."
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
    return flavor


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
