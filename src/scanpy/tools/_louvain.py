from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 15):
    from types import MappingProxyType as frozendict  # noqa: N813

import numpy as np
import pandas as pd
from natsort import natsorted
from scverse_misc import Deprecation, deprecated, deprecated_arg

from .. import _utils
from .. import logging as logg
from .._docs import (
    DEPR_ADJACENCY,
    DEPR_KEY_ADDED,
    DEPR_NEIGHBORS_KEY,
    DEPR_OBSP_GRAPH,
    doc_out_key,
)
from .._utils import _doc_params
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
    from collections.abc import Mapping, Sequence
    from typing import Any, Literal

    from anndata import AnnData

    from .._compat import CSBase
    from .._utils.random import _LegacyRandom
    from ..get._types import GraphAcc, VecRef

    try:  # sphinx-autodoc-typehints + optional dependency
        from louvain.VertexPartition import MutableVertexPartition
    except ImportError:
        if not TYPE_CHECKING:
            MutableVertexPartition = type(
                "MutableVertexPartition", (), dict(__module__="louvain.VertexPartition")
            )


@deprecated(Deprecation("1.12.0", "Use :func:`scanpy.tl.leiden` instead."))
@_doc_params(
    doc_adata=doc_adata,
    restrict_to=doc_restrict_to,
    adjacency=doc_adjacency,
    neighbors_key=doc_neighbors_key.format(method="louvain"),
    obsp=doc_obsp,
    use=doc_use_neighbors(),
    out=doc_out_key(
        "the cluster labels",
        example="A.obs['clusters']",
        key_added="`key_added='k'` is `out=A.obs['k']`.",
        default=(
            "`.obs['louvain']` and `.uns['louvain']` "
            "(`'louvain_R'` if `restrict_to` is passed)"
        ),
    ),
)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("adjacency", DEPR_ADJACENCY)
@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_arg("obsp", DEPR_OBSP_GRAPH)
@deprecated_copy
def louvain(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    resolution: float | None = None,
    *,
    random_state: _LegacyRandom = 0,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    flavor: Literal["vtraag", "igraph"] = "vtraag",
    directed: bool = True,
    use_weights: bool = False,
    partition_type: type[MutableVertexPartition] | None = None,
    partition_kwargs: Mapping[str, Any] = frozendict({}),
    use: GraphAcc | str | None = None,
    out: VecRef | bool = False,
    # deprecated
    adjacency: CSBase | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | pd.Categorical | None:
    """Cluster cells into subgroups :cite:p:`Blondel2008,Levine2015,Traag2017`.

    Cluster cells using the Louvain algorithm :cite:p:`Blondel2008` in the implementation
    of :cite:t:`Traag2017`. The Louvain algorithm was proposed for single-cell
    analysis by :cite:t:`Levine2015`.

    This requires having run :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first,
    or explicitly passing a ``adjacency`` matrix.

    .. array-support:: tl.louvain

    Parameters
    ----------
    {doc_adata}
    resolution
        For the default flavor (``'vtraag'``) or for ```RAPIDS```, you can provide a
        resolution (higher resolution means finding more and smaller clusters),
        which defaults to 1.0.
        See “Time as a resolution parameter” in :cite:t:`Lambiotte2014`.
    random_state
        Change the initialization of the optimization.
    {restrict_to}
    flavor
        Choose between to packages for computing the clustering.

        ``'vtraag'``
            Much more powerful than ``'igraph'``, and the default.
        ``'igraph'``
            Built in ``igraph`` method.
    directed
        Interpret the ``adjacency`` matrix as directed graph?
    use_weights
        Use weights from knn graph.
    partition_type
        Type of partition to use.
        Only a valid argument if ``flavor`` is ``'vtraag'``.
    partition_kwargs
        Key word arguments to pass to partitioning,
        if ``vtraag`` method is being used.
    {use}
    {out}
    {adjacency}
    {neighbors_key}
    {obsp}
    copy
        Copy adata or modify it inplace.

    Returns
    -------
    Returns the cluster labels if `out=True`, else `None`. Sets the following fields:

    `adata.obs['louvain' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['louvain' | key_added]['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.

    """
    key = _out_key(out, key_added, kind="obs", default="louvain")
    # like `key_added`, the default key gets a suffix, keys from `out` are literal
    if restrict_to is not None and key == "louvain" and out is False:
        key += "_R"
    partition_kwargs = dict(partition_kwargs)
    start = logg.info("running Louvain clustering")
    if (flavor != "vtraag") and (partition_type is not None):
        msg = '`partition_type` is only a valid argument when `flavour` is "vtraag"'
        raise ValueError(msg)
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
    if flavor in {"vtraag", "igraph"}:
        if flavor == "igraph" and resolution is not None:
            logg.warning('`resolution` parameter has no effect for flavor "igraph"')
        if directed and flavor == "igraph":
            directed = False
        if not directed:
            logg.debug("    using the undirected graph")
        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
        weights = np.array(g.es["weight"]).astype(np.float64) if use_weights else None
        if flavor == "vtraag":
            import louvain

            if partition_type is None:
                partition_type = louvain.RBConfigurationVertexPartition
            if resolution is not None:
                partition_kwargs["resolution_parameter"] = resolution
            if use_weights:
                partition_kwargs["weights"] = weights
            partition_kwargs["seed"] = random_state
            logg.info('    using the "louvain" package of Traag (2017)')
            part = louvain.find_partition(
                g,
                partition_type,
                **partition_kwargs,
            )
            # adata.uns['louvain_quality'] = part.quality()
        else:
            part = g.community_multilevel(weights=weights)
        groups = np.array(part.membership)
    elif flavor == "taynaud":
        # this is deprecated
        import community
        import networkx as nx

        g = nx.Graph(adjacency)
        partition = community.best_partition(g)
        groups = np.zeros(len(partition), dtype=int)
        for k, v in partition.items():
            groups[k] = v
    else:
        msg = '`flavor` needs to be "vtraag" or "igraph" or "taynaud".'
        raise ValueError(msg)
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
    adata.uns[key] = {}
    adata.uns[key]["params"] = dict(
        resolution=resolution,
        random_state=random_state,
    )
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return None
