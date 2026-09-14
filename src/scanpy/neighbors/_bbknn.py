"""Batch balanced k-nearest neighbors."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 15):
    from types import MappingProxyType as frozendict  # noqa: N813

import numpy as np

from .. import logging as logg
from .._docs import doc_rng
from .._utils import _doc_params
from .._utils._doctests import doctest_needs
from ..get.get import _rep_to_json, _resolve_rep
from ._common import (
    _get_bbknn_metadata,
    _get_indices_distances_from_rect_matrix,
    _get_sparse_matrix_from_indices_distances,
    _make_transformer,
)
from ._connectivity import umap
from ._doc import doc_n_pcs, doc_use_rep

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from anndata import AnnData
    from anndata.acc import AdRef
    from numpy.typing import NDArray

    from .._compat import CSRBase
    from .._utils.random import RNGLike, SeedLike
    from ..get.get import RepAcc
    from ._types import KnnTransformerLike, _KnownTransformer, _Metric, _MetricFn


@doctest_needs("anndata_acc")
@_doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep, rng=doc_rng)
def bbknn(  # noqa: PLR0913
    adata: AnnData,
    neighbors_within_batch: int = 3,
    n_pcs: int | None = None,
    *,
    batches: AdRef | str = "obs.batch",
    use_rep: RepAcc | str | None = None,
    transformer: KnnTransformerLike | _KnownTransformer | None = None,
    metric: _Metric | _MetricFn = "euclidean",
    metric_kwds: Mapping[str, Any] = frozendict({}),
    trim: int | None = None,
    rng: SeedLike | RNGLike | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """Compute a batch balanced neighborhood graph of observations :cite:p:`Polanski2019`.

    Batch balanced kNN alters the kNN procedure to identify each cell’s top neighbors
    in each batch separately instead of the entire cell pool with no accounting for batch.
    The nearest neighbors of each batch are then merged to create a final list of
    neighbors for the cell, which aligns batches in a quick and lightweight manner.

    Use this as an alternative to :func:`~scanpy.pp.neighbors`:
    it writes the same fields, so all downstream steps
    (e.g. :func:`~scanpy.tl.umap` or :func:`~scanpy.tl.leiden`) work unchanged.
    This CPU implementation is based on the rapids-singlecell package.

    .. array-support:: pp.bbknn

    Parameters
    ----------
    adata
        Annotated data matrix.
    neighbors_within_batch
        How many top neighbors to report for each batch.
        The total number of neighbors is this number times the number of batches,
        which then serves as the basis for the construction of a symmetrical
        matrix of connectivities.
    {n_pcs}
    batches
        `adata.obs` column name discriminating between the batches.
    {use_rep}
    transformer
        kNN search backend following the API of
        :class:`~sklearn.neighbors.KNeighborsTransformer`.
        One index is built per batch and queried with all observations,
        so its ``n_neighbors`` is ignored in favor of ``neighbors_within_batch``.
        See :doc:`/how-to/knn-transformers` for more details.
        Also accepts the following known options:

        `None` (the default)
            Behavior depends on data size.
            For small data, we will calculate exact kNN, otherwise we use
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'pynndescent'`
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
    metric
        A known metric’s name or a callable that returns a distance.

        *ignored if* ``transformer`` *is an instance.*
    metric_kwds
        Options for the metric.

        *ignored if* ``transformer`` *is an instance.*
    trim
        Trim each cell’s neighbors to these many top connectivities.
        May help with population independence and improve the tidiness of clustering.
        The lower the value, the more independent the individual populations,
        at the cost of a more conserved batch effect.
        If `None`, this is set to 10 times the total number of neighbors.
        Set to 0 to skip trimming.
    {rng}

        *ignored if* ``transformer`` *is an instance.*
    key_added
        If not specified, the neighbors data is stored in `.uns['neighbors']`,
        distances and connectivities are stored in `.obsp['distances']` and
        `.obsp['connectivities']` respectively.
        If specified, the neighbors data is added to .uns[key_added],
        distances are stored in `.obsp[f'{{key_added}}_distances']` and
        connectivities in `.obsp[f'{{key_added}}_connectivities']`.
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsp['distances' | f'{{key_added}}_distances']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Distance matrix of the batch balanced nearest neighbors search.
        Each row (cell) has ``neighbors_within_batch`` × ``n_batches`` - 1 non-zero entries:
        its nearest neighbors in each batch, excluding the cell itself.
    `adata.obsp['connectivities' | f'{{key_added}}_connectivities']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    `adata.uns['neighbors' | key_added]` : :class:`dict`
        neighbors parameters.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.obs["batch"] = adata.obs["phase"]
    >>> sc.pp.bbknn(adata, batches="obs.batch")
    >>> sc.tl.umap(adata)

    See Also
    --------
    :func:`~scanpy.pp.neighbors`
    :doc:`/how-to/knn-transformers`

    """
    from anndata.acc import A, AdRef

    from ..tools._utils import _choose_representation

    start = logg.info("computing batch balanced neighbors")

    adata = adata.copy() if copy else adata
    if adata.is_view:  # we shouldn’t need this here...
        adata._init_as_actual(adata.copy())

    if not isinstance(batches, AdRef):
        batches = A.resolve(batches, vec=True)
    if use_rep is not None:
        use_rep = _resolve_rep(use_rep)

    if neighbors_within_batch < 1:
        msg = "`neighbors_within_batch` needs to be greater than 0."
        raise ValueError(msg)
    if batches not in adata:
        msg = f"Batch key {batches!r} not found in `adata.obs`."
        raise KeyError(msg)
    rng = np.random.default_rng(rng)

    batch_arr = np.asarray(adata[batches])
    unique_batches, batch_sizes = np.unique(batch_arr, return_counts=True)
    if len(too_small := unique_batches[batch_sizes < neighbors_within_batch]):
        msg = (
            f"Not all batches have at least `neighbors_within_batch = "
            f"{neighbors_within_batch}` cells in them: {list(too_small)}."
        )
        raise ValueError(msg)

    x = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
    knn_indices, knn_distances = _compute_batch_balanced_knn(
        x,
        batches=batch_arr,
        unique_batches=unique_batches,
        batch_sizes=batch_sizes,
        neighbors_within_batch=neighbors_within_batch,
        transformer=transformer,
        metric=metric,
        metric_kwds=metric_kwds,
        rng=rng,
    )
    n_obs, n_neighbors = knn_indices.shape
    if trim is None:
        trim = 10 * n_neighbors

    distances = _get_sparse_matrix_from_indices_distances(
        knn_indices, knn_distances, keep_self=False
    )
    start_connect = logg.debug("computed batch balanced neighbors", time=start)
    connectivities = umap(
        knn_indices, knn_distances, n_obs=n_obs, n_neighbors=n_neighbors
    )
    if trim > 0:
        connectivities = _trim(connectivities, trim)
    logg.debug("computed connectivities", time=start_connect)

    key_added, neighbors_dict = _get_bbknn_metadata(
        key_added,
        n_neighbors=n_neighbors,
        method="umap",
        metric=metric,
        **({} if not metric_kwds else dict(metric_kwds=metric_kwds)),
        **({} if use_rep is None else dict(use_rep=_rep_to_json(use_rep))),
        **({} if n_pcs is None else dict(n_pcs=n_pcs)),
        batches=A.to_json(batches),
        neighbors_within_batch=neighbors_within_batch,
        trim=trim,
    )
    adata.uns[key_added] = neighbors_dict
    adata.obsp[neighbors_dict["distances_key"]] = distances
    adata.obsp[neighbors_dict["connectivities_key"]] = connectivities

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added to `.uns[{key_added!r}]`\n"
            f"    `.obsp[{neighbors_dict['distances_key']!r}]`, distances for each pair of neighbors\n"
            f"    `.obsp[{neighbors_dict['connectivities_key']!r}]`, weighted adjacency matrix"
        ),
    )
    return adata if copy else None


def _compute_batch_balanced_knn(
    x: NDArray[np.float32 | np.float64] | CSRBase,
    /,
    *,
    batches: NDArray[Any],
    unique_batches: NDArray[Any],
    batch_sizes: NDArray[np.int64],
    neighbors_within_batch: int,
    transformer: KnnTransformerLike | _KnownTransformer | None,
    metric: _Metric | _MetricFn,
    metric_kwds: Mapping[str, Any],
    rng: np.random.Generator,
) -> tuple[NDArray[np.int64], NDArray[np.float32 | np.float64]]:
    """Find the `neighbors_within_batch` nearest neighbors of each cell in each batch.

    Returns the merged indices and distances, sorted by distance within each row.
    """
    from sklearn.base import clone

    proto, is_sklearn_shortcut = _handle_transformer(
        transformer,
        n_obs=x.shape[0],
        max_batch_size=int(batch_sizes.max()),
        n_neighbors=neighbors_within_batch,
        metric=metric,
        metric_kwds=metric_kwds,
        rng=rng,
    )

    n_obs = x.shape[0]
    n_neighbors = neighbors_within_batch * len(unique_batches)
    knn_indices = np.empty((n_obs, n_neighbors), dtype=np.int64)
    knn_distances = np.empty((n_obs, n_neighbors), dtype=np.float64)
    obs_indices = np.arange(n_obs)

    for i, batch in enumerate(unique_batches):
        mask = batches == batch
        knn = clone(proto).fit(x[mask])
        d = (
            # ask for exactly `neighbors_within_batch`; `transform` would add one,
            # which fails for batches that only have `neighbors_within_batch` cells
            knn.kneighbors_graph(x, n_neighbors=neighbors_within_batch, mode="distance")
            if is_sklearn_shortcut
            else knn.transform(x)
        )
        indices, distances = _get_indices_distances_from_rect_matrix(
            d, neighbors_within_batch
        )
        cols = slice(i * neighbors_within_batch, (i + 1) * neighbors_within_batch)
        # the transformer’s indices are relative to the batch
        knn_indices[:, cols] = obs_indices[mask][indices]
        knn_distances[:, cols] = distances
        logg.debug(f"    computed neighbors within batch {batch!r}")

    # some backends report a tiny non-zero distance of a cell to itself,
    # which `umap` would mistake for the radius of the cell’s local neighborhood
    is_self = knn_indices == obs_indices[:, None]
    knn_distances[is_self] = 0.0

    # `umap` derives each cell’s local connectivity from its closest neighbors,
    # so the merged rows need to be sorted by distance.
    # Ties are broken in favor of the cell itself, which is dropped from `.obsp['distances']`.
    order = np.lexsort((~is_self, knn_distances), axis=1)
    return (
        np.take_along_axis(knn_indices, order, axis=1),
        np.take_along_axis(knn_distances, order, axis=1),
    )


def _handle_transformer(
    transformer: KnnTransformerLike | _KnownTransformer | None,
    *,
    n_obs: int,
    max_batch_size: int,
    n_neighbors: int,
    metric: _Metric | _MetricFn,
    metric_kwds: Mapping[str, Any],
    rng: np.random.Generator,
) -> tuple[KnnTransformerLike, bool]:
    """Coerce `transformer` to an instance to be cloned for each batch.

    Also returns whether it is a :class:`~sklearn.neighbors.KNeighborsTransformer`
    we created ourselves, i.e. one we can query without going through ``transform``.

    Unlike :func:`~scanpy.pp.neighbors`,
    we build one index per batch and query each with all `n_obs` observations,
    so brute force costs ``n_obs × max_batch_size``,
    while an approximate index’s cost is dominated by building it.
    The cutoff is where the two met in benchmarks on ~50-dimensional data.
    """
    shortcut = transformer == "sklearn" or (
        transformer is None
        and (
            max_batch_size < 4096
            or (metric == "euclidean" and n_obs * max_batch_size < 10**9)
        )
    )
    return _make_transformer(
        transformer,
        shortcut=shortcut,
        n_index=max_batch_size,
        n_neighbors=n_neighbors,
        metric=metric,
        metric_params=metric_kwds,
        rng=rng,
    ), shortcut


def _trim(connectivities: CSRBase, /, trim: int) -> CSRBase:
    """Trim the graph in place to the `trim` strongest connections per cell.

    Following the reference implementation, an edge is dropped if its weight is
    below the `trim`-th largest weight of *either* of the cells it connects,
    which keeps the graph symmetric.
    """
    n_nonzero = np.diff(connectivities.indptr)
    if not (n_nonzero > trim).any():
        return connectivities
    rows = np.repeat(np.arange(connectivities.shape[0]), n_nonzero)
    # sort each row’s weights in descending order to find its `trim`-th largest one.
    # rows with at most `trim` entries have no such weight and keep a cutoff of 0.
    order = np.lexsort((-connectivities.data, rows))
    rank_in_row = np.arange(connectivities.nnz) - np.repeat(
        connectivities.indptr[:-1], n_nonzero
    )
    at_cutoff = rank_in_row == trim - 1
    cutoffs = np.zeros(connectivities.shape[0], dtype=connectivities.data.dtype)
    cutoffs[rows[at_cutoff]] = connectivities.data[order][at_cutoff]

    keep_above = np.maximum(cutoffs[rows], cutoffs[connectivities.indices])
    connectivities.data[connectivities.data < keep_above] = 0
    connectivities.eliminate_zeros()
    return connectivities
