from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from fast_array_utils.stats import is_constant
from scipy import sparse

from .._compat import warn
from .._settings import settings
from .._utils import get_literal_vals
from .._utils.random import _rng_kwds
from ..get._out import _dest_key
from ..get._types import GraphAcc
from ..get._use import _resolve_obs
from ._types import NeighborsDict, NeighborsOut, _KnownTransformer

if TYPE_CHECKING:
    from typing import Unpack

    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import CSBase, CSRBase
    from .._settings import Default
    from .._utils import NeighborsView
    from ..get._out import OutKind
    from ..get._types import _SetRep
    from ._types import (
        BbknnParams,
        KnnTransformerLike,
        KwdsForTransformer,
        NeighborsParams,
    )


def _make_transformer(
    transformer: KnnTransformerLike | _KnownTransformer | None,
    /,
    *,
    shortcut: bool,
    n_index: int,
    **kwds: Unpack[KwdsForTransformer],
) -> KnnTransformerLike:
    """Coerce `transformer` from `None` or a string to an instance.

    `shortcut` requests a brute force :class:`~sklearn.neighbors.KNeighborsTransformer`;
    the caller decides that, as the trade-off depends on how the index is queried.
    Otherwise, `transformer=None` is set up like `umap` does, i.e. to a
    ~`pynndescent.PyNNDescentTransformer` with custom `n_trees` and `n_iters`
    derived from `n_index`, the number of observations in the index.

    `rng` is passed on as whichever of `rng`/`random_state` the class accepts.
    """
    rng = kwds.pop("rng", None)
    if shortcut:
        from sklearn.neighbors import KNeighborsTransformer

        assert transformer in {None, "sklearn"}
        return KNeighborsTransformer(
            algorithm="brute",
            n_jobs=settings.n_jobs,
            n_neighbors=kwds["n_neighbors"],
            metric=kwds["metric"],
            metric_params=dict(kwds["metric_params"]),  # needs dict
            **_rng_kwds(KNeighborsTransformer, rng),
        )
    if transformer is None or transformer == "pynndescent":
        from pynndescent import PyNNDescentTransformer

        kwds["metric_kwds"] = dict(kwds.pop("metric_params"))  # needs to be cloneable
        if transformer is None:
            # Use defaults from UMAP’s `nearest_neighbors` function
            kwds.update(
                n_jobs=settings.n_jobs,
                n_trees=min(64, 5 + round(n_index**0.5 / 20.0)),
                n_iters=max(5, round(np.log2(n_index))),
            )
        return PyNNDescentTransformer(**kwds, **_rng_kwds(PyNNDescentTransformer, rng))
    if isinstance(transformer, str):
        msg = (
            f"Unknown transformer: {transformer}. "
            f"Try passing a class or one of {get_literal_vals(_KnownTransformer)}"
        )
        raise ValueError(msg)
    return transformer  # `transformer` is probably an instance


NEIGHBORS_KINDS: NeighborsOut[set[OutKind]] = NeighborsOut({"obsp"}, {"obsp"})
"""Where `out` can write the outputs of `pp.neighbors` and `pp.bbknn`."""


def _neighbors_dests(key_added: str | Default | None) -> NeighborsOut[_SetRep | None]:
    """Where `key_added` put the outputs of `pp.neighbors` and `pp.bbknn`."""
    prefix = "" if key_added is None else f"{key_added}_"
    return NeighborsOut(
        connectivities=dict(obsp=f"{prefix}connectivities"),
        distances=dict(obsp=f"{prefix}distances"),
    )


def _get_metadata(
    key: str | None,
    dests: NeighborsOut[_SetRep | None],
    /,
    **params: Unpack[NeighborsParams],
) -> tuple[str, NeighborsDict[NeighborsParams]]:
    return _metadata(key, dests, params)


def _get_bbknn_metadata(
    key: str | None,
    dests: NeighborsOut[_SetRep | None],
    /,
    **params: Unpack[BbknnParams],
) -> tuple[str, NeighborsDict[BbknnParams]]:
    return _metadata(key, dests, params)


def _metadata[P: NeighborsParams](
    key: str | None, dests: NeighborsOut[_SetRep | None], params: P
) -> tuple[str, NeighborsDict[P]]:
    """Get the `.uns` key and entry, which names where the matrices were written."""
    neighbors_dict = NeighborsDict(params=params)
    for name, dest in zip(NeighborsOut._fields, dests, strict=True):
        if (dest_key := _dest_key(dest)) is not None:
            neighbors_dict[f"{name}_key"] = dest_key
    return key or "neighbors", neighbors_dict


def _has_self_column(
    indices: NDArray[np.int32 | np.int64],
    distances: NDArray[np.float32 | np.float64],
) -> bool:
    # some algorithms have some messed up reordering.
    return (indices[:, 0] == np.arange(indices.shape[0])).any()


def _remove_self_column(
    indices: NDArray[np.int32 | np.int64],
    distances: NDArray[np.float32 | np.float64],
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    if not _has_self_column(indices, distances):
        msg = "The first neighbor should be the cell itself."
        raise AssertionError(msg)
    return indices[:, 1:], distances[:, 1:]


def _get_sparse_matrix_from_indices_distances(
    indices: NDArray[np.int32 | np.int64],
    distances: NDArray[np.float32 | np.float64],
    *,
    keep_self: bool,
) -> CSRBase:
    """Create a sparse matrix from a pair of indices and distances.

    If keep_self=False, it verifies that the first column is the cell itself,
    then removes it from the explicitly stored zeroes.

    Duplicates in the data are kept as explicitly stored zeroes.
    """
    # instead of calling .eliminate_zeros() on our sparse matrix,
    # we manually handle the nearest neighbor being the cell itself.
    # This allows us to use _ind_dist_shortcut even when the data has duplicates.
    if not keep_self:
        indices, distances = _remove_self_column(indices, distances)
    indptr = np.arange(0, np.prod(indices.shape) + 1, indices.shape[1])
    return sparse.csr_matrix(  # noqa: TID251
        (
            distances.copy().ravel(),  # copy the data, otherwise strange behavior here
            indices.copy().ravel(),
            indptr,
        ),
        shape=(indices.shape[0],) * 2,
    )


def _get_indices_distances_from_dense_matrix(
    d: NDArray[np.float32 | np.float64], /, n_neighbors: int
):
    sample_range = np.arange(d.shape[0])[:, None]
    indices = np.argpartition(d, n_neighbors - 1, axis=1)[:, :n_neighbors]
    indices = indices[sample_range, np.argsort(d[sample_range, indices])]
    distances = d[sample_range, indices]
    return indices, distances


def _get_indices_distances_from_sparse_matrix(
    d: CSRBase, /, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    """Get indices and distances from a sparse matrix.

    Makes sure that for both of the returned matrices:
    1. the first column corresponds to the cell itself as nearest neighbor.
    2. the number of neighbors (`.shape[1]`) is restricted to `n_neighbors`.
    """
    if (shortcut := _ind_dist_shortcut(d)) is not None:
        indices, distances = shortcut
    else:
        indices, distances = _ind_dist_slow(d, n_neighbors)

    # handle RAPIDS style indices_distances lacking the self-column
    if not _has_self_column(indices, distances):
        indices = np.hstack([np.arange(indices.shape[0])[:, None], indices])
        distances = np.hstack([np.zeros(distances.shape[0])[:, None], distances])

    # If using the shortcut or adding the self column resulted in too many neighbors,
    # restrict the output matrices to the correct size
    if indices.shape[1] > n_neighbors:
        indices, distances = indices[:, :n_neighbors], distances[:, :n_neighbors]

    return indices, distances


def _ind_dist_slow(
    d: CSRBase, /, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    indices = np.zeros((d.shape[0], n_neighbors), dtype=int)
    distances = np.zeros((d.shape[0], n_neighbors), dtype=d.dtype)
    n_neighbors_m1 = n_neighbors - 1
    for i in range(indices.shape[0]):
        neighbors = d[i].nonzero()  # 'true' and 'spurious' zeros
        indices[i, 0] = i
        distances[i, 0] = 0
        # account for the fact that there might be more than n_neighbors
        # due to an approximate search
        # [the point itself was not detected as its own neighbor during the search]
        if len(neighbors[1]) > n_neighbors_m1:
            sorted_indices = np.argsort(d[i][neighbors].A1)[:n_neighbors_m1]
            indices[i, 1:] = neighbors[1][sorted_indices]
            distances[i, 1:] = d[i][
                neighbors[0][sorted_indices], neighbors[1][sorted_indices]
            ]
        else:
            indices[i, 1:] = neighbors[1]
            distances[i, 1:] = d[i][neighbors]
    return indices, distances


def _ind_dist_shortcut(
    d: CSRBase, /
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]] | None:
    """Shortcut for scipy or RAPIDS style distance matrices."""
    # Check if each row has the correct number of entries
    nnzs = d.getnnz(axis=1)
    if not is_constant(nnzs):
        msg = (
            "Sparse matrix has no constant number of neighbors per row. "
            "Cannot efficiently get indices and distances."
        )
        warn(msg, RuntimeWarning)
        return None
    n_obs, n_neighbors = d.shape[0], int(nnzs[0])
    return (
        d.indices.reshape(n_obs, n_neighbors),
        d.data.reshape(n_obs, n_neighbors),
    )


def _get_indices_distances_from_rect_matrix(
    d: CSRBase, /, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    """Get the `n_neighbors` nearest neighbors from a rectangular kNN distance matrix.

    In contrast to `_get_indices_distances_from_sparse_matrix`,
    the columns of `d` index a subset of the observations the rows index,
    so there is no self-column to take care of.
    Rows are sorted by distance and truncated to `n_neighbors` entries.
    """
    n_nonzero = np.diff(d.indptr)
    if (n_too_few := int((n_nonzero < n_neighbors).sum())) > 0:
        msg = (
            f"The transformer returned fewer than {n_neighbors} neighbors "
            f"for {n_too_few} of {d.shape[0]} observations."
        )
        raise ValueError(msg)
    if is_constant(n_nonzero):
        n_cols = int(n_nonzero[0])
        indices = d.indices.reshape(d.shape[0], n_cols)
        distances = d.data.reshape(d.shape[0], n_cols)
    else:  # pad the rows to a common width, sorting the padding to the end
        indices = np.zeros((d.shape[0], int(n_nonzero.max())), dtype=d.indices.dtype)
        distances = np.full(indices.shape, np.inf, dtype=d.data.dtype)
        rows = np.repeat(np.arange(d.shape[0]), n_nonzero)
        cols = np.arange(d.nnz) - np.repeat(d.indptr[:-1], n_nonzero)
        indices[rows, cols] = d.indices
        distances[rows, cols] = d.data
    order = np.argsort(distances, axis=1, kind="stable")[:, :n_neighbors]
    return (
        np.take_along_axis(indices, order, axis=1),
        np.take_along_axis(distances, order, axis=1),
    )


def _graph_key(use: GraphAcc | str) -> str:
    """Get the `.obsp` key of a graph given as `use`."""
    graph = _resolve_obs(use)
    if not isinstance(graph, GraphAcc):
        msg = f"`use` must point to `.obsp`, was {use!r}"
        raise TypeError(msg)
    return graph.k


def _neighbors_view(
    adata: AnnData, use: GraphAcc | str | None, neighbors_key: str | None = None
) -> NeighborsView | None:
    """View the neighbors around the graph `use`, or the deprecated `neighbors_key`’s.

    Returns `None` if there are no neighbors at the default place.
    """
    from .._utils import NeighborsView

    if use is not None:
        if neighbors_key is not None:
            msg = "Pass either `use` or `neighbors_key`, not both."
            raise TypeError(msg)
        return NeighborsView.from_graph(adata, _graph_key(use))
    key = neighbors_key or "neighbors"
    return NeighborsView(adata, key) if key in adata.uns else None


def _use_graph(
    adata: AnnData,
    use: GraphAcc | str | None,
    *,
    adjacency: CSBase | None = None,
    obsp: str | None = None,
    neighbors_key: str | None = None,
) -> CSBase:
    """Choose the connectivities from `use`, or the deprecated `adjacency`/`obsp`/`neighbors_key`."""
    from .._utils import _choose_graph

    if (adjacency is not None) + (obsp is not None) + (use is not None) > 1:
        msg = "Pass only one of `use`, `obsp`, or `adjacency`."
        raise TypeError(msg)
    if adjacency is not None:
        return adjacency
    if use is not None:
        if neighbors_key is not None:
            msg = "Pass either `use` or `neighbors_key`, not both."
            raise TypeError(msg)
        return adata.obsp[_graph_key(use)]
    return _choose_graph(adata, obsp, neighbors_key)
