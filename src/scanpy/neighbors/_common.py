from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
from fast_array_utils.stats import is_constant
from scipy import sparse

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .._compat import CSRBase


def _has_self_column(
    indices: NDArray[np.int32 | np.int64],
) -> bool:
    # some algorithms have some messed up reordering.
    return (indices[:, 0] == np.arange(indices.shape[0])).any()


def _remove_self_column(
    indices: NDArray[np.int32 | np.int64],
    distances: NDArray[np.float32 | np.float64],
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    if not _has_self_column(indices):
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
    if not _has_self_column(indices):
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
    d = d.tocsr()
    n_obs = d.shape[0]
    indices = np.zeros((n_obs, n_neighbors), dtype=int)
    indices[:, 0] = np.arange(n_obs)  # set self-indices
    distances = np.zeros((n_obs, n_neighbors), dtype=d.dtype)
    n_neighbors_m1 = n_neighbors - 1

    for i in range(n_obs):
        row = d.getrow(i)
        row_indices = row.indices
        row_data = row.data

        if len(row_indices) == 0:
            continue

        # self-indices are preset, so we filter them out here
        mask = row_indices != i
        non_self_indices = row_indices[mask]
        non_self_data = row_data[mask]

        # account for the fact that there might be more than n_neighbors
        # due to an approximate search
        # [the point itself was not detected as its own neighbor during the search]
        if len(non_self_indices) > n_neighbors_m1:
            # partial sorting to get the n_neighbors - 1 smallest distances
            partition_indices = np.argpartition(non_self_data, n_neighbors_m1 - 1)[
                :n_neighbors_m1
            ]
            sorted_partition = partition_indices[
                np.argsort(non_self_data[partition_indices])
            ]
            indices[i, 1:] = non_self_indices[sorted_partition]
            distances[i, 1:] = non_self_data[sorted_partition]
        else:
            n_actual = len(non_self_indices)
            if n_actual > 0:
                indices[i, 1 : 1 + n_actual] = non_self_indices
                distances[i, 1 : 1 + n_actual] = non_self_data

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
        # 4: caller -> 3: `Neighbors.compute_neighbors` -> 2: `_get_indices_distances_from_sparse_matrix` -> 1: here
        warn(msg, category=RuntimeWarning, stacklevel=4)
        return None
    n_obs, n_neighbors = d.shape[0], int(nnzs[0])
    return (
        d.indices.reshape(n_obs, n_neighbors),
        d.data.reshape(n_obs, n_neighbors),
    )
