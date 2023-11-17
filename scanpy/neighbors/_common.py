from __future__ import annotations
from math import dist

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix


def _has_self_column(
    indices: NDArray[np.int32 | np.int64],
    distances: NDArray[np.float32 | np.float64],
) -> bool:
    return (distances[:, 0] == 0.0).all() and (
        indices[:, 0] == np.arange(indices.shape[0])
    ).all()


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
) -> csr_matrix:
    """\
    Create a sparse matrix from a pair of indices and distances.

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
    return csr_matrix(
        (
            distances.copy().ravel(),  # copy the data, otherwise strange behavior here
            indices.copy().ravel(),
            indptr,
        ),
        shape=(indices.shape[0],) * 2,
    )


def _get_indices_distances_from_dense_matrix(
    D: NDArray[np.float32 | np.float64], n_neighbors: int
):
    sample_range = np.arange(D.shape[0])[:, None]
    indices = np.argpartition(D, n_neighbors - 1, axis=1)[:, :n_neighbors]
    indices = indices[sample_range, np.argsort(D[sample_range, indices])]
    distances = D[sample_range, indices]
    return indices, distances


def _get_indices_distances_from_sparse_matrix(
    D: csr_matrix, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    """\
    Get indices and distances from a sparse matrix.

    Makes sure the first column corresponds to the cell itself as nearest neighbor.
    """
    if (shortcut := _ind_dist_shortcut(D, n_neighbors)) is not None:
        indices, distances = shortcut
    else:
        indices, distances = _ind_dist_slow(D, n_neighbors)

    # handle RAPIDS style indices_distances lacking the self-column
    if not _has_self_column(indices, distances):
        indices = np.hstack([np.arange(indices.shape[0])[:, None], indices])
        distances = np.hstack([np.zeros(distances.shape[0])[:, None], distances])

    return indices, distances


def _ind_dist_slow(
    D: csr_matrix, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]]:
    indices = np.zeros((D.shape[0], n_neighbors), dtype=int)
    distances = np.zeros((D.shape[0], n_neighbors), dtype=D.dtype)
    n_neighbors_m1 = n_neighbors - 1
    for i in range(indices.shape[0]):
        neighbors = D[i].nonzero()  # 'true' and 'spurious' zeros
        indices[i, 0] = i
        distances[i, 0] = 0
        # account for the fact that there might be more than n_neighbors
        # due to an approximate search
        # [the point itself was not detected as its own neighbor during the search]
        if len(neighbors[1]) > n_neighbors_m1:
            sorted_indices = np.argsort(D[i][neighbors].A1)[:n_neighbors_m1]
            indices[i, 1:] = neighbors[1][sorted_indices]
            distances[i, 1:] = D[i][
                neighbors[0][sorted_indices], neighbors[1][sorted_indices]
            ]
        else:
            indices[i, 1:] = neighbors[1]
            distances[i, 1:] = D[i][neighbors]
    return indices, distances


def _ind_dist_shortcut(
    D: csr_matrix, n_neighbors: int
) -> tuple[NDArray[np.int32 | np.int64], NDArray[np.float32 | np.float64]] | None:
    """\
    Shortcut for scipy or RAPIDS style distance matrices.

    These have `n_neighbors` or `n_neighbors + 1` entries per row.
    Therefore, this function will return matrices with either
    `n_neighbors` or `n_neighbors + 1` columns.
    """
    n_obs = D.shape[0]  # shape is square
    # Check if we have a compatible number of entries
    if D.nnz == n_obs * (n_neighbors + 1):
        n_neighbors += 1
    elif D.nnz != n_obs * n_neighbors:
        return None
    # Check if each row has the correct number of entries
    if (D.getnnz(axis=1) != n_neighbors).any():
        return None
    return (
        D.indices.reshape(n_obs, n_neighbors),
        D.data.reshape(n_obs, n_neighbors),
    )
