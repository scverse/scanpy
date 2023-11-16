from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix


def _get_sparse_matrix_from_indices_distances(
    indices: NDArray[np.int32],
    distances: NDArray[np.float32],
    n_obs: int,
    n_neighbors: int,
) -> csr_matrix:
    n_nonzero = n_obs * n_neighbors
    D = csr_matrix(
        (
            distances.copy().ravel(),  # copy the data, otherwise strange behavior here
            indices.copy().ravel(),
            np.arange(0, n_nonzero + 1, n_neighbors),
        ),
        shape=(n_obs, n_obs),
    )
    D.eliminate_zeros()
    return D


def _get_indices_distances_from_dense_matrix(D: NDArray[np.float32], n_neighbors: int):
    sample_range = np.arange(D.shape[0])[:, None]
    indices = np.argpartition(D, n_neighbors - 1, axis=1)[:, :n_neighbors]
    indices = indices[sample_range, np.argsort(D[sample_range, indices])]
    distances = D[sample_range, indices]
    return indices, distances


def _get_indices_distances_from_sparse_matrix(
    D: csr_matrix, n_neighbors: int
) -> tuple[NDArray[np.int64], NDArray]:
    if (shortcut := _ind_dist_shortcut(D, n_neighbors)) is not None:
        return shortcut

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
    distances: csr_matrix, n_neighbors: int
) -> tuple[NDArray[np.int_], NDArray[np.float_]] | None:
    """\
    Shortcut for scipy or RAPIDS style distance matrices.

    These have `n_neighbors` or `n_neighbors + 1` entries per row.
    Therefore, this function will return matrices with either
    `n_neighbors` or `n_neighbors + 1` columns.
    """
    n_obs = distances.shape[0]  # shape is square
    # Check if we have a compatible number of entries
    if distances.nnz == n_obs * (n_neighbors + 1):
        n_neighbors += 1
    elif distances.nnz != n_obs * n_neighbors:
        return None
    # Check if each row has the correct number of entries
    # TODO: Support duplicates
    if (distances.getnnz(axis=1) != n_neighbors).any():
        return None
    return (
        distances.indices.reshape(n_obs, n_neighbors),
        distances.data.reshape(n_obs, n_neighbors),
    )
