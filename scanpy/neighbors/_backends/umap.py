from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, coo_matrix
from umap.utils import fast_knn_indices


def precomputed(
    X: NDArray[np.float32], n_neighbors: int
) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
    """X is the return value of pairwise_distances()"""
    # Note that this does not support sparse distance matrices yet ...
    # Compute indices of n nearest neighbors
    knn_indices = fast_knn_indices(X, n_neighbors)
    # knn_indices = np.argsort(X)[:, :n_neighbors]
    # Compute the nearest neighbor distances
    #   (equivalent to np.sort(X)[:,:n_neighbors])
    knn_dists = X[np.arange(X.shape[0])[:, None], knn_indices].copy()
    # Prune any nearest neighbours that are infinite distance apart.
    disconnected_index = knn_dists == np.inf
    knn_indices[disconnected_index] = -1

    return knn_indices, knn_dists


def compute_connectivities(
    knn_indices: NDArray[np.int32],
    knn_dists: NDArray[np.float32],
    *,
    n_obs: int,
    n_neighbors: int,
    set_op_mix_ratio: float = 1.0,
    local_connectivity: float = 1.0,
) -> csr_matrix:
    """\
    This is from umap.fuzzy_simplicial_set [McInnes18]_.

    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        from umap.umap_ import fuzzy_simplicial_set

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(
        X,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if isinstance(connectivities, tuple):
        # In umap-learn 0.4, this returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    return connectivities.tocsr()
