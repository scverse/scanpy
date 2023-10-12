from __future__ import annotations

import warnings
from typing import Union

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import issparse, csr_matrix, coo_matrix

from ._common import (
    _get_indices_distances_from_dense_matrix,
    _get_indices_distances_from_sparse_matrix,
)


def gauss(
    distances: Union[NDArray[np.float32], csr_matrix], n_neighbors: int, *, knn: bool
):
    """
    Derive gaussian connectivities between data points from their distances.

    Parameters
    ----------
    distances
        The input matrix of distances between data points.
    n_neighbors
        The number of nearest neighbors to consider.
    knn
        Specify if the distances have been restricted to k nearest neighbors.
    """
    # init distances
    if issparse(distances):
        assert isinstance(distances, csr_matrix)
        Dsq = distances.power(2)
        indices, distances_sq = _get_indices_distances_from_sparse_matrix(
            Dsq, n_neighbors
        )
    else:
        assert isinstance(distances, np.ndarray)
        Dsq = np.power(distances, 2)
        indices, distances_sq = _get_indices_distances_from_dense_matrix(
            Dsq, n_neighbors
        )

    # exclude the first point, the 0th neighbor
    indices = indices[:, 1:]
    distances_sq = distances_sq[:, 1:]

    # choose sigma, the heuristic here doesn't seem to make much of a difference,
    # but is used to reproduce the figures of Haghverdi et al. (2016)
    if issparse(distances):
        # as the distances are not sorted
        # we have decay within the n_neighbors first neighbors
        sigmas_sq = np.median(distances_sq, axis=1)
    else:
        # the last item is already in its sorted position through argpartition
        # we have decay beyond the n_neighbors neighbors
        sigmas_sq = distances_sq[:, -1] / 4
    sigmas = np.sqrt(sigmas_sq)

    # compute the symmetric weight matrix
    if not issparse(distances):
        Num = 2 * np.multiply.outer(sigmas, sigmas)
        Den = np.add.outer(sigmas_sq, sigmas_sq)
        W = np.sqrt(Num / Den) * np.exp(-Dsq / Den)
        # make the weight matrix sparse
        if not knn:
            mask = W > 1e-14
            W[~mask] = 0
        else:
            # restrict number of neighbors to ~k
            # build a symmetric mask
            mask = np.zeros(Dsq.shape, dtype=bool)
            for i, row in enumerate(indices):
                mask[i, row] = True
                for j in row:
                    if i not in set(indices[j]):
                        W[j, i] = W[i, j]
                        mask[j, i] = True
            # set all entries that are not nearest neighbors to zero
            W[~mask] = 0
    else:
        assert isinstance(Dsq, csr_matrix)
        W = Dsq.copy()  # need to copy the distance matrix here; what follows is inplace
        for i in range(len(Dsq.indptr[:-1])):
            row = Dsq.indices[Dsq.indptr[i] : Dsq.indptr[i + 1]]
            num = 2 * sigmas[i] * sigmas[row]
            den = sigmas_sq[i] + sigmas_sq[row]
            W.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] = np.sqrt(num / den) * np.exp(
                -Dsq.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] / den
            )
        W = W.tolil()
        for i, row in enumerate(indices):
            for j in row:
                if i not in set(indices[j]):
                    W[j, i] = W[i, j]
        W = W.tocsr()

    return W


def umap(
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
