import numpy as np
from scipy.sparse import issparse

from ..common import _get_indices_distances_from_dense_matrix


def _get_indices_distances_from_sparse_matrix(D, n_neighbors: int):
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


def _compute_connectivities_diffmap(
    distances, n_neighbors, knn, density_normalize=True
):
    # init distances
    if knn:
        Dsq = distances.power(2)
        indices, distances_sq = _get_indices_distances_from_sparse_matrix(
            Dsq, n_neighbors
        )
    else:
        Dsq = np.power(distances, 2)
        indices, distances_sq = _get_indices_distances_from_dense_matrix(
            Dsq, n_neighbors
        )

    # exclude the first point, the 0th neighbor
    indices = indices[:, 1:]
    distances_sq = distances_sq[:, 1:]

    # choose sigma, the heuristic here doesn't seem to make much of a difference,
    # but is used to reproduce the figures of Haghverdi et al. (2016)
    if knn:
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
