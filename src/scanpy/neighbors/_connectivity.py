from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from .._compat import CSRBase
from ._common import (
    _get_indices_distances_from_dense_matrix,
    _get_indices_distances_from_sparse_matrix,
    _get_sparse_matrix_from_indices_distances,
)


def gauss[D: (NDArray[np.float32 | np.float64], CSRBase)](  # noqa: PLR0912
    distances: D, n_neighbors: int, *, knn: bool
) -> D:
    """Derive gaussian connectivities between data points from their distances.

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
    if isinstance(distances, CSRBase):
        d_sq = distances.power(2)
        indices, distances_sq = _get_indices_distances_from_sparse_matrix(
            d_sq, n_neighbors
        )
    else:
        assert isinstance(distances, np.ndarray)
        d_sq = np.power(distances, 2)
        indices, distances_sq = _get_indices_distances_from_dense_matrix(
            d_sq, n_neighbors
        )

    # exclude the first point, the 0th neighbor
    indices = indices[:, 1:]
    distances_sq = distances_sq[:, 1:]

    # choose sigma, the heuristic here doesn't seem to make much of a difference,
    # but is used to reproduce the figures of Haghverdi et al. (2016)
    if isinstance(distances, CSRBase):
        # as the distances are not sorted
        # we have decay within the n_neighbors first neighbors
        sigmas_sq = np.median(distances_sq, axis=1)
    else:
        # the last item is already in its sorted position through argpartition
        # we have decay beyond the n_neighbors neighbors
        sigmas_sq = distances_sq[:, -1] / 4
    sigmas = np.sqrt(sigmas_sq)

    # compute the symmetric weight matrix
    if not isinstance(distances, CSRBase):
        num = 2 * np.multiply.outer(sigmas, sigmas)
        den = np.add.outer(sigmas_sq, sigmas_sq)
        w = np.sqrt(num / den) * np.exp(-d_sq / den)
        # make the weight matrix sparse
        if not knn:
            mask = w > 1e-14
            w[~mask] = 0
        else:
            # restrict number of neighbors to ~k
            # build a symmetric mask
            mask = np.zeros(d_sq.shape, dtype=bool)
            for i, row in enumerate(indices):
                mask[i, row] = True
                for j in row:
                    if i not in set(indices[j]):
                        w[j, i] = w[i, j]
                        mask[j, i] = True
            # set all entries that are not nearest neighbors to zero
            w[~mask] = 0
    else:
        assert isinstance(d_sq, CSRBase)
        # need to copy the distance matrix here; what follows is inplace
        w = d_sq.copy()
        for i in range(len(d_sq.indptr[:-1])):
            row = d_sq.indices[d_sq.indptr[i] : d_sq.indptr[i + 1]]
            num = 2 * sigmas[i] * sigmas[row]
            den = sigmas_sq[i] + sigmas_sq[row]
            w.data[d_sq.indptr[i] : d_sq.indptr[i + 1]] = np.sqrt(num / den) * np.exp(
                -d_sq.data[d_sq.indptr[i] : d_sq.indptr[i + 1]] / den
            )
        w = w.tolil()
        for i, row in enumerate(indices):
            for j in row:
                if i not in set(indices[j]):
                    w[j, i] = w[i, j]
        w = w.tocsr()

    return w


def umap(
    knn_indices: NDArray[np.int32 | np.int64],
    knn_dists: NDArray[np.float32 | np.float64],
    *,
    n_obs: int,
    n_neighbors: int,
    set_op_mix_ratio: float = 1.0,
    local_connectivity: float = 1.0,
) -> CSRBase:
    """Wrap for `umap.fuzzy_simplicial_set` :cite:p:`McInnes2018`.

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

    x = sparse.coo_matrix((n_obs, 1))
    connectivities, _sigmas, _rhos = fuzzy_simplicial_set(
        x,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    return connectivities.tocsr()


def jaccard(
    knn_indices: NDArray[np.int32 | np.int64],
    *,
    n_obs: int,
    n_neighbors: int,
) -> CSRBase:
    """Derive Jaccard connectivities between data points from kNN indices.

    Re-implements the weighting method from Phenograph, :cite:p:`Levine2015`.

    Parameters
    ----------
    knn_indices
        The input matrix of nearest neighbor indices for each cell.
    n_obs
        Number of cells in the data-set.
    n_neighbors
        The number of nearest neighbors to consider.

    """
    # Construct unweighted kNN adjacency matrix (self excluded, as in PhenoGraph)
    adjacency = _get_sparse_matrix_from_indices_distances(
        knn_indices, np.ones_like(knn_indices), keep_self=False
    )

    # Compute |N(i) ∩ N(j)|
    i_idx = np.repeat(np.arange(n_obs), n_neighbors - 1)
    j_idx = knn_indices[:, 1:].ravel()
    rows_i = adjacency[i_idx, :]
    rows_j = adjacency[j_idx, :]
    shared = np.asarray(rows_i.multiply(rows_j).sum(axis=1)).ravel()

    # Jaccard index
    jaccard = shared / (2 * (n_neighbors - 1) - shared)

    # Build connectivity matrix, filtering out zeros
    mask = jaccard != 0
    connectivities = sparse.csr_matrix(  # noqa: TID251
        (jaccard[mask], (i_idx[mask], j_idx[mask])),
        shape=(n_obs, n_obs),
    )

    # Symmetrize by averaging (as default in PhenoGraph)
    connectivities = (connectivities + connectivities.T) / 2

    return connectivities
