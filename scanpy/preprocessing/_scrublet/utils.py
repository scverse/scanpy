from __future__ import annotations

from typing import Literal, Union as _U, overload

import numpy as np
from scipy import sparse
from numpy.typing import NDArray

from ... import logging as logg
from ..._utils import AnyRandom, get_random_state

Scale = _U[Literal["linear", "log", "symlog", "logit"], str]


########## USEFUL SPARSE FUNCTIONS


def sparse_var(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    axis: Literal[0, 1],
) -> NDArray[np.float64]:
    """variance across the specified axis"""

    mean_gene: NDArray[np.float64] = E.mean(axis=axis).A.squeeze()
    tmp: sparse.csc_matrix | sparse.csr_matrix = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene**2


def sparse_multiply(
    E: sparse.csr_matrix | sparse.csc_matrix | NDArray[np.float64],
    a: float | int | NDArray[np.float64],
) -> sparse.csr_matrix | sparse.csc_matrix:
    """multiply each row of E by a scalar"""

    nrow = E.shape[0]
    w = sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    r = w @ E
    if isinstance(r, (np.matrix, np.ndarray)):
        return sparse.csc_matrix(r)
    return r


def sparse_zscore(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    gene_mean: NDArray[np.float64] | None = None,
    gene_stdev: NDArray[np.float64] | None = None,
) -> sparse.csr_matrix | sparse.csc_matrix:
    """z-score normalize each column of E"""

    if gene_mean is None:
        gene_mean = E.mean(0)
    if gene_stdev is None:
        gene_stdev = np.sqrt(sparse_var(E, axis=0))
    return sparse_multiply(np.asarray((E - gene_mean).T), 1 / gene_stdev).T


def subsample_counts(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    rate: float,
    original_totals,
    random_seed: AnyRandom = 0,
) -> tuple[sparse.csr_matrix | sparse.csc_matrix, NDArray[np.int64]]:
    if rate < 1:
        random_seed = get_random_state(random_seed)
        E.data = random_seed.binomial(np.round(E.data).astype(int), rate)
        current_totals = E.sum(1).A.squeeze()
        unsampled_orig_totals = original_totals - current_totals
        unsampled_downsamp_totals = np.random.binomial(
            np.round(unsampled_orig_totals).astype(int), rate
        )
        final_downsamp_totals = current_totals + unsampled_downsamp_totals
    else:
        final_downsamp_totals = original_totals
    return E, final_downsamp_totals


########## GRAPH CONSTRUCTION

AnnoyDist = Literal["angular", "euclidean", "manhattan", "hamming", "dot"]


@overload
def get_knn_graph(
    X,
    k: int = 5,
    *,
    dist_metric: AnnoyDist = "euclidean",
    approx: bool = False,
    return_edges: Literal[True] = True,
    random_seed: AnyRandom = 0,
) -> tuple[set[tuple[int, int]], NDArray[np.int64]]:
    ...


@overload
def get_knn_graph(
    X,
    k: int = 5,
    *,
    dist_metric: AnnoyDist = "euclidean",
    approx: bool = False,
    return_edges: Literal[False],
    random_seed: AnyRandom = 0,
) -> NDArray[np.int64]:
    ...


def get_knn_graph(
    X,
    k: int = 5,
    *,
    dist_metric: AnnoyDist = "euclidean",
    approx: bool = False,
    return_edges: bool = True,
    random_seed: AnyRandom = 0,
):
    """
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    """
    from sklearn.neighbors import NearestNeighbors

    if approx:
        try:
            from annoy import AnnoyIndex
        except ImportError as e:
            raise ValueError(
                'Could not find library "annoy" for approx. nearest neighbor search'
            ) from e
        t0 = logg.info("Using approximate nearest neighbor search")

        if dist_metric == "cosine":
            dist_metric = "angular"
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)
        annoy_index.set_seed(random_seed)

        for i in range(ncell):
            annoy_index.add_item(i, list(X[i, :]))
        annoy_index.build(10)  # 10 trees

        knn = np.array(
            [annoy_index.get_nns_by_item(c, k + 1)[1:] for c in range(ncell)],
            dtype=np.intp,
        )
    else:
        t0 = logg.info("Using sklearn NearestNeighbors")

        if dist_metric == "cosine":
            nbrs = NearestNeighbors(
                n_neighbors=k, metric=dist_metric, algorithm="brute"
            ).fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn: NDArray[np.intp] = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set()
        for i in range(knn.shape[0]):
            for j in knn[i, :]:
                links.add(tuple(sorted((i, j))))
        logg.info("kNN graph built in {time_passed:.3} sec", time=t0)
        return links, knn

    return knn
