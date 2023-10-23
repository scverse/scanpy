from typing import Literal, Union as _U, overload

import numpy as np
import scipy
import scipy.stats
from scipy import sparse
from numpy.typing import NDArray
from sklearn.neighbors import NearestNeighbors

from scanpy._utils import AnyRandom

Scale = _U[Literal["linear", "log", "symlog", "logit"], str]


########## USEFUL SPARSE FUNCTIONS


def sparse_var(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    axis: Literal[0, 1] = 0,
) -> NDArray[np.float64]:
    """variance across the specified axis"""

    mean_gene: NDArray[np.float64] = E.mean(axis=axis).A.squeeze()
    tmp: sparse.csc_matrix | sparse.csr_matrix = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene**2


def sparse_multiply(
    E: sparse.csr_matrix | sparse.csc_matrix, a: float | int | NDArray[np.float64]
) -> sparse.csr_matrix | sparse.csc_matrix:
    """multiply each row of E by a scalar"""

    nrow = E.shape[0]
    w = sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E


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
        gene_stdev = np.sqrt(sparse_var(E))
    return sparse_multiply((E - gene_mean).T, 1 / gene_stdev).T


def subsample_counts(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    rate: float,
    original_totals,
    random_seed: AnyRandom = 0,
) -> tuple[sparse.csr_matrix | sparse.csc_matrix, NDArray[np.int64]]:
    if rate < 1:
        np.random.seed(random_seed)
        E.data = np.random.binomial(np.round(E.data).astype(int), rate)
        current_totals = E.sum(1).A.squeeze()
        unsampled_orig_totals = original_totals - current_totals
        unsampled_downsamp_totals = np.random.binomial(
            np.round(unsampled_orig_totals).astype(int), rate
        )
        final_downsamp_totals = current_totals + unsampled_downsamp_totals
    else:
        final_downsamp_totals = original_totals
    return E, final_downsamp_totals


########## GENE FILTERING


def runningquantile(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    p: float | int,
    nBins: int,
):
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0] + dx / 2, x[-1] - dx / 2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i] - dx / 2) & (x < xOut[i] + dx / 2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i - 1]
            else:
                yOut[i] = np.nan

    return xOut, yOut


def get_vscores(
    E,
    *,
    min_mean: float | int = 0,
    nBins: int = 50,
    fit_percentile: float = 0.1,
    error_wt: float | int = 1,
):
    """
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    """

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:, gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene**2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    def gLog(input):
        return np.log(input[1] * np.exp(-input[0]) + input[2])

    h, b = np.histogram(np.log(FF_gene[mu_gene > 0]), bins=200)
    b = b[:-1] + np.diff(b) / 2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))

    def errFun(b2):
        return np.sum(abs(gLog([x, c, b2]) - y) ** error_wt)

    b0 = 0.1
    b = scipy.optimize.fmin(func=errFun, x0=[b0], disp=False)
    a = c / (1 + b) - 1

    v_scores = FF_gene / ((1 + a) * (1 + b) + b * mu_gene)
    CV_eff = np.sqrt((1 + a) * (1 + b) - 1)
    CV_input = np.sqrt(b)

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b


########## CELL NORMALIZATION


def tot_counts_norm(
    E: sparse.spmatrix,
    total_counts: NDArray[np.intp] | None = None,
    exclude_dominant_frac: float = 1.0,
    included=[],
    target_total: float | None = None,
):
    """
    Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
    Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts
    """

    E = sparse.csc_matrix(E)
    ncell = E.shape[0]
    if total_counts is None:
        if len(included) == 0:
            if exclude_dominant_frac == 1:
                tots_use = E.sum(axis=1)
            else:
                tots = E.sum(axis=1)
                wtmp = sparse.lil_matrix((ncell, ncell))
                wtmp.setdiag(1.0 / tots)
                included = np.asarray(
                    ~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0)
                )[0, :]
                tots_use = E[:, included].sum(axis=1)
                print('Excluded %i genes from normalization' % (np.sum(~included)))
        else:
            tots_use = E[:, included].sum(axis=1)
    else:
        tots_use = total_counts.copy()

    if target_total is None:
        target_total = np.mean(tots_use)

    w = sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_total) / tots_use)
    Enorm = w * E

    return Enorm.tocsc()


########## GRAPH CONSTRUCTION


@overload
def get_knn_graph(
    X,
    k: int = 5,
    *,
    dist_metric: str = 'euclidean',
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
    dist_metric: str = 'euclidean',
    approx: bool = False,
    return_edges: Literal[False],
    random_seed: AnyRandom = 0,
) -> NDArray[np.int64]:
    ...


def get_knn_graph(
    X,
    k: int = 5,
    *,
    dist_metric: str = 'euclidean',
    approx: bool = False,
    return_edges: bool = True,
    random_seed: AnyRandom = 0,
):
    """
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    """

    # t0 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except Exception:
            approx = False
            print('Could not find library "annoy" for approx. nearest neighbor search')
    if approx:
        # print('Using approximate nearest neighbor search')

        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)
        annoy_index.set_seed(random_seed)

        for i in range(ncell):
            annoy_index.add_item(i, list(X[i, :]))
        annoy_index.build(10)  # 10 trees

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

    else:
        # print('Using sklearn NearestNeighbors')

        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(
                n_neighbors=k, metric=dist_metric, algorithm='brute'
            ).fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set()
        for i in range(knn.shape[0]):
            for j in knn[i, :]:
                links.add(tuple(sorted((i, j))))

        # print('kNN graph built in %.3f sec' %(time.time() - t0))

        return links, knn
    return knn
