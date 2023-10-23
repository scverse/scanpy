import os
from typing import Literal, overload

import numpy as np
import scipy
import scipy.stats
from scipy import sparse
from numpy.typing import NDArray
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.neighbors import NearestNeighbors

from scanpy._utils import AnyRandom

Scale = Literal["linear", "log", "symlog", "logit"] | str


########## PREPROCESSING PIPELINE


def print_optional(string, verbose: bool = True) -> None:
    if verbose:
        print(string)


def pipeline_normalize(self, postnorm_total: float | None = None):
    """Total counts normalization"""
    if postnorm_total is None:
        postnorm_total = self._total_counts_obs.mean()

    self._E_obs_norm = tot_counts_norm(
        self._E_obs, target_total=postnorm_total, total_counts=self._total_counts_obs
    )

    if self._E_sim is not None:
        self._E_sim_norm = tot_counts_norm(
            self._E_sim,
            target_total=postnorm_total,
            total_counts=self._total_counts_sim,
        )
    return


def pipeline_get_gene_filter(
    self,
    min_counts: int = 3,
    min_cells: int = 3,
    min_gene_variability_pctl: float | int = 85,
) -> None:
    """Identify highly variable genes expressed above a minimum level"""
    self._gene_filter = filter_genes(
        self._E_obs_norm,
        min_counts=min_counts,
        min_cells=min_cells,
        min_vscore_pctl=min_gene_variability_pctl,
    )


def pipeline_apply_gene_filter(self) -> None:
    if self._E_obs is not None:
        self._E_obs = self._E_obs[:, self._gene_filter]
    if self._E_obs_norm is not None:
        self._E_obs_norm = self._E_obs_norm[:, self._gene_filter]
    if self._E_sim is not None:
        self._E_sim = self._E_sim[:, self._gene_filter]
    if self._E_sim_norm is not None:
        self._E_sim_norm = self._E_sim_norm[:, self._gene_filter]


def pipeline_mean_center(self) -> None:
    gene_means = self._E_obs_norm.mean(0)
    self._E_obs_norm = self._E_obs_norm - gene_means
    if self._E_sim_norm is not None:
        self._E_sim_norm = self._E_sim_norm - gene_means


def pipeline_normalize_variance(self) -> None:
    gene_stdevs = np.sqrt(sparse_var(self._E_obs_norm))
    self._E_obs_norm = sparse_multiply(self._E_obs_norm.T, 1 / gene_stdevs).T
    if self._E_sim_norm is not None:
        self._E_sim_norm = sparse_multiply(self._E_sim_norm.T, 1 / gene_stdevs).T


def pipeline_zscore(self) -> None:
    gene_means = self._E_obs_norm.mean(0)
    gene_stdevs = np.sqrt(sparse_var(self._E_obs_norm))
    self._E_obs_norm = np.array(
        sparse_zscore(self._E_obs_norm, gene_means, gene_stdevs)
    )
    if self._E_sim_norm is not None:
        self._E_sim_norm = np.array(
            sparse_zscore(self._E_sim_norm, gene_means, gene_stdevs)
        )


def pipeline_log_transform(self, pseudocount: int = 1) -> None:
    self._E_obs_norm = log_normalize(self._E_obs_norm, pseudocount)
    if self._E_sim_norm is not None:
        self._E_sim_norm = log_normalize(self._E_sim_norm, pseudocount)


def pipeline_truncated_svd(
    self,
    n_prin_comps: int = 30,
    *,
    random_state: AnyRandom = 0,
    algorithm: Literal['arpack', 'randomized'] = 'arpack',
):
    svd = TruncatedSVD(
        n_components=n_prin_comps, random_state=random_state, algorithm=algorithm
    ).fit(self._E_obs_norm)
    self.set_manifold(svd.transform(self._E_obs_norm), svd.transform(self._E_sim_norm))
    return


def pipeline_pca(
    self,
    n_prin_comps: int = 50,
    *,
    random_state: AnyRandom = 0,
    svd_solver: Literal['auto', 'full', 'arpack', 'randomized'] = 'arpack',
):
    if sparse.issparse(self._E_obs_norm):
        X_obs = self._E_obs_norm.toarray()
    else:
        X_obs = self._E_obs_norm
    if sparse.issparse(self._E_sim_norm):
        X_sim = self._E_sim_norm.toarray()
    else:
        X_sim = self._E_sim_norm

    pca = PCA(
        n_components=n_prin_comps, random_state=random_state, svd_solver=svd_solver
    ).fit(X_obs)
    self.set_manifold(pca.transform(X_obs), pca.transform(X_sim))
    return


def matrix_multiply(X: sparse.spmatrix | NDArray, Y: sparse.spmatrix | NDArray):
    if not type(X) == np.ndarray:
        if sparse.issparse(X):
            X = X.toarray()
        else:
            X = np.array(X)
    if not type(Y) == np.ndarray:
        if sparse.issparse(Y):
            Y = Y.toarray()
        else:
            Y = np.array(Y)
    return np.dot(X, Y)


def log_normalize(X, pseudocount=1):
    X.data = np.log10(X.data + pseudocount)
    return X


########## LOADING DATA
def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []
    gene_dict = {}

    with open(filename) as f:
        for iL in range(skip_rows):
            f.readline()
        for l in f:
            gene = l.strip('\n').split(delimiter)[column]
            if gene in gene_dict:
                gene_dict[gene] += 1
                gene_list.append(gene + '__' + str(gene_dict[gene]))
                if gene_dict[gene] == 2:
                    i = gene_list.index(gene)
                    gene_list[i] = gene + '__1'
            else:
                gene_dict[gene] = 1
                gene_list.append(gene)
    return gene_list


def make_genes_unique(orig_gene_list):
    gene_list = []
    gene_dict = {}

    for gene in orig_gene_list:
        if gene in gene_dict:
            gene_dict[gene] += 1
            gene_list.append(gene + '__' + str(gene_dict[gene]))
            if gene_dict[gene] == 2:
                i = gene_list.index(gene)
                gene_list[i] = gene + '__1'
        else:
            gene_dict[gene] = 1
            gene_list.append(gene)
    return gene_list


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


def sparse_multiply(E: sparse.csr_matrix | sparse.csc_matrix, a):
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
):
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


def filter_genes(
    E,
    base_ix=[],
    min_vscore_pctl: float | int = 85,
    min_counts: int = 3,
    min_cells: int = 3,
    show_vscore_plot: bool = False,
    sample_name: str | None = None,
):
    """
    Filter genes by expression level and variability
    Return list of filtered gene indices
    """

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(
        E[base_ix, :]
    )
    ix2 = Vscores > 0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = ((E[:, gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (
        Vscores >= min_vscore
    )

    if show_vscore_plot:
        import matplotlib.pyplot as plt

        x_min = 0.5 * np.min(mu_gene)
        x_max = 2 * np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max / x_min) * np.linspace(0, 1, 100))
        yTh = (1 + a) * (1 + b) + b * xTh
        plt.figure(figsize=(8, 6))
        plt.scatter(
            np.log10(mu_gene),
            np.log10(FF_gene),
            c=(0.8, 0.8, 0.8),
            alpha=0.3,
            edgecolors=None,
        )
        plt.scatter(
            np.log10(mu_gene)[ix],
            np.log10(FF_gene)[ix],
            c=(0, 0, 0),
            alpha=0.3,
            edgecolors=None,
        )
        plt.plot(np.log10(xTh), np.log10(yTh))
        if sample_name is not None:
            plt.title(sample_name)
        plt.xlabel('log10(mean)')
        plt.ylabel('log10(Fano factor)')
        plt.show()

    return gene_ix[ix]


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


########## DIMENSIONALITY REDUCTION


def get_pca(
    E,
    base_ix=[],
    numpc=50,
    keep_sparse=False,
    normalize=True,
    random_state=0,
    svd_solver='arpack',
):
    """
    Run PCA on the counts matrix E, gene-level normalizing if desired
    Return PCA coordinates
    """
    # If keep_sparse is True, gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    if keep_sparse:
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix, :]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(
            n_components=numpc, random_state=random_state, algorithm=svd_solver
        )

    else:
        if normalize:
            zmean = E[base_ix, :].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix, :]))
            Z = sparse_multiply((E - zmean).T, 1 / zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc, random_state=random_state, svd_solver=svd_solver)

    pca.fit(Z[base_ix, :])
    return pca.transform(Z)


def preprocess_and_pca(
    E,
    *,
    total_counts_normalize: bool = True,
    norm_exclude_abundant_gene_frac: float = 1.0,
    min_counts: int = 3,
    min_cells: int = 5,
    min_vscore_pctl: float | int = 85,
    gene_filter: NDArray[np.intp] | None = None,
    num_pc: int = 50,
    sparse_pca: bool = False,
    zscore_normalize: bool = True,
    show_vscore_plot: bool = False,
):
    """
    Total counts normalize, filter genes, run PCA
    Return PCA coordinates and filtered gene indices
    """

    if total_counts_normalize:
        print('Total count normalizing')
        E = tot_counts_norm(E, exclude_dominant_frac=norm_exclude_abundant_gene_frac)[0]

    if gene_filter is None:
        print('Finding highly variable genes')
        gene_filter = filter_genes(
            E,
            min_vscore_pctl=min_vscore_pctl,
            min_counts=min_counts,
            min_cells=min_cells,
            show_vscore_plot=show_vscore_plot,
        )

    print('Using %i genes for PCA' % len(gene_filter))
    PCdat = get_pca(
        E[:, gene_filter],
        numpc=num_pc,
        keep_sparse=sparse_pca,
        normalize=zscore_normalize,
    )

    return PCdat, gene_filter


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


def build_adj_mat(edges, n_nodes: int) -> sparse.csc_matrix:
    A = sparse.lil_matrix((n_nodes, n_nodes))
    for e in edges:
        i, j = e
        A[i, j] = 1
        A[j, i] = 1
    return A.tocsc()


########## 2-D EMBEDDINGS


def get_umap(
    X,
    *,
    n_neighbors: int = 10,
    min_dist: float = 0.1,
    metric: str = 'euclidean',
    random_state: AnyRandom = 0,
):
    import umap

    return umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    ).fit_transform(X)


def get_tsne(
    X,
    *,
    angle: float = 0.5,
    perplexity: int = 30,
    random_state: AnyRandom = 0,
    verbose: bool = False,
):
    from sklearn.manifold import TSNE

    return TSNE(
        angle=angle, perplexity=perplexity, random_state=random_state, verbose=verbose
    ).fit_transform(X)


def get_force_layout(
    X,
    *,
    n_neighbors: int = 5,
    approx_neighbors: bool = False,
    n_iter: int = 300,
    verbose: bool = False,
):
    edges = get_knn_graph(X, k=n_neighbors, approx=approx_neighbors, return_edges=True)[
        0
    ]
    return run_force_layout(edges, X.shape[0], verbose=verbose)


def run_force_layout(
    links,
    n_cells: int,
    *,
    n_iter: int = 100,
    edgeWeightInfluence=1,
    barnesHutTheta=2,
    scalingRatio=1,
    gravity=0.05,
    jitterTolerance=1,
    verbose=False,
):
    from fa2 import ForceAtlas2
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
        # Behavior alternatives
        outboundAttractionDistribution=False,  # Dissuade hubs
        linLogMode=False,  # NOT IMPLEMENTED
        adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
        edgeWeightInfluence=edgeWeightInfluence,
        # Performance
        jitterTolerance=jitterTolerance,  # Tolerance
        barnesHutOptimize=True,
        barnesHutTheta=barnesHutTheta,
        multiThreaded=False,  # NOT IMPLEMENTED
        # Tuning
        scalingRatio=scalingRatio,
        strongGravityMode=False,
        gravity=gravity,
        # Log
        verbose=verbose,
    )

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    return positions


########## CLUSTERING


def get_spectral_clusters(A, k):
    from sklearn.cluster import SpectralClustering

    spec = SpectralClustering(
        n_clusters=k, random_state=0, affinity='precomputed', assign_labels='discretize'
    )
    return spec.fit_predict(A)


def get_louvain_clusters(nodes, edges):
    import networkx as nx
    import community

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    return np.array(list(community.best_partition(G).values()))


########## GENE ENRICHMENT


def rank_enriched_genes(
    E,
    gene_list,
    cell_mask,
    min_counts: int = 3,
    min_cells: int = 3,
    verbose: bool = False,
):
    gix = (E[cell_mask, :] >= min_counts).sum(0).A.squeeze() >= min_cells
    print_optional('%i cells in group' % (sum(cell_mask)), verbose)
    print_optional('Considering %i genes' % (sum(gix)), verbose)

    gene_list = gene_list[gix]

    z = sparse_zscore(E[:, gix])
    scores = z[cell_mask, :].mean(0).A.squeeze()
    o = np.argsort(-scores)

    return gene_list[o], scores[o]
