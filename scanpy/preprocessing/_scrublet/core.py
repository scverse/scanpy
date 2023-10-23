import time
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from .utils import (
    custom_cmap,
    darken_cmap,
    get_knn_graph,
    pipeline_apply_gene_filter,
    pipeline_get_gene_filter,
    pipeline_log_transform,
    pipeline_mean_center,
    pipeline_normalize,
    pipeline_normalize_variance,
    pipeline_pca,
    pipeline_truncated_svd,
    pipeline_zscore,
    print_optional,
    subsample_counts,
)


class Scrublet:
    def __init__(
        self,
        counts_matrix,
        total_counts=None,
        sim_doublet_ratio=2.0,
        n_neighbors=None,
        expected_doublet_rate=0.1,
        stdev_doublet_rate=0.02,
        random_state=0,
    ):
        '''Initialize Scrublet object with counts matrix and doublet prediction parameters

        Parameters
        ----------
        counts_matrix : scipy sparse matrix or ndarray, shape (n_cells, n_genes)
            Matrix containing raw (unnormalized) UMI-based transcript counts.
            Converted into a scipy.sparse.csc_matrix.

        total_counts : ndarray, shape (n_cells,), optional (default: None)
            Array of total UMI counts per cell. If `None`, this is calculated
            as the row sums of `counts_matrix`.

        sim_doublet_ratio : float, optional (default: 2.0)
            Number of doublets to simulate relative to the number of observed
            transcriptomes.

        n_neighbors : int, optional (default: None)
            Number of neighbors used to construct the KNN graph of observed
            transcriptomes and simulated doublets. If `None`, this is
            set to round(0.5 * sqrt(n_cells))

        expected_doublet_rate : float, optional (default: 0.1)
            The estimated doublet rate for the experiment.

        stdev_doublet_rate : float, optional (default: 0.02)
            Uncertainty in the expected doublet rate.

        random_state : int, optional (default: 0)
            Random state for doublet simulation, approximate
            nearest neighbor search, and PCA/TruncatedSVD.

        Attributes
        ----------
        predicted_doublets_ : ndarray, shape (n_cells,)
            Boolean mask of predicted doublets in the observed
            transcriptomes.

        doublet_scores_obs_ : ndarray, shape (n_cells,)
            Doublet scores for observed transcriptomes.

        doublet_scores_sim_ : ndarray, shape (n_doublets,)
            Doublet scores for simulated doublets.

        doublet_errors_obs_ : ndarray, shape (n_cells,)
            Standard error in the doublet scores for observed
            transcriptomes.

        doublet_errors_sim_ : ndarray, shape (n_doublets,)
            Standard error in the doublet scores for simulated
            doublets.

        threshold_: float
            Doublet score threshold for calling a transcriptome
            a doublet.

        z_scores_ : ndarray, shape (n_cells,)
            Z-score conveying confidence in doublet calls.
            Z = `(doublet_score_obs_ - threhsold_) / doublet_errors_obs_`

        detected_doublet_rate_: float
            Fraction of observed transcriptomes that have been called
            doublets.

        detectable_doublet_fraction_: float
            Estimated fraction of doublets that are detectable, i.e.,
            fraction of simulated doublets with doublet scores above
            `threshold_`

        overall_doublet_rate_: float
            Estimated overall doublet rate,
            `detected_doublet_rate_ / detectable_doublet_fraction_`.
            Should agree (roughly) with `expected_doublet_rate`.

        manifold_obs_: ndarray, shape (n_cells, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for observed transcriptomes. Nearest neighbors are found using
            the union of `manifold_obs_` and `manifold_sim_` (see below).

        manifold_sim_: ndarray, shape (n_doublets, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for simulated doublets. Nearest neighbors are found using
            the union of `manifold_obs_` (see above) and `manifold_sim_`.

        doublet_parents_ : ndarray, shape (n_doublets, 2)
            Indices of the observed transcriptomes used to generate the
            simulated doublets.

        doublet_neighbor_parents_ : list, length n_cells
            A list of arrays of the indices of the doublet neighbors of
            each observed transcriptome (the ith entry is an array of
            the doublet neighbors of transcriptome i).
        '''

        if not scipy.sparse.issparse(counts_matrix):
            counts_matrix = scipy.sparse.csc_matrix(counts_matrix)
        elif not scipy.sparse.isspmatrix_csc(counts_matrix):
            counts_matrix = counts_matrix.tocsc()

        # initialize counts matrices
        self._E_obs = counts_matrix
        self._E_sim = None
        self._E_obs_norm = None
        self._E_sim_norm = None

        if total_counts is None:
            self._total_counts_obs = self._E_obs.sum(1).A.squeeze()
        else:
            self._total_counts_obs = total_counts

        self._gene_filter = np.arange(self._E_obs.shape[1])
        self._embeddings = {}

        self.sim_doublet_ratio = sim_doublet_ratio
        self.n_neighbors = n_neighbors
        self.expected_doublet_rate = expected_doublet_rate
        self.stdev_doublet_rate = stdev_doublet_rate
        self.random_state = random_state

        if self.n_neighbors is None:
            self.n_neighbors = int(round(0.5 * np.sqrt(self._E_obs.shape[0])))

    ######## Core Scrublet functions ########

    def scrub_doublets(
        self,
        synthetic_doublet_umi_subsampling=1.0,
        use_approx_neighbors=True,
        distance_metric='euclidean',
        get_doublet_neighbor_parents=False,
        min_counts=3,
        min_cells=3,
        min_gene_variability_pctl=85,
        log_transform=False,
        mean_center=True,
        normalize_variance=True,
        n_prin_comps=30,
        svd_solver='arpack',
        verbose=True,
    ):
        '''Standard pipeline for preprocessing, doublet simulation, and doublet prediction

        Automatically sets a threshold for calling doublets, but it's best to check
        this by running plot_histogram() afterwards and adjusting threshold
        with call_doublets(threshold=new_threshold) if necessary.

        Arguments
        ---------
        synthetic_doublet_umi_subsampling : float, optional (defuault: 1.0)
            Rate for sampling UMIs when creating synthetic doublets. If 1.0,
            each doublet is created by simply adding the UMIs from two randomly
            sampled observed transcriptomes. For values less than 1, the
            UMI counts are added and then randomly sampled at the specified
            rate.

        use_approx_neighbors : bool, optional (default: True)
            Use approximate nearest neighbor method (annoy) for the KNN
            classifier.

        distance_metric : str, optional (default: 'euclidean')
            Distance metric used when finding nearest neighbors. For list of
            valid values, see the documentation for annoy (if `use_approx_neighbors`
            is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
            is False).

        get_doublet_neighbor_parents : bool, optional (default: False)
            If True, return the parent transcriptomes that generated the
            doublet neighbors of each observed transcriptome. This information can
            be used to infer the cell states that generated a given
            doublet state.

        min_counts : float, optional (default: 3)
            Used for gene filtering prior to PCA. Genes expressed at fewer than
            `min_counts` in fewer than `min_cells` (see below) are excluded.

        min_cells : int, optional (default: 3)
            Used for gene filtering prior to PCA. Genes expressed at fewer than
            `min_counts` (see above) in fewer than `min_cells` are excluded.

        min_gene_variability_pctl : float, optional (default: 85.0)
            Used for gene filtering prior to PCA. Keep the most highly variable genes
            (in the top min_gene_variability_pctl percentile), as measured by
            the v-statistic [Klein et al., Cell 2015].

        log_transform : bool, optional (default: False)
            If True, log-transform the counts matrix (log10(1+TPM)).
            `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
            reduction, unless `mean_center` is True.

        mean_center : bool, optional (default: True)
            If True, center the data such that each gene has a mean of 0.
            `sklearn.decomposition.PCA` will be used for dimensionality
            reduction.

        normalize_variance : bool, optional (default: True)
            If True, normalize the data such that each gene has a variance of 1.
            `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
            reduction, unless `mean_center` is True.

        n_prin_comps : int, optional (default: 30)
            Number of principal components used to embed the transcriptomes prior
            to k-nearest-neighbor graph construction.

        svd_solver : str, optional (default: 'arpack')
            SVD solver to use. See available options for
            `svd_solver` from `sklearn.decomposition.PCA` or
            `algorithm` from `sklearn.decomposition.TruncatedSVD`

        verbose : bool, optional (default: True)
            If True, print progress updates.

        Sets
        ----
        doublet_scores_obs_, doublet_errors_obs_,
        doublet_scores_sim_, doublet_errors_sim_,
        predicted_doublets_, z_scores_
        threshold_, detected_doublet_rate_,
        detectable_doublet_fraction_, overall_doublet_rate_,
        doublet_parents_, doublet_neighbor_parents_
        '''
        t0 = time.time()

        self._E_sim = None
        self._E_obs_norm = None
        self._E_sim_norm = None
        self._gene_filter = np.arange(self._E_obs.shape[1])

        print_optional('Preprocessing...', verbose)
        pipeline_normalize(self)
        pipeline_get_gene_filter(
            self,
            min_counts=min_counts,
            min_cells=min_cells,
            min_gene_variability_pctl=min_gene_variability_pctl,
        )
        pipeline_apply_gene_filter(self)

        print_optional('Simulating doublets...', verbose)
        self.simulate_doublets(
            sim_doublet_ratio=self.sim_doublet_ratio,
            synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling,
        )
        pipeline_normalize(self, postnorm_total=1e6)
        if log_transform:
            pipeline_log_transform(self)
        if mean_center and normalize_variance:
            pipeline_zscore(self)
        elif mean_center:
            pipeline_mean_center(self)
        elif normalize_variance:
            pipeline_normalize_variance(self)

        if mean_center:
            print_optional('Embedding transcriptomes using PCA...', verbose)
            pipeline_pca(
                self,
                n_prin_comps=n_prin_comps,
                random_state=self.random_state,
                svd_solver=svd_solver,
            )
        else:
            print_optional('Embedding transcriptomes using Truncated SVD...', verbose)
            pipeline_truncated_svd(
                self,
                n_prin_comps=n_prin_comps,
                random_state=self.random_state,
                algorithm=svd_solver,
            )

        print_optional('Calculating doublet scores...', verbose)
        self.calculate_doublet_scores(
            use_approx_neighbors=use_approx_neighbors,
            distance_metric=distance_metric,
            get_doublet_neighbor_parents=get_doublet_neighbor_parents,
        )
        self.call_doublets(verbose=verbose)

        t1 = time.time()
        print_optional('Elapsed time: {:.1f} seconds'.format(t1 - t0), verbose)
        return self.doublet_scores_obs_, self.predicted_doublets_

    def simulate_doublets(
        self, sim_doublet_ratio=None, synthetic_doublet_umi_subsampling=1.0
    ):
        '''Simulate doublets by adding the counts of random observed transcriptome pairs.

        Arguments
        ---------
        sim_doublet_ratio : float, optional (default: None)
            Number of doublets to simulate relative to the number of observed
            transcriptomes. If `None`, self.sim_doublet_ratio is used.

        synthetic_doublet_umi_subsampling : float, optional (defuault: 1.0)
            Rate for sampling UMIs when creating synthetic doublets. If 1.0,
            each doublet is created by simply adding the UMIs from two randomly
            sampled observed transcriptomes. For values less than 1, the
            UMI counts are added and then randomly sampled at the specified
            rate.

        Sets
        ----
        doublet_parents_
        '''

        if sim_doublet_ratio is None:
            sim_doublet_ratio = self.sim_doublet_ratio
        else:
            self.sim_doublet_ratio = sim_doublet_ratio

        n_obs = self._E_obs.shape[0]
        n_sim = int(n_obs * sim_doublet_ratio)

        np.random.seed(self.random_state)
        pair_ix = np.random.randint(0, n_obs, size=(n_sim, 2))

        E1 = self._E_obs[pair_ix[:, 0], :]
        E2 = self._E_obs[pair_ix[:, 1], :]
        tots1 = self._total_counts_obs[pair_ix[:, 0]]
        tots2 = self._total_counts_obs[pair_ix[:, 1]]
        if synthetic_doublet_umi_subsampling < 1:
            self._E_sim, self._total_counts_sim = subsample_counts(
                E1 + E2,
                synthetic_doublet_umi_subsampling,
                tots1 + tots2,
                random_seed=self.random_state,
            )
        else:
            self._E_sim = E1 + E2
            self._total_counts_sim = tots1 + tots2
        self.doublet_parents_ = pair_ix
        return

    def set_manifold(self, manifold_obs, manifold_sim):
        '''Set the manifold coordinates used in k-nearest-neighbor graph construction

        Arguments
        ---------
        manifold_obs: ndarray, shape (n_cells, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for observed transcriptomes. Nearest neighbors are found using
            the union of `manifold_obs` and `manifold_sim` (see below).

        manifold_sim: ndarray, shape (n_doublets, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for simulated doublets. Nearest neighbors are found using
            the union of `manifold_obs` (see above) and `manifold_sim`.

        Sets
        ----
        manifold_obs_, manifold_sim_,
        '''

        self.manifold_obs_ = manifold_obs
        self.manifold_sim_ = manifold_sim
        return

    def calculate_doublet_scores(
        self,
        use_approx_neighbors=True,
        distance_metric='euclidean',
        get_doublet_neighbor_parents=False,
    ):
        '''Calculate doublet scores for observed transcriptomes and simulated doublets

        Requires that manifold_obs_ and manifold_sim_ have already been set.

        Arguments
        ---------
        use_approx_neighbors : bool, optional (default: True)
            Use approximate nearest neighbor method (annoy) for the KNN
            classifier.

        distance_metric : str, optional (default: 'euclidean')
            Distance metric used when finding nearest neighbors. For list of
            valid values, see the documentation for annoy (if `use_approx_neighbors`
            is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
            is False).

        get_doublet_neighbor_parents : bool, optional (default: False)
            If True, return the parent transcriptomes that generated the
            doublet neighbors of each observed transcriptome. This information can
            be used to infer the cell states that generated a given
            doublet state.

        Sets
        ----
        doublet_scores_obs_, doublet_scores_sim_,
        doublet_errors_obs_, doublet_errors_sim_,
        doublet_neighbor_parents_

        '''

        self._nearest_neighbor_classifier(
            k=self.n_neighbors,
            exp_doub_rate=self.expected_doublet_rate,
            stdev_doub_rate=self.stdev_doublet_rate,
            use_approx_nn=use_approx_neighbors,
            distance_metric=distance_metric,
            get_neighbor_parents=get_doublet_neighbor_parents,
        )
        return self.doublet_scores_obs_

    def _nearest_neighbor_classifier(
        self,
        k=40,
        use_approx_nn=True,
        distance_metric='euclidean',
        exp_doub_rate=0.1,
        stdev_doub_rate=0.03,
        get_neighbor_parents=False,
    ):
        manifold = np.vstack((self.manifold_obs_, self.manifold_sim_))
        doub_labels = np.concatenate(
            (
                np.zeros(self.manifold_obs_.shape[0], dtype=int),
                np.ones(self.manifold_sim_.shape[0], dtype=int),
            )
        )

        n_obs = np.sum(doub_labels == 0)
        n_sim = np.sum(doub_labels == 1)

        # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
        k_adj = int(round(k * (1 + n_sim / float(n_obs))))

        # Find k_adj nearest neighbors
        neighbors = get_knn_graph(
            manifold,
            k=k_adj,
            dist_metric=distance_metric,
            approx=use_approx_nn,
            return_edges=False,
            random_seed=self.random_state,
        )

        # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
        doub_neigh_mask = doub_labels[neighbors] == 1
        n_sim_neigh = doub_neigh_mask.sum(1)

        rho = exp_doub_rate
        r = n_sim / float(n_obs)
        nd = n_sim_neigh.astype(float)
        N = float(k_adj)

        # Bayesian
        q = (nd + 1) / (N + 2)
        Ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r))

        se_q = np.sqrt(q * (1 - q) / (N + 3))
        se_rho = stdev_doub_rate

        se_Ld = (
            q
            * rho
            / r
            / (1 - rho - q * (1 - rho - rho / r)) ** 2
            * np.sqrt((se_q / q * (1 - rho)) ** 2 + (se_rho / rho * (1 - q)) ** 2)
        )

        self.doublet_scores_obs_ = Ld[doub_labels == 0]
        self.doublet_scores_sim_ = Ld[doub_labels == 1]
        self.doublet_errors_obs_ = se_Ld[doub_labels == 0]
        self.doublet_errors_sim_ = se_Ld[doub_labels == 1]

        # get parents of doublet neighbors, if requested
        neighbor_parents = None
        if get_neighbor_parents:
            parent_cells = self.doublet_parents_
            neighbors = neighbors - n_obs
            neighbor_parents = []
            for iCell in range(n_obs):
                this_doub_neigh = neighbors[iCell, :][neighbors[iCell, :] > -1]
                if len(this_doub_neigh) > 0:
                    this_doub_neigh_parents = np.unique(
                        parent_cells[this_doub_neigh, :].flatten()
                    )
                    neighbor_parents.append(this_doub_neigh_parents)
                else:
                    neighbor_parents.append([])
            self.doublet_neighbor_parents_ = np.array(neighbor_parents)
        return

    def call_doublets(self, threshold=None, verbose=True):
        '''Call trancriptomes as doublets or singlets

        Arguments
        ---------
        threshold : float, optional (default: None)
            Doublet score threshold for calling a transcriptome
            a doublet. If `None`, this is set automatically by looking
            for the minimum between the two modes of the `doublet_scores_sim_`
            histogram. It is best practice to check the threshold visually
            using the `doublet_scores_sim_` histogram and/or based on
            co-localization of predicted doublets in a 2-D embedding.

        verbose : bool, optional (default: True)
            If True, print summary statistics.

        Sets
        ----
        predicted_doublets_, z_scores_, threshold_,
        detected_doublet_rate_, detectable_doublet_fraction,
        overall_doublet_rate_
        '''

        if threshold is None:
            # automatic threshold detection
            # http://scikit-image.org/docs/dev/api/skimage.filters.html
            from skimage.filters import threshold_minimum

            try:
                threshold = threshold_minimum(self.doublet_scores_sim_)
                if verbose:
                    print(
                        "Automatically set threshold at doublet score = {:.2f}".format(
                            threshold
                        )
                    )
            except Exception:
                self.predicted_doublets_ = None
                if verbose:
                    print(
                        "Warning: failed to automatically identify doublet score threshold. Run `call_doublets` with user-specified threshold."
                    )
                return self.predicted_doublets_

        Ld_obs = self.doublet_scores_obs_
        Ld_sim = self.doublet_scores_sim_
        se_obs = self.doublet_errors_obs_
        Z = (Ld_obs - threshold) / se_obs
        self.predicted_doublets_ = Ld_obs > threshold
        self.z_scores_ = Z
        self.threshold_ = threshold
        self.detected_doublet_rate_ = (Ld_obs > threshold).sum() / float(len(Ld_obs))
        self.detectable_doublet_fraction_ = (Ld_sim > threshold).sum() / float(
            len(Ld_sim)
        )
        self.overall_doublet_rate_ = (
            self.detected_doublet_rate_ / self.detectable_doublet_fraction_
        )

        if verbose:
            print(
                'Detected doublet rate = {:.1f}%'.format(
                    100 * self.detected_doublet_rate_
                )
            )
            print(
                'Estimated detectable doublet fraction = {:.1f}%'.format(
                    100 * self.detectable_doublet_fraction_
                )
            )
            print('Overall doublet rate:')
            print('\tExpected   = {:.1f}%'.format(100 * self.expected_doublet_rate))
            print('\tEstimated  = {:.1f}%'.format(100 * self.overall_doublet_rate_))

        return self.predicted_doublets_

    ######## Viz functions ########

    def plot_histogram(
        self, scale_hist_obs='log', scale_hist_sim='linear', fig_size=(8, 3)
    ):
        '''Plot histogram of doublet scores for observed transcriptomes and simulated doublets

        The histogram for simulated doublets is useful for determining the correct doublet
        score threshold. To set threshold to a new value, T, run call_doublets(threshold=T).

        '''

        fig, axs = plt.subplots(1, 2, figsize=fig_size)

        ax = axs[0]
        ax.hist(
            self.doublet_scores_obs_,
            np.linspace(0, 1, 50),
            color='gray',
            linewidth=0,
            density=True,
        )
        ax.set_yscale(scale_hist_obs)
        yl = ax.get_ylim()
        ax.set_ylim(yl)
        ax.plot(self.threshold_ * np.ones(2), yl, c='black', linewidth=1)
        ax.set_title('Observed transcriptomes')
        ax.set_xlabel('Doublet score')
        ax.set_ylabel('Prob. density')

        ax = axs[1]
        ax.hist(
            self.doublet_scores_sim_,
            np.linspace(0, 1, 50),
            color='gray',
            linewidth=0,
            density=True,
        )
        ax.set_yscale(scale_hist_sim)
        yl = ax.get_ylim()
        ax.set_ylim(yl)
        ax.plot(self.threshold_ * np.ones(2), yl, c='black', linewidth=1)
        ax.set_title('Simulated doublets')
        ax.set_xlabel('Doublet score')
        ax.set_ylabel('Prob. density')

        fig.tight_layout()

        return fig, axs

    def set_embedding(self, embedding_name, coordinates):
        '''Add a 2-D embedding for the observed transcriptomes'''
        self._embeddings[embedding_name] = coordinates
        return

    def plot_embedding(
        self,
        embedding_name,
        score='raw',
        marker_size=5,
        order_points=False,
        fig_size=(8, 4),
        color_map=None,
    ):
        '''Plot doublet predictions on 2-D embedding of observed transcriptomes'''

        # from matplotlib.lines import Line2D
        if embedding_name not in self._embeddings:
            print(
                'Cannot find "{}" in embeddings. First add the embedding using `set_embedding`.'.format(
                    embedding_name
                )
            )
            return

        # TO DO: check if self.predicted_doublets exists; plot raw scores only if it doesn't

        fig, axs = plt.subplots(1, 2, figsize=fig_size)

        x = self._embeddings[embedding_name][:, 0]
        y = self._embeddings[embedding_name][:, 1]
        xl = (x.min() - x.ptp() * 0.05, x.max() + x.ptp() * 0.05)
        yl = (y.min() - y.ptp() * 0.05, y.max() + y.ptp() * 0.05)

        ax = axs[1]
        if score == 'raw':
            color_dat = self.doublet_scores_obs_
            vmin = color_dat.min()
            vmax = color_dat.max()
            if color_map is None:
                cmap_use = darken_cmap(plt.cm.Reds, 0.9)
            else:
                cmap_use = color_map
        elif score == 'zscore':
            color_dat = self.z_scores_
            vmin = -color_dat.max()
            vmax = color_dat.max()
            if color_map is None:
                cmap_use = darken_cmap(plt.cm.RdBu_r, 0.9)
            else:
                cmap_use = color_map
        if order_points:
            o = np.argsort(color_dat)
        else:
            o = np.arange(len(color_dat))
        pp = ax.scatter(
            x[o],
            y[o],
            s=marker_size,
            edgecolors=None,
            c=color_dat[o],
            cmap=cmap_use,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Doublet score')
        ax.set_xlabel(embedding_name + ' 1')
        ax.set_ylabel(embedding_name + ' 2')
        fig.colorbar(pp, ax=ax)

        ax = axs[0]
        called_doubs = self.predicted_doublets_
        ax.scatter(
            x[o],
            y[o],
            s=marker_size,
            edgecolors=None,
            c=called_doubs[o],
            cmap=custom_cmap([[0.7, 0.7, 0.7], [0, 0, 0]]),
        )
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Predicted doublets')
        # singlet_marker = Line2D([], [], color=[.7,.7,.7], marker='o', markersize=5, label='Singlet', linewidth=0)
        # doublet_marker = Line2D([], [], color=[.0,.0,.0], marker='o', markersize=5, label='Doublet', linewidth=0)
        # ax.legend(handles = [singlet_marker, doublet_marker])
        ax.set_xlabel(embedding_name + ' 1')
        ax.set_ylabel(embedding_name + ' 2')

        fig.tight_layout()

        return fig, axs
