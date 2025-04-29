from __future__ import annotations

from dataclasses import InitVar, dataclass, field
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from anndata import AnnData, concat
from scipy import sparse
from sklearn.utils import check_random_state

from ... import logging as logg
from ...neighbors import (
    Neighbors,
    _get_indices_distances_from_sparse_matrix,
)
from .._utils import sample_comb
from .sparse_utils import subsample_counts

if TYPE_CHECKING:
    from numpy.random import RandomState
    from numpy.typing import NDArray

    from ..._compat import CSBase, CSCBase
    from ..._utils.random import _LegacyRandom
    from ...neighbors import _Metric, _MetricFn

__all__ = ["Scrublet"]


@dataclass(kw_only=True)
class Scrublet:
    """Initialize Scrublet object with counts matrix and doublet prediction parameters.

    Parameters
    ----------
    counts_obs
        Matrix with shape (n_cells, n_genes) containing raw (unnormalized)
        UMI-based transcript counts.
        Converted into a :class:`scipy.sparse.csc_matrix`.

    total_counts_obs
        Array with shape (n_cells,) of total UMI counts per cell.
        If `None`, this is calculated as the row sums of `counts_obs`.

    sim_doublet_ratio
        Number of doublets to simulate relative to the number of observed
        transcriptomes.

    n_neighbors
        Number of neighbors used to construct the KNN graph of observed
        transcriptomes and simulated doublets.
        If `None`, this is set to round(0.5 * sqrt(n_cells))

    expected_doublet_rate
        The estimated doublet rate for the experiment.

    stdev_doublet_rate
        Uncertainty in the expected doublet rate.

    random_state
        Random state for doublet simulation, approximate
        nearest neighbor search, and PCA/TruncatedSVD.

    """

    # init fields

    counts_obs: InitVar[CSBase | NDArray[np.integer]] = field(kw_only=False)
    total_counts_obs: InitVar[NDArray[np.integer] | None] = None
    sim_doublet_ratio: float = 2.0
    n_neighbors: InitVar[int | None] = None
    expected_doublet_rate: float = 0.1
    stdev_doublet_rate: float = 0.02
    random_state: InitVar[_LegacyRandom] = 0

    # private fields

    _n_neighbors: int = field(init=False, repr=False)
    _random_state: RandomState = field(init=False, repr=False)

    _counts_obs: CSCBase = field(init=False, repr=False)
    _total_counts_obs: NDArray[np.integer] = field(init=False, repr=False)
    _counts_obs_norm: CSBase = field(init=False, repr=False)

    _counts_sim: CSBase = field(init=False, repr=False)
    _total_counts_sim: NDArray[np.integer] = field(init=False, repr=False)
    _counts_sim_norm: CSBase | None = field(default=None, init=False, repr=False)

    # Fields set by methods

    predicted_doublets_: NDArray[np.bool_] | None = field(init=False)
    """(shape: n_cells)
    Boolean mask of predicted doublets in the observed transcriptomes.
    """

    doublet_scores_obs_: NDArray[np.float64] = field(init=False)
    """(shape: n_cells)
    Doublet scores for observed transcriptomes.
    """

    doublet_scores_sim_: NDArray[np.float64] = field(init=False)
    """(shape: n_doublets)
    Doublet scores for simulated doublets.
    """

    doublet_errors_obs_: NDArray[np.float64] = field(init=False)
    """(shape: n_cells)
    Standard error in the doublet scores for observed transcriptomes.
    """

    doublet_errors_sim_: NDArray[np.float64] = field(init=False)
    """(shape: n_doublets)
    Standard error in the doublet scores for simulated doublets.
    """

    threshold_: float = field(init=False)
    """Doublet score threshold for calling a transcriptome a doublet."""

    z_scores_: NDArray[np.float64] = field(init=False)
    """(shape: n_cells)
    Z-score conveying confidence in doublet calls.
    Z = `(doublet_score_obs_ - threhsold_) / doublet_errors_obs_`
    """

    detected_doublet_rate_: float = field(init=False)
    """Fraction of observed transcriptomes that have been called doublets."""

    detectable_doublet_fraction_: float = field(init=False)
    """Estimated fraction of doublets that are detectable, i.e.,
    fraction of simulated doublets with doublet scores above `threshold_`
    """

    overall_doublet_rate_: float = field(init=False)
    """Estimated overall doublet rate,
    `detected_doublet_rate_ / detectable_doublet_fraction_`.
    Should agree (roughly) with `expected_doublet_rate`.
    """

    manifold_obs_: NDArray[np.float64] = field(init=False)
    """(shape: n_cells × n_features)
    The single-cell "manifold" coordinates (e.g., PCA coordinates)
    for observed transcriptomes. Nearest neighbors are found using
    the union of `manifold_obs_` and `manifold_sim_` (see below).
    """

    manifold_sim_: NDArray[np.float64] = field(init=False)
    """shape (n_doublets × n_features)
    The single-cell "manifold" coordinates (e.g., PCA coordinates)
    for simulated doublets. Nearest neighbors are found using
    the union of `manifold_obs_` (see above) and `manifold_sim_`.
    """

    doublet_parents_: NDArray[np.intp] = field(init=False)
    """(shape: n_doublets × 2)
    Indices of the observed transcriptomes used to generate the
    simulated doublets.
    """

    doublet_neighbor_parents_: list[NDArray[np.intp]] = field(init=False)
    """(length: n_cells)
    A list of arrays of the indices of the doublet neighbors of
    each observed transcriptome (the ith entry is an array of
    the doublet neighbors of transcriptome i).
    """

    def __post_init__(
        self,
        counts_obs: CSBase | NDArray[np.integer],
        total_counts_obs: NDArray[np.integer] | None,
        n_neighbors: int | None,
        random_state: _LegacyRandom,
    ) -> None:
        self._counts_obs = sparse.csc_matrix(counts_obs)  # noqa: TID251
        self._total_counts_obs = (
            np.asarray(self._counts_obs.sum(1)).squeeze()
            if total_counts_obs is None
            else total_counts_obs
        )
        self._n_neighbors = (
            round(0.5 * np.sqrt(self._counts_obs.shape[0]))
            if n_neighbors is None
            else n_neighbors
        )
        self._random_state = check_random_state(random_state)

    def simulate_doublets(
        self,
        *,
        sim_doublet_ratio: float | None = None,
        synthetic_doublet_umi_subsampling: float = 1.0,
    ) -> None:
        """Simulate doublets by adding the counts of random observed transcriptome pairs.

        Parameters
        ----------
        sim_doublet_ratio
            Number of doublets to simulate relative to the number of observed
            transcriptomes. If `None`, self.sim_doublet_ratio is used.

        synthetic_doublet_umi_subsampling
            Rate for sampling UMIs when creating synthetic doublets.
            If 1.0, each doublet is created by simply adding the UMIs from two randomly
            sampled observed transcriptomes.
            For values less than 1, the UMI counts are added and then randomly sampled
            at the specified rate.

        Sets
        ----
        doublet_parents_

        """
        if sim_doublet_ratio is None:
            sim_doublet_ratio = self.sim_doublet_ratio
        else:
            self.sim_doublet_ratio = sim_doublet_ratio

        n_obs = self._counts_obs.shape[0]
        n_sim = int(n_obs * sim_doublet_ratio)

        pair_ix = sample_comb((n_obs, n_obs), n_sim, random_state=self._random_state)

        E1 = cast("CSCBase", self._counts_obs[pair_ix[:, 0], :])
        E2 = cast("CSCBase", self._counts_obs[pair_ix[:, 1], :])
        tots1 = self._total_counts_obs[pair_ix[:, 0]]
        tots2 = self._total_counts_obs[pair_ix[:, 1]]
        if synthetic_doublet_umi_subsampling < 1:
            self._counts_sim, self._total_counts_sim = subsample_counts(
                E1 + E2,
                rate=synthetic_doublet_umi_subsampling,
                original_totals=tots1 + tots2,
                random_seed=self._random_state,
            )
        else:
            self._counts_sim = E1 + E2
            self._total_counts_sim = tots1 + tots2
        self.doublet_parents_ = pair_ix

    def set_manifold(
        self, manifold_obs: NDArray[np.float64], manifold_sim: NDArray[np.float64]
    ) -> None:
        """Set the manifold coordinates used in k-nearest-neighbor graph construction.

        Parameters
        ----------
        manifold_obs
            (shape: n_cells × n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for observed transcriptomes. Nearest neighbors are found using
            the union of `manifold_obs` and `manifold_sim` (see below).

        manifold_sim
            (shape: n_doublets × n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for simulated doublets. Nearest neighbors are found using
            the union of `manifold_obs` (see above) and `manifold_sim`.

        Sets
        ----
        manifold_obs_, manifold_sim_,

        """
        self.manifold_obs_ = manifold_obs
        self.manifold_sim_ = manifold_sim

    def calculate_doublet_scores(
        self,
        *,
        use_approx_neighbors: bool | None = None,
        distance_metric: _Metric | _MetricFn = "euclidean",
        get_doublet_neighbor_parents: bool = False,
    ) -> NDArray[np.float64]:
        """Calculate doublet scores for observed transcriptomes and simulated doublets.

        Requires that manifold_obs_ and manifold_sim_ have already been set.

        Parameters
        ----------
        use_approx_neighbors
            Use approximate nearest neighbor method (annoy) for the KNN
            classifier.

        distance_metric
            Distance metric used when finding nearest neighbors. For list of
            valid values, see the documentation for annoy (if `use_approx_neighbors`
            is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
            is False).

        get_doublet_neighbor_parents
            If True, return the parent transcriptomes that generated the
            doublet neighbors of each observed transcriptome. This information can
            be used to infer the cell states that generated a given
            doublet state.

        Sets
        ----
        doublet_scores_obs_, doublet_scores_sim_,
        doublet_errors_obs_, doublet_errors_sim_,
        doublet_neighbor_parents_

        """
        self._nearest_neighbor_classifier(
            k=self._n_neighbors,
            exp_doub_rate=self.expected_doublet_rate,
            stdev_doub_rate=self.stdev_doublet_rate,
            use_approx_neighbors=use_approx_neighbors,
            distance_metric=distance_metric,
            get_neighbor_parents=get_doublet_neighbor_parents,
        )
        return self.doublet_scores_obs_

    def _nearest_neighbor_classifier(
        self,
        k: int = 40,
        *,
        use_approx_neighbors: bool | None = None,
        distance_metric: _Metric | _MetricFn = "euclidean",
        exp_doub_rate: float = 0.1,
        stdev_doub_rate: float = 0.03,
        get_neighbor_parents: bool = False,
    ) -> None:
        adatas = [
            AnnData(
                (arr := getattr(self, f"manifold_{n}_")),
                obs=dict(
                    obs_names=pd.RangeIndex(arr.shape[0]).astype("string") + n,
                    doub_labels=n,
                ),
            )
            for n in ["obs", "sim"]
        ]
        manifold = concat(adatas)

        n_obs: int = (manifold.obs["doub_labels"] == "obs").sum()
        n_sim: int = (manifold.obs["doub_labels"] == "sim").sum()

        # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
        k_adj = round(k * (1 + n_sim / float(n_obs)))

        # Find k_adj nearest neighbors
        knn = Neighbors(manifold)
        transformer = None
        if use_approx_neighbors is not None:
            transformer = "pynndescent" if use_approx_neighbors else "sklearn"
        knn.compute_neighbors(
            k_adj,
            metric=distance_metric,
            knn=True,
            transformer=transformer,
            method=None,
            random_state=self._random_state,
        )
        neighbors, _ = _get_indices_distances_from_sparse_matrix(knn.distances, k_adj)
        if use_approx_neighbors:
            neighbors = neighbors[:, 1:]
        # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
        doub_neigh_mask: NDArray[np.bool_] = (
            manifold.obs["doub_labels"].to_numpy()[neighbors] == "sim"
        )
        n_sim_neigh: NDArray[np.int64] = doub_neigh_mask.sum(axis=1)

        rho = exp_doub_rate
        r = n_sim / float(n_obs)
        nd = n_sim_neigh.astype(np.float64)
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

        self.doublet_scores_obs_ = Ld[manifold.obs["doub_labels"] == "obs"]
        self.doublet_scores_sim_ = Ld[manifold.obs["doub_labels"] == "sim"]
        self.doublet_errors_obs_ = se_Ld[manifold.obs["doub_labels"] == "obs"]
        self.doublet_errors_sim_ = se_Ld[manifold.obs["doub_labels"] == "sim"]

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
                    neighbor_parents.append(np.array([], dtype=np.intp))
            self.doublet_neighbor_parents_ = neighbor_parents

    def call_doublets(
        self, *, threshold: float | None = None, verbose: bool = True
    ) -> NDArray[np.bool_] | None:
        """Call trancriptomes as doublets or singlets.

        Parameters
        ----------
        threshold
            Doublet score threshold for calling a transcriptome
            a doublet. If `None`, this is set automatically by looking
            for the minimum between the two modes of the `doublet_scores_sim_`
            histogram. It is best practice to check the threshold visually
            using the `doublet_scores_sim_` histogram and/or based on
            co-localization of predicted doublets in a 2-D embedding.

        verbose
            If True, log summary statistics.

        Sets
        ----
        predicted_doublets_, z_scores_, threshold_,
        detected_doublet_rate_, detectable_doublet_fraction,
        overall_doublet_rate_

        """
        if threshold is None:
            # automatic threshold detection
            # http://scikit-image.org/docs/dev/api/skimage.filters.html
            from skimage.filters import threshold_minimum

            try:
                threshold = cast("float", threshold_minimum(self.doublet_scores_sim_))
                if verbose:
                    logg.info(
                        f"Automatically set threshold at doublet score = {threshold:.2f}"
                    )
            except Exception:  # noqa: BLE001
                self.predicted_doublets_ = None
                if verbose:
                    logg.warning(
                        "Failed to automatically identify doublet score threshold. "
                        "Run `call_doublets` with user-specified threshold."
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
            logg.info(
                f"Detected doublet rate = {100 * self.detected_doublet_rate_:.1f}%\n"
                f"Estimated detectable doublet fraction = {100 * self.detectable_doublet_fraction_:.1f}%\n"
                "Overall doublet rate:\n"
                f"\tExpected   = {100 * self.expected_doublet_rate:.1f}%\n"
                f"\tEstimated  = {100 * self.overall_doublet_rate_:.1f}%"
            )

        return self.predicted_doublets_
