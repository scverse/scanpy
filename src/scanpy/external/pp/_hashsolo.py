"""A probabilistic cell hashing demultiplexing method.

HashSolo generates a noise distribution and signal distribution
for each hashing barcode from empirically observed counts.
These distributions are updates from the global signal and noise barcode distributions,
which helps in the setting where not many cells are observed.
For a hashing barcode:

Signal distributions
    are estimated from samples where that hashing barcode has the highest count.

Noise distributions
    are estimated from samples where that hashing barcode is one the k-2 lowest barcodes,
    where k is the number of barcodes.

We test each of the following hypotheses in a bayesian fashion,
and select the most probable hypothesis.

A doublet
    should have its two highest barcode counts most likely
    coming from a signal distribution for those barcodes.

A singlet
    should have its highest barcode from a signal distribution,
    and its second highest barcode from a noise distribution.

A negative two highest barcodes
    should come from noise distributions.
"""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.stats import norm

from ..._compat import old_positionals
from ..._utils import check_nonnegative_integers
from ..._utils._doctests import doctest_skip

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData
    from numpy.typing import ArrayLike, NDArray


def _calculate_log_likelihoods(  # noqa: PLR0915
    data: np.ndarray, number_of_noise_barcodes: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], dict[int, str]]:
    """Calculate log likelihoods for each hypothesis, negative, singlet, doublet.

    Parameters
    ----------
    data
        cells by hashing counts matrix
    number_of_noise_barcodes
        number of barcodes to used to calculated noise distribution

    Returns
    -------
    log_likelihoods_for_each_hypothesis
        a 2d np.array log likelihood of each hypothesis
    all_indices
    counter_to_barcode_combo

    """

    def gaussian_updates(
        data: np.ndarray, mu_o: float, std_o: float
    ) -> tuple[float, float]:
        """Update parameters of your gaussian.

        See <https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf>.

        Parameters
        ----------
        data
            1-d array of counts
        mu_o
            global mean for hashing count distribution
        std_o
            global std for hashing count distribution

        Returns
        -------
        mean
            of gaussian
        std
            of gaussian

        """
        lam_o = 1 / (std_o**2)
        n = len(data)
        lam = 1 / np.var(data) if len(data) > 1 else lam_o
        lam_n = lam_o + n * lam
        mu_n = (
            (np.mean(data) * n * lam + mu_o * lam_o) / lam_n if len(data) > 0 else mu_o
        )
        return mu_n, (1 / (lam_n / (n + 1))) ** (1 / 2)

    eps = 1e-15
    # probabilites for negative, singlet, doublets
    log_likelihoods_for_each_hypothesis = np.zeros((data.shape[0], 3))

    all_indices = np.empty(data.shape[0])
    num_of_barcodes = data.shape[1]
    number_of_non_noise_barcodes = (
        num_of_barcodes - number_of_noise_barcodes
        if number_of_noise_barcodes is not None
        else 2
    )

    num_of_noise_barcodes = num_of_barcodes - number_of_non_noise_barcodes

    # assume log normal
    data = np.log(data + 1)
    data_arg = np.argsort(data, axis=1)
    data_sort = np.sort(data, axis=1)

    # global signal and noise counts useful for when we have few cells
    # barcodes with the highest number of counts are assumed to be a true signal
    # barcodes with rank < k are considered to be noise
    global_signal_counts = np.ravel(data_sort[:, -1])
    global_noise_counts = np.ravel(data_sort[:, :-number_of_non_noise_barcodes])
    global_mu_signal_o, global_sigma_signal_o = (
        np.mean(global_signal_counts),
        np.std(global_signal_counts),
    )
    global_mu_noise_o, global_sigma_noise_o = (
        np.mean(global_noise_counts),
        np.std(global_noise_counts),
    )

    noise_params_dict = {}
    signal_params_dict = {}

    # for each barcode get  empirical noise and signal distribution parameterization
    for x in np.arange(num_of_barcodes):
        sample_barcodes = data[:, x]
        sample_barcodes_noise_idx = np.where(data_arg[:, :num_of_noise_barcodes] == x)[
            0
        ]
        sample_barcodes_signal_idx = np.where(data_arg[:, -1] == x)

        # get noise and signal counts
        noise_counts = sample_barcodes[sample_barcodes_noise_idx]
        signal_counts = sample_barcodes[sample_barcodes_signal_idx]

        # get parameters of distribution, assuming lognormal do update from global values
        noise_param = gaussian_updates(
            noise_counts, global_mu_noise_o, global_sigma_noise_o
        )
        signal_param = gaussian_updates(
            signal_counts, global_mu_signal_o, global_sigma_signal_o
        )
        noise_params_dict[x] = noise_param
        signal_params_dict[x] = signal_param

    counter_to_barcode_combo: dict[int, str] = {}
    counter = 0

    # for each combination of noise and signal barcode calculate probiltiy of in silico and real cell hypotheses
    for noise_sample_idx, signal_sample_idx in product(
        np.arange(num_of_barcodes), np.arange(num_of_barcodes)
    ):
        signal_subset = data_arg[:, -1] == signal_sample_idx
        noise_subset = data_arg[:, -2] == noise_sample_idx
        subset = signal_subset & noise_subset
        if sum(subset) == 0:
            continue

        indices = np.where(subset)[0]
        barcode_combo = "_".join([str(noise_sample_idx), str(signal_sample_idx)])
        all_indices[np.where(subset)[0]] = counter
        counter_to_barcode_combo[counter] = barcode_combo
        counter += 1
        noise_params = noise_params_dict[noise_sample_idx]
        signal_params = signal_params_dict[signal_sample_idx]

        # calculate probabilties for each hypothesis for each cell
        data_subset = data[subset]
        log_signal_signal_probs = np.log(
            norm.pdf(
                data_subset[:, signal_sample_idx],
                *signal_params[:-2],
                loc=signal_params[-2],
                scale=signal_params[-1],
            )
            + eps
        )
        signal_noise_params = signal_params_dict[noise_sample_idx]
        log_noise_signal_probs = np.log(
            norm.pdf(
                data_subset[:, noise_sample_idx],
                loc=signal_noise_params[-2],
                scale=signal_noise_params[-1],
            )
            + eps
        )

        log_noise_noise_probs = np.log(
            norm.pdf(
                data_subset[:, noise_sample_idx],
                loc=noise_params[-2],
                scale=noise_params[-1],
            )
            + eps
        )
        log_signal_noise_probs = np.log(
            norm.pdf(
                data_subset[:, signal_sample_idx],
                loc=noise_params[-2],
                scale=noise_params[-1],
            )
            + eps
        )

        probs_of_negative = np.sum(
            [log_noise_noise_probs, log_signal_noise_probs], axis=0
        )
        probs_of_singlet = np.sum(
            [log_noise_noise_probs, log_signal_signal_probs], axis=0
        )
        probs_of_doublet = np.sum(
            [log_noise_signal_probs, log_signal_signal_probs], axis=0
        )
        log_probs_list = [probs_of_negative, probs_of_singlet, probs_of_doublet]

        # each cell and each hypothesis probability
        for prob_idx, log_prob in enumerate(log_probs_list):
            log_likelihoods_for_each_hypothesis[indices, prob_idx] = log_prob
    return (
        log_likelihoods_for_each_hypothesis,
        all_indices,
        counter_to_barcode_combo,
    )


def _calculate_bayes_rule(
    data: np.ndarray, priors: ArrayLike, number_of_noise_barcodes: int
) -> dict[str, np.ndarray]:
    """Calculate bayes rule from log likelihoods.

    Parameters
    ----------
    data
        Anndata object filled only with hashing counts
    priors
        a list of your prior for each hypothesis
        first element is your prior for the negative hypothesis
        second element is your prior for the singlet hypothesis
        third element is your prior for the doublet hypothesis
        We use [0.01, 0.8, 0.19] by default because we assume the barcodes
        in your cell hashing matrix are those cells which have passed QC
        in the transcriptome space, e.g. UMI counts, pct mito reads, etc.
    number_of_noise_barcodes
        number of barcodes to used to calculated noise distribution

    Returns
    -------
    A dict of bayes key results with the following entries:

    `"most_likely_hypothesis"`
        A 1d np.array of the most likely hypothesis
    `"probs_hypotheses"`
        A 2d np.array probability of each hypothesis
    `"log_likelihoods_for_each_hypothesis"`
        A 2d np.array log likelihood of each hypothesis

    """
    priors = np.array(priors)
    log_likelihoods_for_each_hypothesis, _, _ = _calculate_log_likelihoods(
        data, number_of_noise_barcodes
    )
    probs_hypotheses = (
        np.exp(log_likelihoods_for_each_hypothesis)
        * priors
        / np.sum(
            np.multiply(np.exp(log_likelihoods_for_each_hypothesis), priors),
            axis=1,
        )[:, None]
    )
    most_likely_hypothesis = np.argmax(probs_hypotheses, axis=1)
    return {
        "most_likely_hypothesis": most_likely_hypothesis,
        "probs_hypotheses": probs_hypotheses,
        "log_likelihoods_for_each_hypothesis": log_likelihoods_for_each_hypothesis,
    }


@old_positionals(
    "priors", "pre_existing_clusters", "number_of_noise_barcodes", "inplace"
)
@doctest_skip("Illustrative but not runnable doctest code")
def hashsolo(
    adata: AnnData,
    cell_hashing_columns: Sequence[str],
    *,
    priors: tuple[float, float, float] = (0.01, 0.8, 0.19),
    pre_existing_clusters: str | None = None,
    number_of_noise_barcodes: int | None = None,
    inplace: bool = True,
) -> AnnData | None:
    """Probabilistic demultiplexing of cell hashing data using HashSolo :cite:p:`Bernstein2020`.

    .. note::
        More information and bug reports `here <https://github.com/calico/solo>`__.

    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    cell_hashing_columns
        `.obs` columns that contain cell hashing counts.
    priors
        Prior probabilities of each hypothesis, in
        the order `[negative, singlet, doublet]`. The default is set to
        `[0.01, 0.8, 0.19]` assuming barcode counts are from cells that
        have passed QC in the transcriptome space, e.g. UMI counts, pct
        mito reads, etc.
    pre_existing_clusters
        The column in `.obs` containing pre-existing cluster assignments
        (e.g. Leiden clusters or cell types, but not batch assignments).
        If provided, demultiplexing will be performed separately for each
        cluster.
    number_of_noise_barcodes
        The number of barcodes used to create the noise distribution.
        Defaults to `len(cell_hashing_columns) - 2`.
    inplace
        Whether to update `adata` in-place or return a copy.

    Returns
    -------
    A copy of the input `adata` if `inplace=False`, otherwise the input
    `adata`. The following fields are added:

    `.obs["most_likely_hypothesis"]`
        Index of the most likely hypothesis, where `0` corresponds to negative,
        `1` to singlet, and `2` to doublet.
    `.obs["cluster_feature"]`
        The cluster assignments used for demultiplexing.
    `.obs["negative_hypothesis_probability"]`
        Probability of the negative hypothesis.
    `.obs["singlet_hypothesis_probability"]`
        Probability of the singlet hypothesis.
    `.obs["doublet_hypothesis_probability"]`
        Probability of the doublet hypothesis.
    `.obs["Classification"]`:
        Classification of the cell, one of the barcodes in `cell_hashing_columns`,
        `"Negative"`, or `"Doublet"`.

    Examples
    --------
    >>> import anndata
    >>> import scanpy.external as sce
    >>> adata = anndata.read_h5ad("data.h5ad")
    >>> sce.pp.hashsolo(adata, ["Hash1", "Hash2", "Hash3"])
    >>> adata.obs.head()

    """
    print(
        "Please cite HashSolo paper:\nhttps://www.cell.com/cell-systems/fulltext/S2405-4712(20)30195-2"
    )
    adata = adata.copy() if not inplace else adata
    data = adata.obs[cell_hashing_columns].values
    if not check_nonnegative_integers(data):
        msg = "Cell hashing counts must be non-negative"
        raise ValueError(msg)
    if (number_of_noise_barcodes is not None) and (
        number_of_noise_barcodes >= len(cell_hashing_columns)
    ):
        msg = "number_of_noise_barcodes must be at least one less \
        than the number of samples you have as determined by the number of \
        cell_hashing_columns you've given as input  "
        raise ValueError(msg)
    num_of_cells = adata.shape[0]
    results = pd.DataFrame(
        np.zeros((num_of_cells, 6)),
        columns=[
            "most_likely_hypothesis",
            "probs_hypotheses",
            "cluster_feature",
            "negative_hypothesis_probability",
            "singlet_hypothesis_probability",
            "doublet_hypothesis_probability",
        ],
        index=adata.obs_names,
    )
    if pre_existing_clusters is not None:
        cluster_features = pre_existing_clusters
        unique_cluster_features = np.unique(adata.obs[cluster_features])
        for cluster_feature in unique_cluster_features:
            cluster_feature_bool_vector = adata.obs[cluster_features] == cluster_feature
            posterior_dict = _calculate_bayes_rule(
                data[cluster_feature_bool_vector],
                priors,
                number_of_noise_barcodes,
            )
            results.loc[cluster_feature_bool_vector, "most_likely_hypothesis"] = (
                posterior_dict["most_likely_hypothesis"]
            )
            results.loc[cluster_feature_bool_vector, "cluster_feature"] = (
                cluster_feature
            )
            results.loc[
                cluster_feature_bool_vector, "negative_hypothesis_probability"
            ] = posterior_dict["probs_hypotheses"][:, 0]
            results.loc[
                cluster_feature_bool_vector, "singlet_hypothesis_probability"
            ] = posterior_dict["probs_hypotheses"][:, 1]
            results.loc[
                cluster_feature_bool_vector, "doublet_hypothesis_probability"
            ] = posterior_dict["probs_hypotheses"][:, 2]
    else:
        posterior_dict = _calculate_bayes_rule(data, priors, number_of_noise_barcodes)
        results.loc[:, "most_likely_hypothesis"] = posterior_dict[
            "most_likely_hypothesis"
        ]
        results.loc[:, "cluster_feature"] = 0
        results.loc[:, "negative_hypothesis_probability"] = posterior_dict[
            "probs_hypotheses"
        ][:, 0]
        results.loc[:, "singlet_hypothesis_probability"] = posterior_dict[
            "probs_hypotheses"
        ][:, 1]
        results.loc[:, "doublet_hypothesis_probability"] = posterior_dict[
            "probs_hypotheses"
        ][:, 2]

    adata.obs["most_likely_hypothesis"] = results.loc[
        adata.obs_names, "most_likely_hypothesis"
    ]
    adata.obs["cluster_feature"] = results.loc[adata.obs_names, "cluster_feature"]
    adata.obs["negative_hypothesis_probability"] = results.loc[
        adata.obs_names, "negative_hypothesis_probability"
    ]
    adata.obs["singlet_hypothesis_probability"] = results.loc[
        adata.obs_names, "singlet_hypothesis_probability"
    ]
    adata.obs["doublet_hypothesis_probability"] = results.loc[
        adata.obs_names, "doublet_hypothesis_probability"
    ]

    adata.obs["Classification"] = None
    adata.obs.loc[adata.obs["most_likely_hypothesis"] == 2, "Classification"] = (
        "Doublet"
    )
    adata.obs.loc[adata.obs["most_likely_hypothesis"] == 0, "Classification"] = (
        "Negative"
    )
    all_sings = adata.obs["most_likely_hypothesis"] == 1
    singlet_sample_index = np.argmax(
        adata.obs.loc[all_sings, cell_hashing_columns].values, axis=1
    )
    adata.obs.loc[all_sings, "Classification"] = adata.obs[
        cell_hashing_columns
    ].columns[singlet_sample_index]

    return adata if not inplace else None
