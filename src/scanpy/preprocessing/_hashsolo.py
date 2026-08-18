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
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from .._utils import check_nonnegative_integers
from .._utils._doctests import doctest_skipif
from ..get.get import _get_vec

if TYPE_CHECKING:
    from collections.abc import Collection, Sequence

    from anndata import AnnData
    from anndata.acc import AdRef
    from numpy.typing import ArrayLike, NDArray


class GaussianParams(NamedTuple):
    """Parameters of a gaussian, named after :func:`scipy.stats.norm.pdf`’s parameters."""

    loc: float
    """Mean of the gaussian."""
    scale: float
    """Standard deviation of the gaussian."""


def _calculate_log_likelihoods(
    data: np.ndarray, number_of_noise_barcodes: int | None
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

    def gaussian_updates(data: np.ndarray, mu_o: float, std_o: float) -> GaussianParams:
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
        The updated parameters of the gaussian.

        """
        lam_o = 1 / (std_o**2)
        n = len(data)
        lam = 1 / np.var(data) if len(data) > 1 else lam_o
        lam_n = lam_o + n * lam
        mu_n = (
            (np.mean(data) * n * lam + mu_o * lam_o) / lam_n if len(data) > 0 else mu_o
        )
        return GaussianParams(mu_n, (1 / (lam_n / (n + 1))) ** (1 / 2))

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

    noise_params_dict: dict[int, GaussianParams] = {}
    signal_params_dict: dict[int, GaussianParams] = {}

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
        noise_params_dict[x] = gaussian_updates(
            noise_counts, global_mu_noise_o, global_sigma_noise_o
        )
        signal_params_dict[x] = gaussian_updates(
            signal_counts, global_mu_signal_o, global_sigma_signal_o
        )

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
        data_noise = data_subset[:, noise_sample_idx]
        data_signal = data_subset[:, signal_sample_idx]
        # the distributions the 2nd-highest and highest barcode’s counts
        # are assumed to come from under each hypothesis
        hypotheses = [
            (noise_params, noise_params),  # negative: neither barcode is signal
            (noise_params, signal_params),  # singlet: only the highest barcode is
            (signal_params_dict[noise_sample_idx], signal_params),  # doublet: both are
        ]

        # each cell and each hypothesis probability
        for prob_idx, log_prob in enumerate(
            (
                np.log(norm.pdf(data_noise, loc=d_noise.loc, scale=d_noise.scale) + eps)
                + np.log(norm.pdf(data_signal, loc=d_sig.loc, scale=d_sig.scale) + eps)
                for d_noise, d_sig in hypotheses
            )
        ):
            log_likelihoods_for_each_hypothesis[indices, prob_idx] = log_prob
    return (
        log_likelihoods_for_each_hypothesis,
        all_indices,
        counter_to_barcode_combo,
    )


def _calculate_bayes_rule(
    data: np.ndarray, priors: ArrayLike, number_of_noise_barcodes: int | None
) -> NDArray[np.float64]:
    """Calculate the posterior probability of each hypothesis from log likelihoods.

    Parameters
    ----------
    data
        cells by hashing counts matrix
    priors
        prior for each hypothesis, in the order `[negative, singlet, doublet]`
    number_of_noise_barcodes
        number of barcodes to used to calculated noise distribution

    Returns
    -------
    A 2d array of shape `(n_cells, 3)` with the probability of each hypothesis.

    """
    log_likelihoods, _, _ = _calculate_log_likelihoods(data, number_of_noise_barcodes)
    likelihoods = np.exp(log_likelihoods) * np.asarray(priors)
    return likelihoods / likelihoods.sum(axis=1)[:, None]


def _hashsolo(
    data: NDArray,
    *,
    priors: ArrayLike,
    clusters: NDArray | None,
    number_of_noise_barcodes: int | None,
) -> NDArray[np.float64]:
    """Validate counts and run the bayes rule, optionally per cluster."""
    if not check_nonnegative_integers(data):
        msg = "Cell hashing counts must be non-negative"
        raise ValueError(msg)
    if number_of_noise_barcodes is not None and number_of_noise_barcodes >= (
        n_barcodes := data.shape[1]
    ):
        msg = (
            f"`number_of_noise_barcodes` ({number_of_noise_barcodes}) must be smaller "
            f"than the number of hashing barcodes ({n_barcodes})."
        )
        raise ValueError(msg)

    if clusters is None:
        return _calculate_bayes_rule(data, priors, number_of_noise_barcodes)

    probs = np.zeros((data.shape[0], 3))
    for cluster in pd.unique(clusters):
        mask = clusters == cluster
        probs[mask] = _calculate_bayes_rule(
            data[mask], priors, number_of_noise_barcodes
        )
    return probs


def _ref_name(ref: AdRef | str) -> str:
    """Get the name of the barcode a `ref` points to, e.g. `A.X[:, "Hash1"]` → `Hash1`."""
    if isinstance(ref, str):
        return ref
    idx = ref.idx[-1] if isinstance(ref.idx, tuple) else ref.idx
    return str(idx)


@doctest_skipif(reason="Illustrative but not runnable doctest code")
def hashsolo(
    adata: AnnData,
    hashes: Collection[AdRef | str],
    *,
    priors: tuple[float, float, float] = (0.01, 0.8, 0.19),
    pre_existing_clusters: AdRef | str | None = None,
    number_of_noise_barcodes: int | None = None,
    key_added: str = "hashsolo",
    copy: bool = False,
) -> AnnData | None:
    r"""Probabilistic demultiplexing of cell hashing data using HashSolo :cite:p:`Bernstein2020`.

    .. array-support:: pp.hashsolo

    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    hashes
        References to the vectors holding the cell hashing counts\ [#ref]_,
        e.g. `A.obs[["Hash1", "Hash2"]]` for columns in :attr:`~anndata.AnnData.obs`
        or `A.X[:, ["Hash1", "Hash2"]]` for hashing barcodes among the
        :attr:`~anndata.AnnData.var_names`.
    priors
        Prior probabilities of each hypothesis, in
        the order `[negative, singlet, doublet]`. The default is set to
        `[0.01, 0.8, 0.19]` assuming barcode counts are from cells that
        have passed QC in the transcriptome space, e.g. UMI counts, pct
        mito reads, etc.
    pre_existing_clusters
        Reference to a vector of pre-existing cluster assignments\ [#ref]_
        (e.g. Leiden clusters or cell types, but not batch assignments).
        If provided, demultiplexing is performed separately for each cluster.
    number_of_noise_barcodes
        The number of barcodes used to create the noise distribution.
        Defaults to `len(hashes) - 2`.
    key_added
        Key under which to add the demultiplexing results.
    copy
        Whether to modify a copy of `adata` instead of `adata` itself.

    Returns
    -------
    Returns `None` if `copy=False`, else the modified `adata`.
    Sets the following fields:

    `adata.obs[key_added]` : :class:`~pandas.Categorical` (shape `(n_obs,)`)
        Classification of each cell: the name of one of the `hashes`,
        `"Negative"`, or `"Doublet"`.
    `adata.obsm[key_added]` : :class:`~pandas.DataFrame` (shape `(n_obs, 3)`)
        Probability of the `negative`, `singlet`, and `doublet` hypothesis.

    Examples
    --------
    Hashing counts stored as :attr:`~anndata.AnnData.var_names`
    (e.g. the “Multiplexing Capture” features of a Cell Ranger run):

    >>> import scanpy as sc
    >>> from anndata.acc import A
    >>> adata = sc.io.read_h5ad("data.h5ad")
    >>> sc.pp.hashsolo(adata, A.X[:, ["Hash1", "Hash2", "Hash3"]])
    >>> adata.obs["hashsolo"].value_counts()

    Hashing counts stored as :attr:`~anndata.AnnData.obs` columns:

    >>> sc.pp.hashsolo(adata, A.obs[["Hash1", "Hash2", "Hash3"]])

    .. [#ref] If :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`,
       :class:`str`\ s are :meth:`anndata.acc.AdAcc.resolve`\ d to :class:`~anndata.acc.AdRef`\ s,
       otherwise interpreted as :attr:`anndata.AnnData.obs` columns.

    """
    adata = adata.copy() if copy else adata
    refs = list(hashes)
    data = np.column_stack(_get_vec(adata, refs, dim="obs"))
    clusters = (
        None
        if pre_existing_clusters is None
        else np.asarray(_get_vec(adata, pre_existing_clusters, dim="obs"))
    )
    probs = _hashsolo(
        data,
        priors=priors,
        clusters=clusters,
        number_of_noise_barcodes=number_of_noise_barcodes,
    )

    names = np.array([_ref_name(ref) for ref in refs])
    most_likely_hypothesis = np.argmax(probs, axis=1)
    classification = pd.Series(
        names[np.argmax(data, axis=1)], index=adata.obs_names, dtype="string"
    )
    classification[most_likely_hypothesis == 0] = "Negative"
    classification[most_likely_hypothesis == 2] = "Doublet"

    adata.obs[key_added] = classification.astype("category")
    adata.obsm[key_added] = pd.DataFrame(
        probs, columns=["negative", "singlet", "doublet"], index=adata.obs_names
    )
    return adata if copy else None


def _legacy_hashsolo(
    adata: AnnData,
    cell_hashing_columns: Sequence[str],
    *,
    priors: tuple[float, float, float],
    pre_existing_clusters: str | None,
    number_of_noise_barcodes: int | None,
    inplace: bool,
) -> AnnData | None:
    """Implement the pre-1.13 `scanpy.external.pp.hashsolo` API."""
    adata = adata if inplace else adata.copy()
    cell_hashing_columns = list(cell_hashing_columns)
    data = adata.obs[cell_hashing_columns].to_numpy()
    clusters = (
        None
        if pre_existing_clusters is None
        else adata.obs[pre_existing_clusters].to_numpy()
    )
    probs = _hashsolo(
        data,
        priors=priors,
        clusters=clusters,
        number_of_noise_barcodes=number_of_noise_barcodes,
    )
    most_likely_hypothesis = np.argmax(probs, axis=1)

    adata.obs["most_likely_hypothesis"] = most_likely_hypothesis.astype(float)
    adata.obs["cluster_feature"] = 0.0 if clusters is None else clusters
    for i, hypothesis in enumerate(["negative", "singlet", "doublet"]):
        adata.obs[f"{hypothesis}_hypothesis_probability"] = probs[:, i]

    classification = np.asarray(
        np.array(cell_hashing_columns)[np.argmax(data, axis=1)], dtype=object
    )
    classification[most_likely_hypothesis == 0] = "Negative"
    classification[most_likely_hypothesis == 2] = "Doublet"
    adata.obs["Classification"] = classification

    return adata if not inplace else None
