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
from fast_array_utils.conv import to_dense
from scipy.stats import norm

from .._compat import CSBase
from .._utils import check_nonnegative_integers
from .._utils._doctests import doctest_needs
from ..get.get import MultiAcc, _get_arr, _get_vec

if TYPE_CHECKING:
    from collections.abc import Collection

    from anndata import AnnData
    from anndata.acc import AdRef
    from numpy.typing import ArrayLike, NDArray


# smaller than what log-counts can resolve
# large enough that its inverse doesn’t get close to overflowing
_VAR_EPS = 1e-8
"""Variance floor, so degenerate spreads (e.g. a one-cell cluster) stay finite."""


class Gaussian(NamedTuple):
    """A gaussian, its fields named after :func:`scipy.stats.norm.pdf`’s parameters."""

    loc: float | np.floating
    """Mean of the gaussian."""
    scale: float | np.floating
    """Standard deviation of the gaussian."""

    def log_pdf(self, counts: NDArray[np.floating]) -> NDArray[np.float64]:
        """Log of the probability density at each of `counts`."""
        eps = 1e-15  # avoid log(0)
        return np.log(norm.pdf(counts, loc=self.loc, scale=self.scale) + eps)

    def posterior(self, data: NDArray[np.floating]) -> Gaussian:
        """Update this gaussian, used as a prior, with the observed 1-d `data`.

        See <https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf>.
        """
        n = len(data)
        lam_o = 1 / max(self.scale**2, _VAR_EPS)
        lam = 1 / max(np.var(data), _VAR_EPS) if n > 1 else lam_o
        lam_n = lam_o + n * lam
        mu_n = (np.mean(data) * n * lam + self.loc * lam_o) / lam_n if n else self.loc
        return Gaussian(mu_n, np.sqrt((n + 1) / lam_n))


def _calculate_log_likelihoods(
    data: NDArray[np.integer], n_barcodes_noise: int
) -> NDArray[np.float64]:
    """Calculate log likelihoods for each hypothesis, negative, singlet, doublet.

    Parameters
    ----------
    data
        cells by hashing counts matrix
    n_barcodes_noise
        number of barcodes to used to calculated noise distribution

    Returns
    -------
    A 2d array of shape `(n_cells, 3)` with the log likelihood of each hypothesis.

    """
    # probabilites for negative, singlet, doublets
    log_likelihoods = np.zeros((data.shape[0], 3))

    n_barcodes = data.shape[1]

    # assume log normal
    data: NDArray[np.floating] = np.log(data + 1)
    # per cell, the barcode indices ordered by ascending count
    data_arg = np.argsort(data, axis=1)
    data_sort = np.sort(data, axis=1)

    # global signal and noise counts useful for when we have few cells
    # barcodes with the highest number of counts are assumed to be a true signal
    # barcodes with rank < k are considered to be noise
    global_signal = data_sort[:, -1]
    global_noise = data_sort[:, :n_barcodes_noise]
    prior_sig = Gaussian(np.mean(global_signal), np.std(global_signal))
    prior_noise = Gaussian(np.mean(global_noise), np.std(global_noise))

    # for each barcode get empirical noise and signal distribution parameterization,
    # assuming lognormal, as an update from the global values
    p_noise: list[Gaussian] = []
    p_sig: list[Gaussian] = []
    for x in range(n_barcodes):
        is_noise = (data_arg[:, :n_barcodes_noise] == x).any(axis=1)
        is_sig = data_arg[:, -1] == x
        p_noise.append(prior_noise.posterior(data[is_noise, x]))
        p_sig.append(prior_sig.posterior(data[is_sig, x]))

    # for each combination of noise and signal barcode calculate probiltiy of in silico and real cell hypotheses
    for i_noise, i_sig in product(range(n_barcodes), repeat=2):
        mask = (data_arg[:, -1] == i_sig) & (data_arg[:, -2] == i_noise)
        if not mask.any():
            continue

        # the distributions the 2nd-highest and highest barcode’s counts
        # are assumed to come from under each hypothesis
        hypotheses = [
            (p_noise[i_noise], p_noise[i_noise]),  # negative: neither barcode is signal
            (p_noise[i_noise], p_sig[i_sig]),  # singlet: only the highest barcode is
            (p_sig[i_noise], p_sig[i_sig]),  # doublet: both are
        ]
        for i_prob, (p_2nd, p_top) in enumerate(hypotheses):
            log_likelihoods[mask, i_prob] = p_2nd.log_pdf(
                data[mask, i_noise]
            ) + p_top.log_pdf(data[mask, i_sig])
    return log_likelihoods


def _calculate_bayes_rule(
    data: NDArray[np.integer], priors: ArrayLike, n_barcodes_noise: int
) -> NDArray[np.float64]:
    """Calculate the posterior probability of each hypothesis from log likelihoods.

    Parameters
    ----------
    data
        cells by hashing counts matrix
    priors
        prior for each hypothesis, in the order `[negative, singlet, doublet]`
    n_barcodes_noise
        number of barcodes to used to calculated noise distribution

    Returns
    -------
    A 2d array of shape `(n_cells, 3)` with the probability of each hypothesis.

    """
    log_likelihoods = _calculate_log_likelihoods(data, n_barcodes_noise)
    likelihoods = np.exp(log_likelihoods) * np.asarray(priors)
    return likelihoods / likelihoods.sum(axis=1)[:, None]


def _hashsolo(
    data: NDArray[np.integer],
    *,
    priors: ArrayLike,
    clusters: NDArray | None,
    n_barcodes_noise: int | None,
) -> NDArray[np.float64]:
    """Validate counts and run the bayes rule, optionally per cluster."""
    if not check_nonnegative_integers(data):
        msg = "Cell hashing counts must be non-negative integers"
        raise ValueError(msg)
    n_barcodes = data.shape[1]
    if n_barcodes_noise is None:
        n_barcodes_noise = n_barcodes - 2
    if not 1 <= n_barcodes_noise < n_barcodes:
        msg = (
            f"The number of noise barcodes ({n_barcodes_noise}) must be at least 1 and smaller "
            f"than the number of hashing barcodes ({n_barcodes}). "
            f"Pass at least 3 `hashes` or set `n_noise_barcodes` explicitly."
        )
        raise ValueError(msg)

    if clusters is None:
        return _calculate_bayes_rule(data, priors, n_barcodes_noise)

    probs = np.zeros((data.shape[0], 3))
    # `factorize` then `np.unique` gives cells with no cluster the code -1
    codes = pd.factorize(clusters)[0]
    for code in np.unique(codes):
        mask = codes == code
        probs[mask] = _calculate_bayes_rule(data[mask], priors, n_barcodes_noise)
    return probs


def _ref_name(ref: AdRef | str) -> str:
    """Get the name of the barcode a `ref` points to, e.g. `A.X[:, "Hash1"]` → `Hash1`."""
    if isinstance(ref, str):
        return ref
    idx = ref.idx[-1] if isinstance(ref.idx, tuple) else ref.idx
    return str(idx)


def _get_hashes(
    adata: AnnData, hashes: MultiAcc | Collection[AdRef] | Collection[str]
) -> tuple[NDArray, NDArray[np.str_]]:
    """Get the hashing count matrix and each hashing barcode’s name."""
    if not isinstance(hashes, MultiAcc):
        data = np.column_stack(_get_vec(adata, hashes, dim="obs"))
        return data, np.array([_ref_name(ref) for ref in hashes])
    # a `MultiAcc` such as `A.obsm["hto"]` refers to all of its columns
    data = _get_arr(adata, hashes, dim="obs")
    if isinstance(data, pd.DataFrame):
        return data.to_numpy(), data.columns.to_numpy(dtype=str)
    if isinstance(data, CSBase):
        data = to_dense(data)
    return data, np.arange(data.shape[1]).astype(str)


@doctest_needs("anndata_acc")
def hashsolo(
    adata: AnnData,
    hashes: MultiAcc | Collection[AdRef] | Collection[str],
    *,
    priors: tuple[float, float, float] = (0.01, 0.8, 0.19),
    pre_existing_clusters: AdRef | str | None = None,
    n_noise_barcodes: int | None = None,
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
        References to the vectors holding the cell hashing counts,
        e.g. `A.obs[["Hash1", "Hash2"]]` for columns in :attr:`~anndata.AnnData.obs`
        or `A.X[:, ["Hash1", "Hash2"]]` for hashing barcodes among the
        :attr:`~anndata.AnnData.var_names`.
        A :class:`~anndata.acc.MultiAcc` such as `A.obsm["hto"]` refers to
        *all* columns of what it points to,
        named after the :class:`~pandas.DataFrame`’s columns or,
        for a plain array, its column positions.
    priors
        Prior probabilities of each hypothesis, in
        the order `[negative, singlet, doublet]`. The default is set to
        `[0.01, 0.8, 0.19]` assuming barcode counts are from cells that
        have passed QC in the transcriptome space, e.g. UMI counts, pct
        mito reads, etc.
    pre_existing_clusters
        Reference to a vector of pre-existing cluster assignments
        (e.g. Leiden clusters or cell types, but not batch assignments).
        If provided, demultiplexing is performed separately for each cluster.
    n_noise_barcodes
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
    Simulate 300 cells, each carrying one of 3 hashtag oligos:

    >>> import numpy as np
    >>> import scanpy as sc
    >>> from anndata import AnnData
    >>> from anndata.acc import A
    >>>
    >>> rng = np.random.default_rng(0)
    >>> hto = rng.poisson(20, size=(300, 3))  # ambient background
    >>> hto[np.arange(300), np.arange(300) % 3] = rng.poisson(1000, 300)  # signal
    >>> adata = AnnData(rng.poisson(1, (300, 5)).astype("f4"), obsm=dict(hto=hto))

    A :class:`~anndata.acc.MultiAcc` demultiplexes using every column it points to.
    A plain array has no column names, so the barcodes are named after their positions:

    >>> sc.pp.hashsolo(adata, A.obsm["hto"])
    >>> adata.obs["hashsolo"].cat.categories.astype("string")
    Index(['0', '1', '2'], dtype='string')
    >>> adata.obs["hashsolo"].value_counts().tolist()
    [100, 100, 100]

    The same counts as :attr:`~anndata.AnnData.obs` columns
    (as in a Cell Ranger run’s “Multiplexing Capture” features),
    referenced one by one:

    >>> adata.obs[["Hash1", "Hash2", "Hash3"]] = hto
    >>> sc.pp.hashsolo(adata, A.obs[["Hash1", "Hash2", "Hash3"]], key_added="hs2")
    >>> adata.obs["hs2"].cat.categories.astype("string")
    Index(['Hash1', 'Hash2', 'Hash3'], dtype='string')

    """
    adata = adata.copy() if copy else adata
    data, names = _get_hashes(adata, hashes)
    clusters = (
        None
        if pre_existing_clusters is None
        else np.asarray(_get_vec(adata, pre_existing_clusters, dim="obs"))
    )
    probs = _hashsolo(
        data, priors=priors, clusters=clusters, n_barcodes_noise=n_noise_barcodes
    )

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
