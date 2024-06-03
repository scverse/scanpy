from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

from scanpy.preprocessing._utils import _get_mean_var

from ..._utils import get_random_state

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from ..._utils import AnyRandom


def sparse_multiply(
    E: sparse.csr_matrix | sparse.csc_matrix | NDArray[np.float64],
    a: float | int | NDArray[np.float64],
) -> sparse.csr_matrix | sparse.csc_matrix:
    """multiply each row of E by a scalar"""

    nrow = E.shape[0]
    w = sparse.dia_matrix((a, 0), shape=(nrow, nrow), dtype=a.dtype)
    r = w @ E
    if isinstance(r, np.ndarray):
        return sparse.csc_matrix(r)
    return r


def sparse_zscore(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    gene_mean: NDArray[np.float64] | None = None,
    gene_stdev: NDArray[np.float64] | None = None,
) -> sparse.csr_matrix | sparse.csc_matrix:
    """z-score normalize each column of E"""
    if gene_mean is None or gene_stdev is None:
        gene_means, gene_stdevs = _get_mean_var(E, axis=0)
        gene_stdevs = np.sqrt(gene_stdevs)
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
        current_totals = np.asarray(E.sum(1)).squeeze()
        unsampled_orig_totals = original_totals - current_totals
        unsampled_downsamp_totals = np.random.binomial(
            np.round(unsampled_orig_totals).astype(int), rate
        )
        final_downsamp_totals = current_totals + unsampled_downsamp_totals
    else:
        final_downsamp_totals = original_totals
    return E, final_downsamp_totals
