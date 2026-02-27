from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from fast_array_utils.stats import mean_var
from scipy import sparse

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from ..._compat import CSBase


def sparse_multiply(
    e: CSBase | NDArray[np.float64],
    a: float | NDArray[np.float64],
    /,
) -> CSBase:
    """Multiply each row of E by a scalar."""
    nrow = e.shape[0]
    w = sparse.dia_matrix((a, 0), shape=(nrow, nrow), dtype=a.dtype)
    r = w @ e
    if isinstance(r, np.ndarray):
        return sparse.csc_matrix(r)  # noqa: TID251
    return r


def sparse_zscore(
    e: CSBase,
    /,
    *,
    gene_mean: NDArray[np.float64] | None = None,
    gene_stdev: NDArray[np.float64] | None = None,
) -> CSBase:
    """z-score normalize each column of E."""
    if gene_mean is None or gene_stdev is None:
        _gene_means, gene_stdevs = mean_var(e, axis=0, correction=1)
        gene_stdevs = np.sqrt(gene_stdevs)
    return sparse_multiply(np.asarray((e - gene_mean).T), 1 / gene_stdev).T


def subsample_counts(
    e: CSBase,
    /,
    *,
    rate: float,
    original_totals,
    rng: np.random.Generator,
) -> tuple[CSBase, NDArray[np.int64]]:
    if rate < 1:
        e.data = rng.binomial(np.round(e.data).astype(int), rate)
        current_totals = np.asarray(e.sum(1)).squeeze()
        unsampled_orig_totals = original_totals - current_totals
        unsampled_downsamp_totals = rng.binomial(
            np.round(unsampled_orig_totals).astype(int), rate
        )
        final_downsamp_totals = current_totals + unsampled_downsamp_totals
    else:
        final_downsamp_totals = original_totals
    return e, final_downsamp_totals
