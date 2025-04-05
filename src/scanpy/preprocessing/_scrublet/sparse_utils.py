from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

from scanpy.preprocessing._utils import _get_mean_var

from ..._utils import _get_legacy_random

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from ..._compat import CSBase, _LegacyRandom


def sparse_multiply(
    E: CSBase | NDArray[np.float64],
    a: float | NDArray[np.float64],
) -> CSBase:
    """Multiply each row of E by a scalar."""
    nrow = E.shape[0]
    w = sparse.dia_matrix((a, 0), shape=(nrow, nrow), dtype=a.dtype)
    r = w @ E
    if isinstance(r, np.ndarray):
        return sparse.csc_matrix(r)  # noqa: TID251
    return r


def sparse_zscore(
    E: CSBase,
    *,
    gene_mean: NDArray[np.float64] | None = None,
    gene_stdev: NDArray[np.float64] | None = None,
) -> CSBase:
    """z-score normalize each column of E."""
    if gene_mean is None or gene_stdev is None:
        gene_means, gene_stdevs = _get_mean_var(E, axis=0)
        gene_stdevs = np.sqrt(gene_stdevs)
    return sparse_multiply(np.asarray((E - gene_mean).T), 1 / gene_stdev).T


def subsample_counts(
    E: CSBase,
    *,
    rate: float,
    original_totals,
    random_seed: _LegacyRandom = 0,
) -> tuple[CSBase, NDArray[np.int64]]:
    if rate < 1:
        random_seed = _get_legacy_random(random_seed)
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
