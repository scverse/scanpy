from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy import sparse

from ..._utils import AnyRandom, get_random_state

if TYPE_CHECKING:
    from numpy.typing import NDArray


def sparse_var(
    E: sparse.csr_matrix | sparse.csc_matrix,
    *,
    axis: Literal[0, 1],
) -> NDArray[np.float64]:
    """variance across the specified axis"""

    mean_gene: NDArray[np.float64] = E.mean(axis=axis).A.squeeze()
    tmp: sparse.csc_matrix | sparse.csr_matrix = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene**2


def sparse_multiply(
    E: sparse.csr_matrix | sparse.csc_matrix | NDArray[np.float64],
    a: float | int | NDArray[np.float64],
) -> sparse.csr_matrix | sparse.csc_matrix:
    """multiply each row of E by a scalar"""

    nrow = E.shape[0]
    w = sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    r = w @ E
    if isinstance(r, (np.matrix, np.ndarray)):
        return sparse.csc_matrix(r)
    return r


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
        gene_stdev = np.sqrt(sparse_var(E, axis=0))
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
        current_totals = E.sum(1).A.squeeze()
        unsampled_orig_totals = original_totals - current_totals
        unsampled_downsamp_totals = np.random.binomial(
            np.round(unsampled_orig_totals).astype(int), rate
        )
        final_downsamp_totals = current_totals + unsampled_downsamp_totals
    else:
        final_downsamp_totals = original_totals
    return E, final_downsamp_totals
