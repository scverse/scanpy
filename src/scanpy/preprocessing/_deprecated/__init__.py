from __future__ import annotations

import numpy as np
from scipy import sparse

from ..._compat import CSBase, old_positionals
from ..._utils import dematrix


@old_positionals("max_fraction", "mult_with_mean")
def normalize_per_cell_weinreb16_deprecated(
    x: np.ndarray | CSBase,
    *,
    max_fraction: float = 1,
    mult_with_mean: bool = False,
) -> np.ndarray:
    """Normalize each cell :cite:p:`Weinreb2017`.

    This is a deprecated version. See `normalize_per_cell` instead.

    Normalize each cell by UMI count, so that every cell has the same total
    count.

    Parameters
    ----------
    X
        Expression matrix. Rows correspond to cells and columns to genes.
    max_fraction
        Only use genes that make up more than max_fraction of the total
        reads in every cell.
    mult_with_mean
        Multiply the result with the mean of total counts.

    Returns
    -------
    Normalized version of the original expression matrix.

    """
    if max_fraction < 0 or max_fraction > 1:
        msg = "Choose max_fraction between 0 and 1."
        raise ValueError(msg)

    counts_per_cell = dematrix(x.sum(1)).ravel()
    gene_subset = dematrix(
        np.all(x <= counts_per_cell[:, None] * max_fraction, axis=0)
    ).ravel()
    tc_include = dematrix(x[:, gene_subset].sum(1)).ravel()

    x_norm = (
        x.multiply(sparse.csr_matrix(1 / tc_include[:, None]))  # noqa: TID251
        if isinstance(x, CSBase)
        else x / tc_include[:, None]
    )
    if mult_with_mean:
        x_norm *= np.mean(counts_per_cell)

    return x_norm


def zscore_deprecated(X: np.ndarray) -> np.ndarray:
    """Z-score standardize each variable/gene in X :cite:p:`Weinreb2017`.

    Use `scale` instead.

    Parameters
    ----------
    X
        Data matrix. Rows correspond to cells and columns to genes.

    Returns
    -------
    Z-score standardized version of the data matrix.

    """
    means = np.tile(np.mean(X, axis=0)[None, :], (X.shape[0], 1))
    stds = np.tile(np.std(X, axis=0)[None, :], (X.shape[0], 1))
    return (X - means) / (stds + 0.0001)
