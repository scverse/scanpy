from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix, issparse

from ..._compat import old_positionals


@old_positionals("max_fraction", "mult_with_mean")
def normalize_per_cell_weinreb16_deprecated(
    x: np.ndarray,
    *,
    max_fraction: float = 1,
    mult_with_mean: bool = False,
) -> np.ndarray:
    """\
    Normalize each cell :cite:p:`Weinreb2017`.

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
        raise ValueError("Choose max_fraction between 0 and 1.")

    counts_per_cell = x.sum(1).A1 if issparse(x) else x.sum(1)
    gene_subset = np.all(x <= counts_per_cell[:, None] * max_fraction, axis=0)
    if issparse(x):
        gene_subset = gene_subset.A1
    tc_include = (
        x[:, gene_subset].sum(1).A1 if issparse(x) else x[:, gene_subset].sum(1)
    )

    x_norm = (
        x.multiply(csr_matrix(1 / tc_include[:, None]))
        if issparse(x)
        else x / tc_include[:, None]
    )
    if mult_with_mean:
        x_norm *= np.mean(counts_per_cell)

    return x_norm


def zscore_deprecated(X: np.ndarray) -> np.ndarray:
    """\
    Z-score standardize each variable/gene in X :cite:p:`Weinreb2017`.

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
