import numpy as np
from scipy.sparse import issparse, csr_matrix


def normalize_per_cell_weinreb16_deprecated(
    X: np.ndarray,
    max_fraction: float = 1,
    mult_with_mean: bool = False,
) -> np.ndarray:
    """\
    Normalize each cell [Weinreb17]_.

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
        raise ValueError('Choose max_fraction between 0 and 1.')

    counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
    gene_subset = np.all(X <= counts_per_cell[:, None] * max_fraction, axis=0)
    if issparse(X):
        gene_subset = gene_subset.A1
    tc_include = (
        X[:, gene_subset].sum(1).A1 if issparse(X) else X[:, gene_subset].sum(1)
    )

    X_norm = (
        X.multiply(csr_matrix(1 / tc_include[:, None]))
        if issparse(X)
        else X / tc_include[:, None]
    )
    if mult_with_mean:
        X_norm *= np.mean(counts_per_cell)

    return X_norm


def zscore_deprecated(X: np.ndarray) -> np.ndarray:
    """\
    Z-score standardize each variable/gene in X.

    Use `scale` instead.

    Reference: Weinreb et al. (2017).

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
