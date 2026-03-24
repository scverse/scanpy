from __future__ import annotations

import numpy as np
from scipy import sparse

from ..._compat import CSBase
from ..._utils import dematrix

__all__ = [
    "filter_genes_cv",
    "filter_genes_fano",
    "normalize_per_cell",
    "zscore",
]


def normalize_per_cell(
    x: np.ndarray | CSBase,
    *,
    max_fraction: float = 1,
    mult_with_mean: bool = False,
) -> np.ndarray:
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


def zscore(x: np.ndarray, /) -> np.ndarray:
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
    means = np.tile(np.mean(x, axis=0)[None, :], (x.shape[0], 1))
    stds = np.tile(np.std(x, axis=0)[None, :], (x.shape[0], 1))
    return (x - means) / (stds + 0.0001)


def filter_genes_cv(x, /, e_cutoff, cv_filter):
    """Filter genes by coefficient of variance and mean."""
    return _filter_genes(x, e_cutoff, cv_filter, np.std)


def filter_genes_fano(x, /, e_cutoff, v_cutoff):
    """Filter genes by fano factor and mean."""
    return _filter_genes(x, e_cutoff, v_cutoff, np.var)


def _filter_genes(x, /, e_cutoff, v_cutoff, meth):
    if isinstance(x, CSBase):
        msg = "Not defined for sparse input."
        raise ValueError(msg)
    mean_filter = np.mean(x, axis=0) > e_cutoff
    var_filter = meth(x, axis=0) / (np.mean(x, axis=0) + 0.0001) > v_cutoff
    gene_subset = np.nonzero(np.all([mean_filter, var_filter], axis=0))[0]
    return gene_subset
