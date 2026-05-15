from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from fast_array_utils.numba import njit


if TYPE_CHECKING:
    from numpy.typing import DTypeLike, NDArray
    from scanpy._compat import CSRBase

@njit
def _csr_gram_upper_triangular(
    mat: CSRBase,
    out: NDArray,
) -> None:
    """Accumulate upper triangle of ``X.T @ X`` into ``out``, then write over the lower."""
    n_rows = mat.indptr.shape[0] - 1
    for i in range(n_rows):
        row_start = mat.indptr[i]
        row_end = mat.indptr[i + 1]
        for k in range(row_start, row_end):
            jk = mat.indices[k]
            vk = mat.data[k]
            for l in range(k, row_end):
                out[jk, mat.indices[l]] += vk * mat.data[l]
    for i in range(out.shape[0]):
        for j in range(i):
            out[i, j] = out[j, i]

@njit
def _csr_gram_full(
    mat: CSRBase,
    out: NDArray,
) -> None:
    """Accumulate ``X.T @ X`` into ``out``."""
    n_rows = mat.indptr.shape[0] - 1
    for i in range(n_rows):
        row_start = mat.indptr[i]
        row_end = mat.indptr[i + 1]
        for k in range(row_start, row_end):
            jk = mat.indices[k]
            vk = mat.data[k]
            for l in range(row_start, row_end):
                out[jk, mat.indices[l]] += vk * mat.data[l]


def csr_gram_dense(mat: CSRBase) -> NDArray:
    """Compute the full dense gram matrix ``X.T @ X`` for a CSR ``X``."""
    out = np.zeros((mat.shape[1], mat.shape[1]), dtype=mat.dtype)
    if mat.has_sorted_indices:
        _csr_gram_upper_triangular(mat, out)
    else:
        _csr_gram_full(mat, out)
    return out
