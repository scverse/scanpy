from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING

import numba
import numpy as np
from scipy import sparse
from sklearn.random_projection import sample_without_replacement

from .._utils import axis_sum, elem_mul

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike, NDArray

    from .._compat import DaskArray
    from .._utils import AnyRandom, _SupportedArray


@singledispatch
def axis_mean(X: DaskArray, *, axis: Literal[0, 1], dtype: DTypeLike) -> DaskArray:
    total = axis_sum(X, axis=axis, dtype=dtype)
    return total / X.shape[axis]


@axis_mean.register(np.ndarray)
def _(X: np.ndarray, *, axis: Literal[0, 1], dtype: DTypeLike) -> np.ndarray:
    return X.mean(axis=axis, dtype=dtype)


def _get_mean_var(
    X: _SupportedArray, *, axis: Literal[0, 1] = 0
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if isinstance(X, sparse.spmatrix):
        mean, var = sparse_mean_variance_axis(X, axis=axis)
    else:
        mean = axis_mean(X, axis=axis, dtype=np.float64)
        mean_sq = axis_mean(elem_mul(X, X), axis=axis, dtype=np.float64)
        var = mean_sq - mean**2
    # enforce R convention (unbiased estimator) for variance
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var


def sparse_mean_variance_axis(mtx: sparse.spmatrix, axis: int):
    """
    This code and internal functions are based on sklearns
    `sparsefuncs.mean_variance_axis`.

    Modifications:
    * allow deciding on the output type, which can increase accuracy when calculating the mean and variance of 32bit floats.
    * This doesn't currently implement support for null values, but could.
    * Uses numba not cython
    """
    assert axis in (0, 1)
    if isinstance(mtx, sparse.csr_matrix):
        ax_minor = 1
        shape = mtx.shape
    elif isinstance(mtx, sparse.csc_matrix):
        ax_minor = 0
        shape = mtx.shape[::-1]
    else:
        raise ValueError("This function only works on sparse csr and csc matrices")
    if axis == ax_minor:
        return sparse_mean_var_major_axis(
            mtx.data,
            mtx.indptr,
            major_len=shape[0],
            minor_len=shape[1],
            n_threads=numba.get_num_threads(),
        )
    else:
        return sparse_mean_var_minor_axis(
            mtx.data,
            mtx.indices,
            mtx.indptr,
            major_len=shape[0],
            minor_len=shape[1],
            n_threads=numba.get_num_threads(),
        )


@numba.njit(cache=True, parallel=True)
def sparse_mean_var_minor_axis(
    data, indices, indptr, *, major_len, minor_len, n_threads
):
    """
    Computes mean and variance for a sparse matrix for the minor axis.

    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    rows = len(indptr) - 1
    sums_minor = np.zeros((n_threads, minor_len))
    squared_sums_minor = np.zeros((n_threads, minor_len))
    means = np.zeros(minor_len)
    variances = np.zeros(minor_len)
    for i in numba.prange(n_threads):
        for r in range(i, rows, n_threads):
            for j in range(indptr[r], indptr[r + 1]):
                minor_index = indices[j]
                if minor_index >= minor_len:
                    continue
                value = data[j]
                sums_minor[i, minor_index] += value
                squared_sums_minor[i, minor_index] += value * value
    for c in numba.prange(minor_len):
        sum_minor = sums_minor[:, c].sum()
        means[c] = sum_minor / major_len
        variances[c] = (
            squared_sums_minor[:, c].sum() / major_len - (sum_minor / major_len) ** 2
        )
    return means, variances


@numba.njit(cache=True, parallel=True)
def sparse_mean_var_major_axis(data, indptr, *, major_len, minor_len, n_threads):
    """
    Computes mean and variance for a sparse array for the major axis.

    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    rows = len(indptr) - 1
    means = np.zeros(major_len)
    variances = np.zeros_like(means)

    for i in numba.prange(n_threads):
        for r in range(i, rows, n_threads):
            sum_major = 0.0
            squared_sum_minor = 0.0
            for j in range(indptr[r], indptr[r + 1]):
                value = np.float64(data[j])
                sum_major += value
                squared_sum_minor += value * value
            means[r] = sum_major
            variances[r] = squared_sum_minor
    for c in numba.prange(major_len):
        mean = means[c] / minor_len
        means[c] = mean
        variances[c] = variances[c] / minor_len - mean * mean
    return means, variances


def sample_comb(
    dims: tuple[int, ...],
    nsamp: int,
    *,
    random_state: AnyRandom = None,
    method: Literal[
        "auto", "tracking_selection", "reservoir_sampling", "pool"
    ] = "auto",
) -> NDArray[np.int64]:
    """Randomly sample indices from a grid, without repeating the same tuple."""
    idx = sample_without_replacement(
        np.prod(dims), nsamp, random_state=random_state, method=method
    )
    return np.vstack(np.unravel_index(idx, dims)).T
