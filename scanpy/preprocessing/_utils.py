from __future__ import annotations

from functools import singledispatch
from typing import TYPE_CHECKING, Literal

import numba
import numpy as np
from scipy import sparse
from sklearn.random_projection import sample_without_replacement

from .._utils import AnyRandom, _SupportedArray, axis_sum, elem_mul

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .._compat import DaskArray


@singledispatch
def axis_mean(
    X: DaskArray, *, axis: Literal[0, 1], dtype: np.typing.DTypeLike
) -> DaskArray:
    total = axis_sum(X, axis=axis, dtype=dtype)
    return total / X.shape[axis]


@axis_mean.register(np.ndarray)
def _(X: np.ndarray, *, axis: Literal[0, 1], dtype: np.typing.DTypeLike) -> np.ndarray:
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
            mtx.indices,
            mtx.indptr,
            major_len=shape[0],
            minor_len=shape[1],
            dtype=np.float64,
        )
    else:
        return sparse_mean_var_minor_axis(mtx.data, mtx.indices, *shape, np.float64)


@numba.njit(cache=True)
def sparse_mean_var_minor_axis(data, indices, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse matrix for the minor axis.

    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    non_zero = indices.shape[0]

    means = np.zeros(minor_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    counts = np.zeros(minor_len, dtype=np.int64)

    for i in range(non_zero):
        col_ind = indices[i]
        means[col_ind] += data[i]

    for i in range(minor_len):
        means[i] /= major_len

    for i in range(non_zero):
        col_ind = indices[i]
        diff = data[i] - means[col_ind]
        variances[col_ind] += diff * diff
        counts[col_ind] += 1

    for i in range(minor_len):
        variances[i] += (major_len - counts[i]) * means[i] ** 2
        variances[i] /= major_len

    return means, variances


@numba.njit(cache=True)
def sparse_mean_var_major_axis(data, indices, indptr, *, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse array for the major axis.

    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    means = np.zeros(major_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    for i in range(major_len):
        startptr = indptr[i]
        endptr = indptr[i + 1]
        counts = endptr - startptr

        for j in range(startptr, endptr):
            means[i] += data[j]
        means[i] /= minor_len

        for j in range(startptr, endptr):
            diff = data[j] - means[i]
            variances[i] += diff * diff

        variances[i] += (minor_len - counts) * means[i] ** 2
        variances[i] /= minor_len

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
