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
            mtx.indptr,
            major_len=shape[0],
            minor_len=shape[1],
            nthr=numba.get_num_threads(),
        )
    else:
        return sparse_mean_var_minor_axis(
            mtx.data,
            mtx.indices,
            mtx.indtpr,
            major_len=shape[0],
            minor_len=shape[1],
            nthr=numba.get_num_threads(),
        )


@numba.njit(cache=True, parallel=True)
def sparse_mean_var_minor_axis(data, indices, indptr, *, major_len, minor_len, nthr):
    """
    Computes mean and variance for a sparse matrix for the minor axis.

    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    nr = len(indptr) - 1
    s = np.zeros((nthr, minor_len))
    ss = np.zeros((nthr, minor_len))
    means = np.zeros(minor_len)
    variances = np.zeros(minor_len)
    for i in numba.prange(nthr):
        for r in range(i, nr, nthr):
            for j in range(indptr[r], indptr[r + 1]):
                c = indices[j]
                if c >= minor_len:
                    continue
                v = data[j]
                s[i, c] += v
                ss[i, c] += v * v
    for c in numba.prange(minor_len):
        s0 = s[:, c].sum()
        means[c] = s0 / major_len
        variances[c] = (
            (ss[:, c].sum() / major_len - (s0 / major_len) ** 2)
            * major_len
            / (major_len - 1)
        )
    return means, variances


@numba.njit(cache=True, parallel=True)
def sparse_mean_var_major_axis(data, indptr, *, major_len, minor_len, nthr):
    """
    Computes mean and variance for a sparse array for the major axis.

    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    nr = len(indptr) - 1
    means = np.zeros(major_len)
    variances = np.zeros_like(means)

    for i in numba.prange(nthr):
        for r in range(i, nr, nthr):
            s = 0.0
            ss = 0.0
            for j in range(indptr[r], indptr[r + 1]):
                v = np.float64(data[j])
                s += v
                ss += v * v
            means[r] = s
            variances[r] = ss
    for c in numba.prange(major_len):
        m0 = means[c] / minor_len
        means[c] = m0
        variances[c] = (
            (variances[c] / minor_len - m0 * m0) * minor_len / (minor_len - 1)
        )
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
