from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils.numba import njit

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .._compat import CSCBase, CSRBase


@njit
def agg_sum_csr(indicator: CSRBase, data: CSRBase, out: NDArray) -> None:
    for cat_num in numba.prange(indicator.shape[0]):
        start_cat_idx = indicator.indptr[cat_num]
        stop_cat_idx = indicator.indptr[cat_num + 1]
        for row_num in range(start_cat_idx, stop_cat_idx):
            obs_per_cat = indicator.indices[row_num]

            start_obs = data.indptr[obs_per_cat]
            end_obs = data.indptr[obs_per_cat + 1]

            for j in range(start_obs, end_obs):
                col = data.indices[j]
                out[cat_num, col] += data.data[j]


@njit
def agg_sum_csc(indicator: CSRBase, data: CSCBase, out: np.ndarray) -> None:
    obs_to_cat = np.full(data.shape[0], -1, dtype=np.int64)

    for cat in range(indicator.shape[0]):
        for k in range(indicator.indptr[cat], indicator.indptr[cat + 1]):
            obs_to_cat[indicator.indices[k]] = cat

    for col in numba.prange(data.shape[1]):
        start = data.indptr[col]
        end = data.indptr[col + 1]

        for j in range(start, end):
            obs = data.indices[j]
            cat = obs_to_cat[obs]

            if cat != -1:
                out[cat, col] += data.data[j]


@njit
def mean_var_dense(
    indicator: CSRBase, data: NDArray
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    # Welford's online algorithm, parallelized over categories. The indicator
    # CSR lists which observations belong to each category, allowing mask
    # handling to be folded in naturally.
    n_cats = indicator.shape[0]
    n_features = data.shape[1]
    mean = np.zeros((n_cats, n_features), dtype="float64")
    var = np.zeros((n_cats, n_features), dtype="float64")

    for cat in numba.prange(n_cats):
        start = indicator.indptr[cat]
        stop = indicator.indptr[cat + 1]
        n = 0
        for row_num in range(start, stop):
            obs = indicator.indices[row_num]
            n += 1
            for col in range(n_features):
                value = np.float64(data[obs, col])
                delta = value - mean[cat, col]
                mean[cat, col] += delta / n
                delta2 = value - mean[cat, col]
                var[cat, col] += delta * delta2
        if n > 0:
            for col in range(n_features):
                var[cat, col] /= n
    return mean, var


@njit
def mean_var_csr(
    indicator: CSRBase,
    data: CSCBase,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    # Welford's online algorithm over nonzeros, then merge with the block of
    # implicit zeros per (category, feature). Merging a Welford accumulator
    # (n_A, mean_A, M2_A) with k zeros gives:
    #   mean   = mean_A * n_A / (n_A + k)
    #   M2_new = M2_A + mean_A^2 * n_A * k / (n_A + k)
    n_cats = indicator.shape[0]
    n_features = data.shape[1]
    mean = np.zeros((n_cats, n_features), dtype="float64")
    var = np.zeros((n_cats, n_features), dtype="float64")

    for cat_num in numba.prange(n_cats):
        start_cat_idx = indicator.indptr[cat_num]
        stop_cat_idx = indicator.indptr[cat_num + 1]
        n_obs = stop_cat_idx - start_cat_idx
        if n_obs == 0:
            continue

        n_nonzero = np.zeros(n_features, dtype=np.int64)

        for row_num in range(start_cat_idx, stop_cat_idx):
            obs_per_cat = indicator.indices[row_num]
            start_obs = data.indptr[obs_per_cat]
            end_obs = data.indptr[obs_per_cat + 1]

            for j in range(start_obs, end_obs):
                col = data.indices[j]
                value = np.float64(data.data[j])
                n_nonzero[col] += 1
                n = n_nonzero[col]
                delta = value - mean[cat_num, col]
                mean[cat_num, col] += delta / n
                delta2 = value - mean[cat_num, col]
                var[cat_num, col] += delta * delta2

        for col in range(n_features):
            n_nz = n_nonzero[col]
            k = n_obs - n_nz
            if k > 0 and n_nz > 0:
                mean_a = mean[cat_num, col]
                mean[cat_num, col] = mean_a * n_nz / n_obs
                var[cat_num, col] += mean_a * mean_a * n_nz * k / n_obs
            var[cat_num, col] /= n_obs
    return mean, var


@njit
def mean_var_csc(
    indicator: CSRBase, data: CSCBase
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    # Welford's online algorithm, parallelized over columns. For each column
    # we accumulate per-category over the explicit nonzeros, then merge each
    # category's accumulator with its block of implicit zeros (see merge
    # formula in `mean_var_csr`).
    n_cats = indicator.shape[0]
    n_features = data.shape[1]
    obs_to_cat = np.full(data.shape[0], -1, dtype=np.int64)
    n_obs_per_cat = np.zeros(n_cats, dtype=np.int64)
    for cat in range(n_cats):
        n_obs_per_cat[cat] = indicator.indptr[cat + 1] - indicator.indptr[cat]
        for k in range(indicator.indptr[cat], indicator.indptr[cat + 1]):
            obs_to_cat[indicator.indices[k]] = cat

    mean = np.zeros((n_cats, n_features), dtype="float64")
    var = np.zeros((n_cats, n_features), dtype="float64")

    for col in numba.prange(n_features):
        n_nonzero = np.zeros(n_cats, dtype=np.int64)
        start = data.indptr[col]
        end = data.indptr[col + 1]

        for j in range(start, end):
            obs = data.indices[j]
            cat = obs_to_cat[obs]
            if cat == -1:
                continue
            value = np.float64(data.data[j])
            n_nonzero[cat] += 1
            n = n_nonzero[cat]
            delta = value - mean[cat, col]
            mean[cat, col] += delta / n
            delta2 = value - mean[cat, col]
            var[cat, col] += delta * delta2

        for cat in range(n_cats):
            n_obs = n_obs_per_cat[cat]
            if n_obs == 0:
                continue
            n_nz = n_nonzero[cat]
            k = n_obs - n_nz
            if k > 0 and n_nz > 0:
                mean_a = mean[cat, col]
                mean[cat, col] = mean_a * n_nz / n_obs
                var[cat, col] += mean_a * mean_a * n_nz * k / n_obs
            var[cat, col] /= n_obs
    return mean, var
