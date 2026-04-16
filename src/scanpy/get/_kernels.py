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
def mean_var_csr(
    indicator: CSRBase,
    data: CSCBase,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    mean = np.zeros((indicator.shape[0], data.shape[1]), dtype="float64")
    var = np.zeros((indicator.shape[0], data.shape[1]), dtype="float64")

    for cat_num in numba.prange(indicator.shape[0]):
        start_cat_idx = indicator.indptr[cat_num]
        stop_cat_idx = indicator.indptr[cat_num + 1]
        for row_num in range(start_cat_idx, stop_cat_idx):
            obs_per_cat = indicator.indices[row_num]

            start_obs = data.indptr[obs_per_cat]
            end_obs = data.indptr[obs_per_cat + 1]

            for j in range(start_obs, end_obs):
                col = data.indices[j]
                value = np.float64(data.data[j])
                value = data.data[j]
                mean[cat_num, col] += value
                var[cat_num, col] += value * value

        n_obs = stop_cat_idx - start_cat_idx
        mean_cat = mean[cat_num, :] / n_obs
        mean[cat_num, :] = mean_cat
        var[cat_num, :] = (var[cat_num, :] / n_obs) - (mean_cat * mean_cat)
    return mean, var


@njit
def mean_var_csc(
    indicator: CSRBase, data: CSCBase
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    obs_to_cat = np.full(data.shape[0], -1, dtype=np.int64)

    mean = np.zeros((indicator.shape[0], data.shape[1]), dtype="float64")
    var = np.zeros((indicator.shape[0], data.shape[1]), dtype="float64")

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
                value = np.float64(data.data[j])
                value = data.data[j]
                mean[cat, col] += value
                var[cat, col] += value * value

    for cat_num in numba.prange(indicator.shape[0]):
        start_cat_idx = indicator.indptr[cat_num]
        stop_cat_idx = indicator.indptr[cat_num + 1]
        n_obs = stop_cat_idx - start_cat_idx
        mean_cat = mean[cat_num, :] / n_obs
        mean[cat_num, :] = mean_cat
        var[cat_num, :] = (var[cat_num, :] / n_obs) - (mean_cat * mean_cat)
    return mean, var
