from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils.numba import njit

if TYPE_CHECKING:
    from .._compat import CSCBase, CSRBase


@njit
def agg_sum_csr(
    indicator: CSRBase,
    data: CSRBase,
):
    out = np.zeros((indicator.shape[0], data.shape[1]))
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
    return out


@njit
def agg_sum_csc(
    indicator: CSRBase,
    data: CSCBase,
):
    out = np.zeros((indicator.shape[0], data.shape[1]), dtype=data.data.dtype)

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

    return out
