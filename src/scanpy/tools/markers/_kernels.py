from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils.numba import njit
from scipy import sparse

from ..._compat import CSBase
from ...get._aggregated import _chan_combine

if TYPE_CHECKING:
    from collections.abc import Generator

    from numpy.typing import NDArray


_CONST_MAX_SIZE: int = 10_000_000


@njit
def rankdata(data: NDArray[np.number]) -> NDArray[np.float64]:
    """Parallelized version of scipy.stats.rankdata."""
    ranked = np.empty(data.shape, dtype=np.float64)
    for j in numba.prange(data.shape[1]):
        arr = np.ravel(data[:, j])
        sorter = np.argsort(arr)

        arr = arr[sorter]
        obs = np.concatenate((np.array([True]), arr[1:] != arr[:-1]))

        dense = np.empty(obs.size, dtype=np.int64)
        dense[sorter] = obs.cumsum()

        # cumulative counts of each unique value
        count = np.concatenate((np.flatnonzero(obs), np.array([len(obs)])))
        ranked[:, j] = 0.5 * (count[dense] + count[dense - 1] + 1)

    return ranked


@njit
def _tiecorrect(rankvals: NDArray[np.number]) -> NDArray[np.float64]:
    """Parallelized version of scipy.stats.tiecorrect."""
    tc = np.ones(rankvals.shape[1], dtype=np.float64)
    for j in numba.prange(rankvals.shape[1]):
        arr = np.sort(np.ravel(rankvals[:, j]))
        idx = np.flatnonzero(
            np.concatenate((np.array([True]), arr[1:] != arr[:-1], np.array([True])))
        )
        cnt = np.diff(idx).astype(np.float64)

        size = np.float64(arr.size)
        if size >= 2:
            tc[j] = 1.0 - (cnt**3 - cnt).sum() / (size**3 - size)

    return tc


def _dense_chunk(
    x: NDArray[np.number] | CSBase,
    left: int,
    right: int,
    mask_obs: NDArray[np.bool] | None,
    mask_obs_rest: NDArray[np.bool] | None,
) -> NDArray[np.number]:
    if mask_obs is None or mask_obs_rest is None:
        chunk = x[:, left:right]
        return chunk.toarray() if isinstance(chunk, CSBase) else chunk
    blocks = (x[mask_obs, left:right], x[mask_obs_rest, left:right])
    if isinstance(x, CSBase):
        return sparse.vstack(blocks).toarray()
    return np.vstack(blocks)


def _ranks(
    x: NDArray[np.number] | CSBase,
    /,
    mask_obs: NDArray[np.bool] | None = None,
    mask_obs_rest: NDArray[np.bool] | None = None,
) -> Generator[tuple[NDArray[np.float64], int, int]]:
    n_genes = x.shape[1]

    if mask_obs is not None and mask_obs_rest is not None:
        n_cells = np.count_nonzero(mask_obs) + np.count_nonzero(mask_obs_rest)
    else:
        n_cells = x.shape[0]

    # Calculate chunk frames
    max_chunk = max(_CONST_MAX_SIZE // n_cells, 1)

    for left in range(0, n_genes, max_chunk):
        right = min(left + max_chunk, n_genes)

        ranks = rankdata(_dense_chunk(x, left, right, mask_obs, mask_obs_rest))
        yield ranks, left, right


@numba.njit  # noqa: TID251  (inner kernel called from _vars_rest's nopython loop)
def _chan_accumulate(
    group_counts: NDArray[np.float64],
    mean: NDArray[np.float64],
    m2: NDArray[np.float64],
    j: int,
    direction: int,
) -> NDArray[np.float64]:
    """Accumulate a running Chan combine of the groups for gene ``j``.

    ``acc[i]`` holds the combined ``(count, mean, M2)`` of group ``i`` together with
    every group toward ``direction``: ``+1`` gives groups ``0..i`` (forward),
    ``-1`` gives groups ``i..end`` (backward).
    """
    n_groups = group_counts.shape[0]
    acc = np.empty((n_groups, 3))
    # accumulated (count, mean, M2) of the groups visited so far
    acc_count = acc_mean = acc_m2 = 0.0
    # visit groups forward (direction +1) or backward (-1)
    group_order = range(n_groups) if direction == 1 else range(n_groups - 1, -1, -1)
    for i in group_order:
        acc_count, acc_mean, acc_m2 = _chan_combine(
            acc_count, acc_mean, acc_m2, group_counts[i], mean[i, j], m2[i, j]
        )
        acc[i, 0], acc[i, 1], acc[i, 2] = acc_count, acc_mean, acc_m2
    return acc


@njit
def _vars_rest(
    group_counts: NDArray[np.float64],
    mean: NDArray[np.float64],
    m2: NDArray[np.float64],
    k: int,
) -> NDArray[np.float64]:
    """Leave-one-out variance for each selected group, parallel over genes.

    Group ``g``'s "rest" is every other group combined — the groups up to ``g - 1``
    pooled with the groups from ``g + 1`` — so variances are never subtracted
    (Chan's cancellation-free combine).
    """
    n_genes = mean.shape[1]
    vars_rest = np.zeros((k, n_genes))
    for j in numba.prange(n_genes):
        # combined (count, mean, M2) of groups 0..i (forward) and i..end (backward)
        combined_stats_upto_group = _chan_accumulate(group_counts, mean, m2, j, 1)
        combined_stats_from_group = _chan_accumulate(group_counts, mean, m2, j, -1)

        # each group g's "rest" stats = the groups before g pooled with the groups after g
        for g in range(k):
            stats_after_g = combined_stats_from_group[g + 1]
            if g >= 1:
                stats_before_g = combined_stats_upto_group[g - 1]
                # each stats row is (count, mean, M2)
                n_r, _, m2_r = _chan_combine(
                    n_a=stats_before_g[0],
                    mean_a=stats_before_g[1],
                    m2_a=stats_before_g[2],
                    n_b=stats_after_g[0],
                    mean_b=stats_after_g[1],
                    m2_b=stats_after_g[2],
                )
            else:
                # g == 0 has no groups before it, so its rest is just the groups after
                n_r, m2_r = stats_after_g[0], stats_after_g[2]
            denom = n_r - 1.0
            v = m2_r / denom
            vars_rest[g, j] = v
    return vars_rest
