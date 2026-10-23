"""Fused numba kernels for Harmony. Float arrays share one dtype (float32 or float64), C-contiguous.

Layout: ``z_norm``/``x`` cells x d, ``r`` cells x K, ``o`` levels x K,
``y_t`` d x K (centroids transposed, so every inner loop runs over K and vectorises
without float reassociation), ``codes`` cells x C int32 with per-covariate offsets.
``n_chunks`` is the number of static cell chunks per parallel sweep; the driver passes a fixed value,
so results do not depend on the number of numba threads. They differ at rounding level between values of ``n_chunks``.
"""

from __future__ import annotations

import numpy as np
from fast_array_utils.numba import njit
from numba import prange
from numba.extending import register_jitable


@njit
def _cluster_pass(  # noqa: PLR0912, PLR0913, PLR0915, PLR0917
    z_norm,
    y_t,
    r,
    codes,
    o,
    r_sum,
    pr_b,
    theta,
    perm,
    bounds,
    term,
    stabilized,
    n_chunks,
    scatter,
    sum_penalties,
    want_y,
):
    """One clustering iteration: block-wise soft assignment, objective parts, next centroids.

    Mutates ``r``, ``o`` (levels x K) and ``r_sum`` (K,) in place.
    Returns ``(dist, ent, o_fresh, y_t_next)``: ``dist = sum r*2(1-z.y)``, ``ent = sum r*log(r)``
    (float64), ``o`` rebuilt from the new ``r``, and the normalised centroids of the new ``r``.

    ``scatter``: move a block's counts by serial scatters instead of per-chunk accumulators (see ``_use_scatter``).
    ``sum_penalties``: add the keys' penalty factors (over levels with ``theta != 0``) instead of multiplying them.
    ``want_y=False``: skip the centroids and return the input ``y_t``.
    """
    ft = z_norm.dtype.type
    n_pc = z_norm.shape[1]
    n_levels, n_clusters = o.shape
    n_cov = codes.shape[1]
    one = ft(1.0)
    inv_cov = one / ft(n_cov)
    combine_sum = sum_penalties and n_cov > 1
    neg_inf = ft(-np.inf)
    logpen = np.empty((n_levels, n_clusters), dtype=z_norm.dtype)
    acc = np.empty(
        (1 if scatter else n_chunks, n_levels, n_clusters), dtype=z_norm.dtype
    )
    present = np.ones(n_levels, dtype=np.bool)  # levels whose penalty row is needed
    y_acc = np.zeros((n_chunks, n_pc, n_clusters), dtype=z_norm.dtype)
    o_fresh = np.zeros((n_levels, n_clusters), dtype=z_norm.dtype)
    dist = np.zeros(n_chunks, dtype=np.float64)
    ent = np.zeros(n_chunks, dtype=np.float64)
    for blk in range(bounds.size - 1):
        lo = bounds[blk]
        size = bounds[blk + 1] - lo
        # 1. hold the block out of the counts: o -= sum of the block's old r, same for r_sum
        if scatter:
            present[:] = False
            for j in range(lo, lo + size):
                i = perm[j]
                for c in range(n_cov):
                    b = codes[i, c]
                    present[b] = True
                    for k in range(n_clusters):
                        o[b, k] -= r[i, k]
                for k in range(n_clusters):
                    r_sum[k] -= r[i, k]
        else:
            acc[:] = 0
            for t in prange(n_chunks):
                for j in range(
                    lo + t * size // n_chunks, lo + (t + 1) * size // n_chunks
                ):
                    i = perm[j]
                    for c in range(n_cov):
                        b = codes[i, c]
                        for k in range(n_clusters):
                            acc[t, b, k] += r[i, k]
            for t in range(n_chunks):
                for b in range(n_levels):
                    for k in range(n_clusters):
                        o[b, k] -= acc[t, b, k]
                        r_sum[k] -= acc[t, b, k] * inv_cov
        # 2. log penalty table on held-out counts, e[b, k] = pr_b[b] * r_sum[k]
        for b in range(n_levels):
            if present[b]:
                for k in range(n_clusters):
                    e = pr_b[b] * r_sum[k]
                    den = o[b, k] + e + one if stabilized else o[b, k] + one
                    logpen[b, k] = theta[b] * (np.log(e + one) - np.log(den))
        # 3. new assignments; cells of a block are independent given the table
        if not scatter:
            acc[:] = 0
        for t in prange(n_chunks):
            dots = np.empty(n_clusters, dtype=z_norm.dtype)
            logit = np.empty(n_clusters, dtype=z_norm.dtype)
            y_loc = y_acc[t]
            dist_t = 0.0
            ent_t = 0.0
            for j in range(lo + t * size // n_chunks, lo + (t + 1) * size // n_chunks):
                i = perm[j]
                for k in range(n_clusters):
                    dots[k] = 0
                for dd in range(n_pc):
                    zv = z_norm[i, dd]
                    for k in range(n_clusters):
                        dots[k] += zv * y_t[dd, k]
                for k in range(n_clusters):
                    logit[k] = term * (one - dots[k])
                if combine_sum:
                    # log of the summed factors of the keys with theta > 0, shifted by their maximum
                    for k in range(n_clusters):
                        top_pen = neg_inf
                        for c in range(n_cov):
                            b = codes[i, c]
                            if theta[b] != 0:
                                top_pen = max(top_pen, logpen[b, k])
                        if top_pen > neg_inf:
                            summed = ft(0.0)
                            for c in range(n_cov):
                                b = codes[i, c]
                                if theta[b] != 0:
                                    summed += np.exp(logpen[b, k] - top_pen)
                            logit[k] += top_pen + np.log(summed)
                else:
                    for c in range(n_cov):
                        b = codes[i, c]
                        for k in range(n_clusters):
                            logit[k] += logpen[b, k]
                top = logit[0]
                for k in range(1, n_clusters):
                    top = max(top, logit[k])
                total = ft(0.0)
                for k in range(n_clusters):
                    logit[k] -= top
                    p = np.exp(logit[k])
                    r[i, k] = p
                    total += p
                inv = one / total
                cell_dist = ft(0.0)
                cell_ent = ft(0.0)
                for k in range(n_clusters):
                    rk = r[i, k] * inv
                    r[i, k] = rk
                    cell_dist += rk * (one - dots[k])
                    cell_ent += rk * logit[k]  # log r_k = logit_k - log(total)
                dist_t += 2.0 * cell_dist
                ent_t += cell_ent - np.log(total)
                if not scatter:
                    for c in range(n_cov):
                        b = codes[i, c]
                        for k in range(n_clusters):
                            acc[t, b, k] += r[i, k]
                if want_y:
                    for dd in range(n_pc):
                        zv = z_norm[i, dd]
                        for k in range(n_clusters):
                            y_loc[dd, k] += zv * r[i, k]
            dist[t] += dist_t
            ent[t] += ent_t
        # 4. put the block back
        if scatter:
            for j in range(lo, lo + size):
                i = perm[j]
                for c in range(n_cov):
                    b = codes[i, c]
                    for k in range(n_clusters):
                        o[b, k] += r[i, k]
                        o_fresh[b, k] += r[i, k]
                for k in range(n_clusters):
                    r_sum[k] += r[i, k]
        else:
            for t in range(n_chunks):
                for b in range(n_levels):
                    for k in range(n_clusters):
                        o[b, k] += acc[t, b, k]
                        o_fresh[b, k] += acc[t, b, k]
                        r_sum[k] += acc[t, b, k] * inv_cov
    # serial sums: a reduction over the parallel chunks would depend on the thread count
    dist_total = ent_total = 0.0
    for t in range(n_chunks):
        dist_total += dist[t]
        ent_total += ent[t]
    if not want_y:
        return dist_total, ent_total, o_fresh, y_t
    # 5. centroids of the new r: y_k = normalise(sum_i r_ik z_i)
    return dist_total, ent_total, o_fresh, _centroids(y_acc)


@register_jitable
def _centroids(y_acc):
    """Sum per-chunk centroid accumulators (chunks x d x K) serially and L2-normalise every column."""
    ft = y_acc.dtype.type
    n_chunks, n_pc, n_clusters = y_acc.shape
    y_t = np.empty((n_pc, n_clusters), dtype=y_acc.dtype)
    for k in range(n_clusters):
        sq = ft(0.0)
        for dd in range(n_pc):
            v = ft(0.0)
            for t in range(n_chunks):
                v += y_acc[t, dd, k]
            y_t[dd, k] = v
            sq += v * v
        scale = ft(1.0) / max(np.sqrt(sq), ft(1e-12))
        for dd in range(n_pc):
            y_t[dd, k] *= scale
    return y_t


@njit
def _joint_stats(x, r, joint, order, n_joint, n_chunks):  # noqa: PLR0912, PLR0917
    """Sufficient statistics of the per-cluster regression, one row per joint batch category.

    ``joint[i]`` is the index of cell i's combination of batch levels (with one key: its batch) and ``order``
    lists the cells sorted by it. Returns ``s[j, k, :] = sum of r_ik * x_i`` and ``o[j, k] = sum of r_ik`` over
    the cells of category j.

    The sorted cells are cut into ``n_chunks`` equal pieces. A category that lies strictly inside a piece is
    written by that piece alone; only the first and last category of a piece can continue in a neighbour, and
    those two go through a small per-piece buffer that is merged serially. Scratch is 2 * n_chunks * K * d floats
    whatever the number of categories.
    """
    n_cells, n_pc = x.shape
    n_clusters = r.shape[1]
    s = np.zeros((n_joint, n_clusters, n_pc), dtype=x.dtype)
    o = np.zeros((n_joint, n_clusters), dtype=x.dtype)
    s_edge = np.zeros((n_chunks, 2, n_clusters, n_pc), dtype=x.dtype)
    o_edge = np.zeros((n_chunks, 2, n_clusters), dtype=x.dtype)
    first = np.full(n_chunks, -1, dtype=np.int64)
    last = np.full(n_chunks, -1, dtype=np.int64)
    for t in range(n_chunks):
        lo = t * n_cells // n_chunks
        hi = (t + 1) * n_cells // n_chunks
        if hi > lo:
            first[t] = joint[order[lo]]
            last[t] = joint[order[hi - 1]]
    for t in prange(n_chunks):
        for p in range(t * n_cells // n_chunks, (t + 1) * n_cells // n_chunks):
            i = order[p]
            j = joint[i]
            if j == first[t] or j == last[t]:
                side = 0 if j == first[t] else 1
                for k in range(n_clusters):
                    rk = r[i, k]
                    o_edge[t, side, k] += rk
                    for dd in range(n_pc):
                        s_edge[t, side, k, dd] += rk * x[i, dd]
            else:
                for k in range(n_clusters):
                    rk = r[i, k]
                    o[j, k] += rk
                    for dd in range(n_pc):
                        s[j, k, dd] += rk * x[i, dd]
    for t in range(n_chunks):
        if first[t] < 0:
            continue
        for side in range(2):
            j = first[t] if side == 0 else last[t]
            for k in range(n_clusters):
                o[j, k] += o_edge[t, side, k]
                for dd in range(n_pc):
                    s[j, k, dd] += s_edge[t, side, k, dd]
    return s, o


@njit
def _apply_correction(x, r, cat, w, z_hat, z_norm, n_chunks, normalise):  # noqa: PLR0917
    """Write ``z_hat = x_i - sum_k r_ik * sum_c w[cat_ic, k, :]`` in place.

    With ``normalise``, also write its row-normalised copy into ``z_norm`` and return the normalised centroids
    ``y_t`` (d x K) of ``r`` on it; otherwise return an empty array. Each row is written in full before it is read,
    so ``z_hat`` and ``z_norm`` may be the same array (the low_memory path).
    """
    ft = x.dtype.type
    n_cells, n_pc = x.shape
    n_clusters = r.shape[1]
    y_acc = np.zeros((n_chunks, n_pc, n_clusters), dtype=x.dtype)
    for t in prange(n_chunks):
        y_loc = y_acc[t]
        for i in range(t * n_cells // n_chunks, (t + 1) * n_cells // n_chunks):
            for dd in range(n_pc):
                z_hat[i, dd] = x[i, dd]
            for k in range(n_clusters):
                rk = r[i, k]
                for c in range(cat.shape[1]):
                    j = cat[i, c]
                    for dd in range(n_pc):
                        z_hat[i, dd] -= rk * w[j, k, dd]
            if not normalise:
                continue
            sq = ft(0.0)
            for dd in range(n_pc):
                sq += z_hat[i, dd] * z_hat[i, dd]
            scale = ft(1.0) / max(np.sqrt(sq), ft(1e-12))
            for dd in range(n_pc):
                zv = z_hat[i, dd] * scale
                z_norm[i, dd] = zv
                for k in range(n_clusters):
                    y_loc[dd, k] += zv * r[i, k]
    if not normalise:
        return np.empty((0, 0), dtype=x.dtype)
    return _centroids(y_acc)
