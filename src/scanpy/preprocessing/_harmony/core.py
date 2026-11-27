from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans, kmeans_plusplus
from tqdm.auto import tqdm

from ... import logging as log
from ..._compat import warn
from ..._settings import settings
from ..._settings.verbosity import Verbosity
from ..._utils.random import _legacy_random_state
from ._kernels import _apply_correction, _cluster_pass, _joint_stats

if TYPE_CHECKING:
    from collections.abc import Sequence

    import pandas as pd
    from numpy.typing import NDArray

    from ..._utils.random import RNGLike, SeedLike


__all__ = ["_SUPPRESS_PENALTY", "_compute_lambda_kb", "harmonize"]


_SUPPRESS_PENALTY = 1e30

# flavor -> (max_iter_clustering, tol_clustering, tol_harmony); the "paper" rules use harmony1's for both
_DEFAULTS = {"harmony1": (200, 1e-5, 1e-4), "harmony2": (4, 1e-3, 1e-2)}
# stop_rule -> (window of the clustering test, record every round's objective)
STOP_RULES = {
    # the round objective is only recorded if clustering converged
    "reference": (3, False),
    # every round's objective is recorded, as in rapids-singlecell
    "reference_record": (3, True),
    "paper": (1, True),
    "paper_window3": (3, True),
}
SEEDINGS = ("afkmc2", "kmeans", "kmeans++", "paper_nearest")
# Fixed number of static chunks per parallel sweep, so results do not depend on the thread count.
N_CHUNKS = 32
# clusters per batched solve of the correction; bounds its dense systems
_SOLVE_CHUNK = 8


def harmonize(  # noqa: PLR0913
    x: NDArray[np.floating],
    batch_df: pd.DataFrame,
    batch_key: str | Sequence[str],
    *,
    flavor: str = "harmony2",
    n_clusters: int | None = None,
    sigma: float = 0.1,
    theta: float | Sequence[float] = 2.0,
    tau: int = 0,
    ridge_lambda: float = 1.0,
    alpha: float = 0.2,
    batch_prune_threshold: float | None = 1e-5,
    block_proportion: float = 0.05,
    max_iter_harmony: int = 10,
    max_iter_clustering: int | None = None,
    tol_harmony: float | None = None,
    tol_clustering: float | None = None,
    stop_rule: str = "reference_record",
    seeding: str | None = None,
    penalty_combine: str = "product",
    seeding_lloyd_steps: int = 0,
    seeding_cells_per_cluster: int | None = None,
    recompute_counts: bool = True,
    low_memory: bool = False,
    rng: SeedLike | RNGLike | None = None,
) -> tuple[NDArray[np.floating], dict]:
    """Run Harmony on ``x`` (cells x components) in the dtype of ``x``.

    Returns the corrected embedding and a dict with the Harmony objective per round
    (``objectives``) and the clustering iterations per round (``iterations``).
    """
    params = _HarmonyParams(
        flavor=flavor,
        n_clusters=n_clusters,
        sigma=sigma,
        theta=theta,
        tau=tau,
        ridge_lambda=ridge_lambda,
        alpha=alpha,
        batch_prune_threshold=batch_prune_threshold,
        block_proportion=block_proportion,
        max_iter_harmony=max_iter_harmony,
        max_iter_clustering=max_iter_clustering,
        tol_harmony=tol_harmony,
        tol_clustering=tol_clustering,
        stop_rule=stop_rule,
        seeding=seeding,
        penalty_combine=penalty_combine,
        seeding_lloyd_steps=seeding_lloyd_steps,
        seeding_cells_per_cluster=seeding_cells_per_cluster,
        recompute_counts=recompute_counts,
        low_memory=low_memory,
    )
    x = _check_input(x)
    rng = np.random.default_rng(rng)

    codes, n_levels = _get_batch_codes(batch_df, batch_key)
    # all-zero cells have no direction on the unit sphere: leave them out and return them unchanged
    keep = np.any(x != 0, axis=1)
    harmony = _Harmony(*_drop_zero_cells(x, codes, keep), n_levels, rng, params)
    z_hat, info = harmony.run()
    if not keep.all():
        x_all = x.copy()
        x_all[keep] = z_hat
        return x_all, info
    return z_hat, info


@dataclass(frozen=True, kw_only=True)
class _HarmonyParams:
    """Validated parameters of :func:`harmonize`, see there."""

    flavor: str
    n_clusters: int | None
    sigma: float
    theta: float | Sequence[float]
    tau: int
    ridge_lambda: float
    alpha: float
    batch_prune_threshold: float | None
    block_proportion: float
    max_iter_harmony: int
    max_iter_clustering: int | None
    tol_harmony: float | None
    tol_clustering: float | None
    stop_rule: str
    seeding: str | None
    penalty_combine: str
    seeding_lloyd_steps: int
    seeding_cells_per_cluster: int | None
    recompute_counts: bool
    low_memory: bool

    def __post_init__(self) -> None:
        """Raise a :class:`ValueError` for the first invalid parameter."""
        for name, value, allowed in (
            ("flavor", self.flavor, ("harmony1", "harmony2")),
            ("stop_rule", self.stop_rule, tuple(STOP_RULES)),
            ("seeding", self.seeding_method, SEEDINGS),
            ("penalty_combine", self.penalty_combine, ("product", "sum")),
        ):
            if value not in allowed:
                msg = f"{name} must be one of {allowed}, got {value!r}."
                raise ValueError(msg)
        finite_pos = lambda v: np.isfinite(v) and v > 0  # noqa: E731
        for ok, msg in (
            (self.max_iter_harmony >= 1, "max_iter_harmony must be >= 1"),
            (
                self.max_iter_clustering is None or self.max_iter_clustering >= 1,
                "max_iter_clustering must be >= 1",
            ),
            (
                finite_pos(self.sigma),
                f"sigma must be a finite positive number, got {self.sigma}.",
            ),
            (
                0 < self.block_proportion <= 1,
                f"block_proportion must be in (0, 1], got {self.block_proportion}.",
            ),
            (
                not self.stabilized or finite_pos(self.alpha),
                f"alpha must be a finite positive number when flavor='harmony2', got {self.alpha}.",
            ),
            (
                not self.stabilized
                or self.batch_prune_threshold is None
                or 0 <= self.batch_prune_threshold <= 1,
                f"batch_prune_threshold must be in [0, 1] or None, got {self.batch_prune_threshold}.",
            ),
            (
                self.stabilized or finite_pos(self.ridge_lambda),
                f"ridge_lambda must be a finite positive number when flavor='harmony1', got {self.ridge_lambda}.",
            ),
            (
                self.seeding_lloyd_steps >= 0,
                f"seeding_lloyd_steps must be >= 0, got {self.seeding_lloyd_steps}.",
            ),
            (
                self.seeding_cells_per_cluster is None
                or self.seeding_cells_per_cluster >= 1,
                f"seeding_cells_per_cluster must be >= 1 or None, got {self.seeding_cells_per_cluster}.",
            ),
        ):
            if not ok:
                raise ValueError(msg)

    @property
    def stabilized(self) -> bool:
        return self.flavor == "harmony2"

    @property
    def seeding_method(self) -> str:
        if self.seeding is not None:
            return self.seeding
        return "afkmc2" if self.stabilized else "kmeans"

    @property
    def sum_penalties(self) -> bool:
        return self.penalty_combine == "sum"

    @property
    def window(self) -> int:
        """Window of the clustering convergence test."""
        return STOP_RULES[self.stop_rule][0]

    @property
    def record_always(self) -> bool:
        """Whether to record every round’s objective or only those where clustering converged."""
        return STOP_RULES[self.stop_rule][1]

    @property
    def iteration_limits(self) -> tuple[int, float, float]:
        """``max_iter_clustering``, ``tol_clustering``, and ``tol_harmony``, unset ones filled from ``_DEFAULTS``."""
        d_iter, d_tol_c, d_tol_h = _DEFAULTS[
            "harmony1" if self.stop_rule.startswith("paper") else self.flavor
        ]
        return (
            d_iter if self.max_iter_clustering is None else self.max_iter_clustering,
            d_tol_c if self.tol_clustering is None else self.tol_clustering,
            d_tol_h if self.tol_harmony is None else self.tol_harmony,
        )


class _Harmony:
    """State of one Harmony run on the cells that are not all-zero.

    ``r`` (cells x K) holds the soft cluster assignments, ``o`` (levels x K) and ``r_sum`` (K,) the cluster
    sizes per batch level and overall, ``e`` their expected values, ``y_t`` (d x K) the centroids,
    and ``z_hat`` the corrected embedding.
    """

    y_t: np.ndarray
    e: np.ndarray
    w: np.ndarray

    def __init__(
        self,
        x: NDArray[np.floating],
        codes: NDArray[np.int32],
        n_levels: NDArray[np.int32],
        rng: np.random.Generator,
        params: _HarmonyParams,
    ) -> None:
        ft = x.dtype.type
        self.x, self.codes, self.n_levels, self.rng = x, codes, n_levels, rng
        self.p = params

        n_cells = x.shape[0]
        n_batches = int(n_levels.sum())
        self.plan = _correction_plan(codes, n_batches)
        self.n_b = np.bincount(codes.ravel(), minlength=n_batches).astype(ft)
        self.pr_b = self.n_b / ft(n_cells)
        self.theta_b = _get_theta_array(params.theta, n_levels, x.dtype).ravel()
        self.k = _n_clusters(params.n_clusters, n_cells)
        if params.tau > 0:
            self.theta_b = (
                self.theta_b * (1 - np.exp(-self.n_b / (self.k * params.tau)) ** 2)
            ).astype(ft)

        self.z_norm = _normalize_rows_l2(x)
        self.term = ft(-2.0 / params.sigma)
        self.r = np.zeros((n_cells, self.k), dtype=ft)
        self.o = np.zeros((n_batches, self.k), dtype=ft)
        self.r_sum = np.zeros(self.k, dtype=ft)
        self.bounds = _block_bounds(n_cells, params.block_proportion)
        self.scatter = _use_scatter(n_batches, codes.shape[1], int(self.bounds[1]))
        # low_memory: z_hat shares z_norm's buffer and is rebuilt once at the end
        self.z_hat = self.z_norm if params.low_memory else np.empty_like(x)

    def run(self) -> tuple[NDArray[np.floating], dict]:
        """Alternate clustering and correction until convergence; return the corrected embedding and info."""
        max_iter_clustering, tol_clustering, tol_harmony = self.p.iteration_limits
        objs_harmony = [self.init_clusters()]
        info = dict(objectives=objs_harmony, iterations=[])
        with tqdm(
            range(self.p.max_iter_harmony),
            disable=settings.verbosity < Verbosity.info,
        ) as bar:
            for i in bar:
                objs, done = self.cluster(
                    max_iter=max_iter_clustering, tol=tol_clustering
                )
                info["iterations"].append(len(objs))
                if done or self.p.record_always:
                    objs_harmony.append(objs[-1])
                self.correct()
                if _converged(objs_harmony, tol_harmony, window=1, signed=True):
                    log.info(f"Harmony converged in {i + 1} iterations")
                    break
            else:
                log.info(
                    f"Harmony did not converge after {self.p.max_iter_harmony} iterations."
                )
        return self.result(), info

    def init_clusters(self) -> float:
        """Seed the centroids and assign the cells in one un-penalised pass as a single block; return the objective."""
        ft = self.x.dtype.type
        n_cells, n_batches = self.x.shape[0], self.pr_b.size
        y = _seed_centroids(
            self.z_norm,
            self.k,
            self.rng,
            seeding=self.p.seeding_method,
            lloyd_steps=self.p.seeding_lloyd_steps,
            cells_per_cluster=self.p.seeding_cells_per_cluster,
        )
        dist, ent, _, self.y_t = _cluster_pass(
            self.z_norm,
            np.ascontiguousarray(y.T),
            self.r,
            self.codes,
            self.o,
            self.r_sum,
            self.pr_b,
            np.zeros(n_batches, dtype=ft),
            np.arange(n_cells),
            np.array([0, n_cells]),
            self.term,
            self.p.stabilized,
            N_CHUNKS,
            _use_scatter(n_batches, self.codes.shape[1], n_cells),
            self.p.sum_penalties,
            True,  # noqa: FBT003
        )
        return self._objective(dist, ent)

    def cluster(self, *, max_iter: int, tol: float) -> tuple[list[float], bool]:
        """Iterate the block-wise clustering; return the objective per iteration and whether it converged."""
        window = self.p.window
        objs: list[float] = []
        done = False
        for it in range(max_iter):
            perm = self.rng.permutation(self.x.shape[0])
            # sort each block's cells so the kernel reads z_norm and r in ascending order;
            # block membership is unchanged, only the order within a block differs
            for lo, hi in pairwise(self.bounds):
                perm[lo:hi].sort()
            dist, ent, o_fresh, self.y_t = _cluster_pass(
                self.z_norm,
                self.y_t,
                self.r,
                self.codes,
                self.o,
                self.r_sum,
                self.pr_b,
                self.theta_b,
                perm,
                self.bounds,
                self.term,
                self.p.stabilized,
                N_CHUNKS,
                self.scatter,
                self.p.sum_penalties,
                # the correction recomputes the centroids after the last iteration
                it < max_iter - 1,
            )
            if self.p.recompute_counts:
                self.o[:] = o_fresh
                self.r_sum[:] = o_fresh[: self.n_levels[0]].sum(axis=0)
            objs.append(self._objective(dist, ent))
            if done := _converged(objs, tol, window=window, signed=window > 1):
                break
        return objs, done

    def _objective(self, dist: float, ent: float) -> float:
        """Update the expected counts ``e`` and return the objective with the diversity penalty."""
        self.e = self.pr_b[:, None] * self.r_sum[None, :]
        diversity = _diversity(
            self.o, self.e, self.theta_b, self.p.sigma, stabilized=self.p.stabilized
        )
        return dist + self.p.sigma * ent + diversity

    def correct(self) -> None:
        """Regress the batch effects out of ``x`` into ``z_hat`` and recompute the centroids on it."""
        lam = _compute_lambda_kb(
            self.e,
            o=self.o,
            n_b=self.n_b,
            alpha=self.p.alpha,
            threshold=self.p.batch_prune_threshold,
            ridge_lambda=self.p.ridge_lambda,
            dynamic_lambda=self.p.stabilized,
        )
        self.w = _correction_weights(self.x, self.r, self.plan, self.n_levels, lam)
        self.y_t = _apply_correction(
            self.x,
            self.r,
            self.plan.cat,
            self.w,
            self.z_hat,
            self.z_norm,
            N_CHUNKS,
            True,  # noqa: FBT003
        )

    def result(self) -> NDArray[np.floating]:
        """Return the corrected embedding of the cells, building it if ``low_memory``."""
        if self.p.low_memory:
            _apply_correction(
                self.x,
                self.r,
                self.plan.cat,
                self.w,
                self.z_hat,
                self.z_norm,
                N_CHUNKS,
                False,  # noqa: FBT003
            )
        return self.z_hat


def _n_clusters(n_clusters: int | None, n_cells: int) -> int:
    """Validate ``n_clusters`` or derive it from the number of cells."""
    k = max(2, int(min(100, n_cells / 30))) if n_clusters is None else n_clusters
    if not 1 <= k <= n_cells:
        msg = f"n_clusters must be between 1 and the number of cells ({n_cells}), got {k}."
        raise ValueError(msg)
    return k


def _check_input(x: NDArray[np.floating]) -> NDArray[np.floating]:
    """Validate the embedding and return it C-contiguous."""
    x = np.ascontiguousarray(x)
    if x.dtype not in (np.float32, np.float64):
        msg = f"dtype must be float32 or float64, got {x.dtype}."
        raise ValueError(msg)
    if x.ndim != 2:
        msg = f"x must be a 2-dimensional array (cells x components), got {x.ndim} dimension(s)."
        raise ValueError(msg)
    if not np.isfinite(x).all():
        msg = "Input data contains non-finite values."
        raise ValueError(msg)
    return x


def _drop_zero_cells(
    x: NDArray[np.floating], codes: NDArray[np.int32], keep: NDArray[np.bool]
) -> tuple[NDArray[np.floating], NDArray[np.int32]]:
    """Leave out the all-zero cells (not ``keep``): they have no direction on the unit sphere."""
    if codes.shape[0] != x.shape[0]:  # the kernels index codes without bounds checks
        msg = f"batch_df has {codes.shape[0]} rows but x has {x.shape[0]} cells."
        raise ValueError(msg)
    if not keep.all():
        warn(
            f"{(~keep).sum()} cells have an all-zero embedding; "
            "they are excluded from Harmony and returned unchanged.",
            UserWarning,
        )
        x, codes = np.ascontiguousarray(x[keep]), np.ascontiguousarray(codes[keep])
    if x.shape[0] == 0:
        msg = "No cell has a non-zero embedding."
        raise ValueError(msg)
    return x, codes


def _seed_centroids(
    z_norm: np.ndarray,
    k: int,
    rng: np.random.Generator,
    *,
    seeding: str,
    lloyd_steps: int,
    cells_per_cluster: int | None,
) -> np.ndarray:
    """Seed the centroids (k x d) by ``seeding``, optionally polished by ``lloyd_steps`` Lloyd steps."""
    ft = z_norm.dtype.type
    # optional subsample for the parts of seeding that iterate over cells (KMeans, Lloyd polish);
    # the other seeders read all cells and ignore the cap
    z_seed = z_norm
    uses_cap = seeding == "kmeans" or lloyd_steps > 0
    n_cells = z_norm.shape[0]
    if uses_cap and cells_per_cluster is not None and n_cells > cells_per_cluster * k:
        cells = rng.choice(n_cells, cells_per_cluster * k, replace=False)
        z_seed = z_norm[np.sort(cells)]
    if seeding == "kmeans":
        km = KMeans(
            n_clusters=k,
            random_state=_legacy_random_state(rng, always_state=True),
            max_iter=25,
        )
        y = _normalize_rows_l2(km.fit(z_seed).cluster_centers_.astype(ft))
    elif seeding == "kmeans++":
        y = _seed_kmeanspp(z_norm, k, rng)
    elif seeding == "afkmc2":
        y = _seed_afkmc2(z_norm, k, rng)
    else:
        y = _seed_paper_nearest(z_norm, k, rng)
    if lloyd_steps:
        y = _lloyd_spherical(z_seed, y, lloyd_steps)
    return y


def _get_batch_codes(
    batch_df: pd.DataFrame,
    batch_key: str | Sequence[str],
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Encode each batch variable into a disjoint range of marginal codes."""
    keys = [batch_key] if isinstance(batch_key, str) else list(batch_key)
    if not keys:
        msg = "batch_key must contain at least one column name"
        raise ValueError(msg)

    codes = np.empty((len(batch_df), len(keys)), dtype=np.int32)
    n_levels = np.empty(len(keys), dtype=np.int32)
    offset = 0

    for covariate, key in enumerate(keys):
        batch_vec = batch_df[key].astype("category")
        local_codes = batch_vec.cat.codes.to_numpy(dtype=np.int32, copy=False)
        if np.any(local_codes < 0):
            msg = f"Batch variable {key!r} contains missing values"
            raise ValueError(msg)

        n_categories = batch_vec.cat.categories.size
        n_levels[covariate] = n_categories
        codes[:, covariate] = local_codes + offset
        offset += n_categories

    return codes, n_levels


def _get_theta_array[S: np.generic](
    theta: float | Sequence[float],
    n_levels: NDArray[np.integer],
    dtype: np.dtype[S],
) -> NDArray[S]:
    """Normalize scalar, per-variable, or per-category theta values."""
    levels = n_levels.astype(np.int64, copy=False)
    n_covariates = levels.size
    n_categories = int(levels.sum())

    try:
        theta_array = np.asarray(theta, dtype=dtype)
    except (TypeError, ValueError) as e:
        msg = (
            "theta must be a scalar or an array-like collection of numeric values, "
            f"got {type(theta).__name__}"
        )
        raise ValueError(msg) from e

    if theta_array.ndim == 0:
        return np.full((1, n_categories), theta_array.item(), dtype=dtype)

    theta_array = theta_array.ravel()
    if theta_array.size == n_covariates:
        theta_array = np.repeat(theta_array, levels)
    elif theta_array.size != n_categories:
        msg = (
            f"theta array size ({theta_array.size}) must match the number of batch "
            f"variables ({n_covariates}) or categorical levels ({n_categories})"
        )
        raise ValueError(msg)

    return theta_array.reshape(1, -1)


def _normalize_rows_l2(x: np.ndarray) -> np.ndarray:
    """L2 normalize each row of x."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    return x / norms


def _converged(objs: list[float], tol: float, *, window: int, signed: bool) -> bool:
    """Relative change of the summed last `window` objectives against the `window` before."""
    if len(objs) < window + 1:
        return False
    old, new = sum(objs[-window - 1 : -1]), sum(objs[-window:])
    diff = old - new
    return (diff if signed else abs(diff)) < tol * abs(old)


def _seed_paper_nearest(
    z_norm: np.ndarray, k: int, rng: np.random.Generator
) -> np.ndarray:
    """Seed as :cite:t:`Patikas2026` describe, read so that the seeds spread.

    k distinct cells; each new one is the cell with the highest ``log(u_j) / D_j``, ``u ~ U(0, 1)``
    (the paper prints "lowest", which would pick the closest cells), which draws it with probability
    proportional to its Euclidean distance D_j to the nearest cell chosen so far.
    """
    n = z_norm.shape[0]
    chosen = [int(rng.integers(n))]
    taken = np.zeros(n, dtype=bool)
    dist = np.full(n, np.inf, dtype=z_norm.dtype)
    while len(chosen) < k:
        taken[chosen[-1]] = True
        d = np.sqrt(_sqdist(z_norm, z_norm[chosen[-1]]))
        np.minimum(dist, d, out=dist)
        with np.errstate(divide="ignore"):
            s_j = np.log(rng.random(n)) / dist
        s_j[taken] = -np.inf
        chosen.append(int(np.argmax(s_j)))
    return z_norm[chosen].copy()


def _seed_kmeanspp(z_norm: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """k-means++ seeding, scikit-learn's with one candidate per step (the classic rule)."""
    rs = _legacy_random_state(rng, always_state=True)
    return z_norm[kmeans_plusplus(z_norm, k, random_state=rs, n_local_trials=1)[1]]


def _sqdist(z_norm: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Squared Euclidean distances between unit rows: ``2 - 2 a.b``, clipped at 0."""
    return np.maximum(2.0 - 2.0 * (z_norm @ y.T), 0.0)


def _seed_afkmc2(
    z_norm: np.ndarray, k: int, rng: np.random.Generator, *, chain: int = 200
) -> np.ndarray:
    """AFK-MC² :cite:p:`Bachem2016`: k-means++ approximated by one Metropolis chain per seed.

    Proposal q = d²(x, first seed) / (2 Σ d²) + 1 / (2 n); one pass over the cells, then cost O(k² · chain · d).
    """
    n = z_norm.shape[0]
    first = int(rng.integers(n))
    d2 = _sqdist(z_norm, z_norm[first]).astype(np.float64)
    total = d2.sum()
    q = 0.5 * d2 / total + 0.5 / n if total > 0 else np.full(n, 1.0 / n)
    candidates = rng.choice(n, size=(k - 1, chain), p=q)
    uniforms = rng.random((k - 1, chain))
    seeds = np.empty((k, z_norm.shape[1]), dtype=z_norm.dtype)
    seeds[0] = z_norm[first]
    for c in range(1, k):
        cand = z_norm[candidates[c - 1]]
        # d² to the nearest seed so far
        dist = _sqdist(cand, seeds[:c]).min(axis=1)
        weight = dist / q[candidates[c - 1]]
        cur = 0
        for j in range(1, chain):
            if weight[cur] == 0 or weight[j] > uniforms[c - 1, j] * weight[cur]:
                cur = j
        seeds[c] = cand[cur]
    return seeds


def _lloyd_spherical(z_norm: np.ndarray, y: np.ndarray, steps: int) -> np.ndarray:
    """Hard spherical k-means steps: assign each cell to its most similar centroid, re-normalise the sums."""
    cells = np.arange(z_norm.shape[0])
    for _ in range(steps):
        labels = np.argmax(z_norm @ y.T, axis=1)
        onehot = sparse.csr_matrix(  # noqa: TID251
            (np.ones(cells.size, dtype=z_norm.dtype), (labels, cells)),
            shape=(y.shape[0], cells.size),
        )
        sums = np.asarray(onehot @ z_norm)
        # empty clusters keep their centroid
        y = np.where(
            sums.any(axis=1, keepdims=True), _normalize_rows_l2(sums), y
        ).astype(y.dtype)
    return y


class _Plan(NamedTuple):
    levels: np.ndarray  # J x C, the levels of every joint category
    joint: np.ndarray  # joint category of every cell
    joint_order: np.ndarray
    cols: list[np.ndarray]  # category columns summed by the statistics passes
    orders: list[np.ndarray]
    cat: np.ndarray  # n x C', the categories whose coefficients a cell receives
    by_joint: bool


def _correction_plan(codes: np.ndarray, n_batches: int) -> _Plan:
    """How the correction groups cells: by joint category when there are few, else by each key's level.

    Few joint categories (any single key, nested keys) take one statistics pass and one coefficient row per
    category. Many (crossed keys) take one pass per key, so nothing larger than levels x K x d is held.
    """
    levels, joint = np.unique(codes, axis=0, return_inverse=True)
    joint = np.ascontiguousarray(joint.ravel(), dtype=np.int32)
    joint_order = np.argsort(joint, kind="stable")
    # the cut-off keeps the per-category arrays within twice the per-level ones
    by_joint = levels.shape[0] <= 2 * n_batches
    cols = [joint] if by_joint else [np.ascontiguousarray(c) for c in codes.T]
    orders = [joint_order] if by_joint else [np.argsort(c, kind="stable") for c in cols]
    cat = np.ascontiguousarray(joint[:, None] if by_joint else codes)
    return _Plan(levels, joint, joint_order, cols, orders, cat, by_joint)


def _correction_weights(
    x: np.ndarray,
    r: np.ndarray,
    plan: _Plan,
    n_levels: np.ndarray,
    lam: np.ndarray,
) -> np.ndarray:
    """Coefficients per category of ``plan`` (categories x K x d) for the apply kernels."""
    levels, joint, joint_order, cols, orders, _, by_joint = plan
    n_joint, n_batches = levels.shape[0], int(n_levels.sum())
    if by_joint:
        s_joint, o_joint = _joint_stats(x, r, joint, joint_order, n_joint, N_CHUNKS)
        fold = sparse.csr_matrix(  # noqa: TID251
            (
                np.ones(levels.size),
                (levels.ravel(), np.repeat(np.arange(n_joint), levels.shape[1])),
            ),
            shape=(n_batches, n_joint),
        )
        s = (fold @ s_joint.reshape(n_joint, -1).astype(np.float64)).reshape(
            n_batches, *s_joint.shape[1:]
        )
    else:
        s = sum(
            _joint_stats(x, r, c, order, n_batches, N_CHUNKS)[0]
            for c, order in zip(cols, orders, strict=True)
        )
        # masses per joint category: a pass without the PC dimension
        no_pc = np.empty((x.shape[0], 0), dtype=x.dtype)
        _, o_joint = _joint_stats(no_pc, r, joint, joint_order, n_joint, N_CHUNKS)
    w = _solve_correction(s, o_joint, levels, n_levels, lam)
    if by_joint:
        w = w[levels].sum(axis=1)
    return np.ascontiguousarray(w, dtype=x.dtype)


def _solve_correction(
    s: np.ndarray,
    o_joint: np.ndarray,
    levels: np.ndarray,
    n_levels: np.ndarray,
    lam: np.ndarray,
) -> np.ndarray:
    """Ridge coefficients of every cluster per batch level (L x K x d, float64), without the intercept.

    ``s[l, k] = sum of r_ik * x_i`` over the cells with level l; ``o_joint`` holds the masses per joint category.
    Per cluster the normal equations are ``(Phi* diag(r_k) Phi*^T + Lambda_k) W_k = Phi* diag(r_k) X^T``
    with an unpenalised intercept. The block of the covariate with the most levels is diagonal,
    so it is eliminated in closed form and only the Schur complement over the intercept and the remaining
    levels is solved; with one covariate that is a 1 x 1 system. ``lam[b, k] >= _SUPPRESS_PENALTY`` removes level b
    from cluster k (its coefficient is 0). Runs in float64, ``_SOLVE_CHUNK`` clusters at a time,
    so the dense systems stay at chunk x m x n_big.
    """
    layout = _SchurLayout.from_levels(levels, n_levels)
    o_all = o_joint.astype(np.float64)  # J x K
    s = s.astype(np.float64)  # L x K x d
    lam = lam.astype(np.float64)
    w = np.zeros(s.shape)
    n_clusters = s.shape[1]
    for k0 in range(0, n_clusters, _SOLVE_CHUNK):
        ks = slice(k0, min(k0 + _SOLVE_CHUNK, n_clusters))
        w_big, w_small = _solve_chunk(s[:, ks], o_all[:, ks], lam[:, ks], layout)
        w[layout.big_lv, ks] = w_big.transpose(1, 0, 2)
        # the intercept is not part of the correction
        w[layout.small, ks] = w_small[:, 1:].transpose(1, 0, 2)
    return w


class _SchurLayout(NamedTuple):
    """How :func:`_solve_correction` splits the batch levels.

    The levels of the covariate with the most levels (``big_lv``) are eliminated in closed form;
    the system that remains has a row for the intercept and one per other (``small``) level.
    """

    big_lv: np.ndarray
    big_of_joint: np.ndarray  # each joint category's index among `big_lv`
    small: np.ndarray
    rows: (
        np.ndarray
    )  # J x C, the rows of the small system that each joint category belongs to

    @classmethod
    def from_levels(cls, levels: np.ndarray, n_levels: np.ndarray) -> _SchurLayout:
        offsets = np.concatenate([[0], np.cumsum(n_levels)])
        big = int(np.argmax(n_levels))
        big_lv = np.arange(offsets[big], offsets[big + 1])
        others = [c for c in range(len(n_levels)) if c != big]
        big_of_joint = levels[:, big] - offsets[big]
        small = np.concatenate(
            [np.arange(offsets[c], offsets[c + 1]) for c in others]
            + [np.empty(0, dtype=np.int64)]
        ).astype(np.int64)
        row_of_level = np.full(int(offsets[-1]), -1)
        row_of_level[small] = 1 + np.arange(small.size)
        # the intercept and one level per other key
        rows = np.concatenate(
            [
                np.zeros((len(levels), 1), dtype=np.int64),
                row_of_level[levels[:, others]],
            ],
            axis=1,
        )
        return cls(big_lv, big_of_joint, small, rows)

    @property
    def n_big(self) -> int:
        return self.big_lv.size

    @property
    def m(self) -> int:
        """Size of the small system: the intercept and the levels that are not eliminated."""
        return 1 + self.small.size


def _solve_chunk(
    s: np.ndarray, o: np.ndarray, lam: np.ndarray, layout: _SchurLayout
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the normal equations of a chunk of c clusters (see :func:`_solve_correction`).

    Takes ``s`` (L x c x d), the joint masses ``o`` (J x c) and ``lam`` (L x c), all float64.
    Returns the coefficients of the eliminated levels (c x n_big x d) and of the small system (c x m x d).
    """
    big_lv, big_of_joint, small, rows = layout
    n_big, m = layout.n_big, layout.m
    c, n_pc = o.shape[1], s.shape[2]
    kl = np.arange(c)
    pruned = lam >= _SUPPRESS_PENALTY
    lam_big, pruned_big = lam[big_lv].T, pruned[big_lv].T  # c x n_big
    lam_small, pruned_small = lam[small].T, pruned[small].T  # c x small
    row_k = kl * m + rows[:, :, None]  # J x C x c
    o_big = _masses(kl * n_big + big_of_joint[:, None], o, (c, n_big))
    gram = _masses(row_k[:, :, None] * m + rows[:, None, :, None], o, (c, m, m))
    coupling = _masses(  # mass shared by a small row and a big level
        row_k * n_big + big_of_joint[:, None, None], o, (c, m, n_big)
    )
    rhs_big = s[big_lv].transpose(1, 0, 2)  # c x n_big x d
    rhs = np.empty((c, m, n_pc))
    rhs[:, 0] = rhs_big.sum(axis=1)
    rhs[:, 1:] = s[small].transpose(1, 0, 2)

    diag = o_big + lam_big
    # 0 drops the level
    d_inv = np.where(pruned_big | (diag <= 0), 0.0, 1.0 / np.where(diag > 0, diag, 1.0))
    idx = np.arange(1, m)
    gram[:, idx, idx] += np.where(pruned_small, 0.0, lam_small)
    # a pruned small level: empty row and column, unit diagonal after elimination
    kk, ss = np.nonzero(pruned_small)
    gram[kk, ss + 1, :] = 0
    gram[kk, :, ss + 1] = 0
    coupling[kk, ss + 1, :] = 0
    rhs[kk, ss + 1] = 0
    scaled = coupling * d_inv[:, None, :]  # C D^-1, c x m x n_big
    schur = gram - scaled @ coupling.transpose(0, 2, 1)
    rhs -= scaled @ rhs_big
    schur[kk, ss + 1, ss + 1] = 1
    empty = gram[:, 0, 0] <= 1e-12  # cluster without cells: no correction
    schur[empty] = np.eye(m)
    # positive definite: every unpruned lambda is positive
    w_small = np.linalg.solve(schur, rhs)  # c x m x d
    w_big = d_inv[:, :, None] * (rhs_big - coupling.transpose(0, 2, 1) @ w_small)
    w_small[empty] = 0
    w_big[empty] = 0
    return w_big, w_small


def _masses(idx: np.ndarray, o: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    """Sum the joint masses ``o`` (J x c) into a dense array of ``shape`` at flat indices ``idx`` (J x ... x c)."""
    weights = np.broadcast_to(
        o.reshape(o.shape[0], *[1] * (idx.ndim - 2), -1), idx.shape
    )
    return np.bincount(idx.ravel(), weights.ravel(), int(np.prod(shape))).reshape(shape)


def _compute_lambda_kb(
    e: np.ndarray,
    *,
    o: np.ndarray,
    n_b: np.ndarray,
    alpha: float,
    threshold: float | None,
    ridge_lambda: float,
    dynamic_lambda: bool,
) -> np.ndarray:
    """Compute per-(k,b) ridge regularization array."""
    sentinel = e.dtype.type(_SUPPRESS_PENALTY)
    if not dynamic_lambda:
        lambda_kb = np.full_like(e, ridge_lambda)
    else:
        lambda_kb = (alpha * e).astype(e.dtype)
        if threshold is not None:
            safe_n_b = np.where(n_b > 0, n_b, np.ones_like(n_b))
            prune_mask = (o / safe_n_b[:, None]) < threshold
            prune_mask |= n_b[:, None] == 0
            lambda_kb[prune_mask] = sentinel
    # Where both O and lambda_kb are zero, the kernel computes 1/(O+lambda)
    # which would divide by zero.
    lambda_kb[(o + lambda_kb) == 0] = sentinel
    return lambda_kb


def _use_scatter(n_batches: int, n_cov: int, block_size: int) -> bool:
    """Whether the clustering kernel moves a block's counts by serial scatters instead of per-chunk accumulators.

    Accumulators cost about N_CHUNKS * levels * K per block, scatters about cells * keys * K, but the scatters
    are serial: they win only at about four times the accumulators' work (measured at 200k to 1M cells).
    """
    return N_CHUNKS * n_batches > 4 * block_size * n_cov


def _block_bounds(n_cells: int, block_proportion: float) -> np.ndarray:
    """Block boundaries within a permutation; the split of `np.array_split` into `1 // block_proportion` blocks."""
    n_blocks = int(min(n_cells, 1 // block_proportion))
    sizes = np.full(n_blocks, n_cells // n_blocks)
    sizes[: n_cells % n_blocks] += 1
    return np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)


def _diversity(
    o: np.ndarray, e: np.ndarray, theta: np.ndarray, sigma: float, *, stabilized: bool
) -> float:
    o, e = o.astype(np.float64), e.astype(np.float64)
    ratio = (o + e + 1) / (e + 1) if stabilized else (o + 1) / (e + 1)
    return sigma * float(np.sum(theta[:, None] * o * np.log(ratio)))
