from __future__ import annotations

from dataclasses import KW_ONLY, InitVar, dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from sklearn.cluster import KMeans
from tqdm.auto import tqdm

from ... import logging as log
from ..._settings import settings
from ..._settings.verbosity import Verbosity
from ..._utils.random import _legacy_random_state

if TYPE_CHECKING:
    from collections.abc import Sequence

    import pandas as pd

    from ..._utils.random import RNGLike, SeedLike


__all__ = ["_SUPPRESS_PENALTY", "Harmony", "_compute_lambda_kb"]


_SUPPRESS_PENALTY = 1e30


@dataclass
class Harmony:
    """Harmony batch correction algorithm.

    Parameters
    ----------
    x
        Data matrix (n_cells x d) - typically PCA embeddings.
    batch_df
        DataFrame containing batch information.
    batch_key
        Column name(s) in batch_df containing batch labels.
    """

    batch_df: InitVar[pd.DataFrame]
    batch_key: InitVar[str | Sequence[str]]
    _: KW_ONLY
    theta: float | Sequence[float]
    sigma: float
    n_clusters: int | None
    max_iter_harmony: int
    max_iter_clustering: int
    tol_harmony: float
    tol_clustering: float
    ridge_lambda: float
    block_proportion: float
    tau: int
    rng: InitVar[SeedLike | RNGLike | None]
    stabilized_penalty: bool = True
    dynamic_lambda: bool = True
    alpha: float = 0.2
    batch_prune_threshold: float | None = 1e-5

    batch_codes: np.ndarray = field(init=False)
    n_levels: np.ndarray = field(init=False)
    n_batches: int = field(init=False)
    n_covariates: int = field(init=False)
    _rng: np.random.Generator = field(init=False)

    def __post_init__(
        self,
        batch_df: pd.DataFrame,
        batch_key: str | Sequence[str],
        rng: SeedLike | RNGLike | None,
    ) -> None:
        self._rng = np.random.default_rng(rng)

        if self.max_iter_harmony < 1:
            msg = "max_iter_harmony must be >= 1"
            raise ValueError(msg)

        if self.dynamic_lambda:
            if not np.isfinite(self.alpha) or self.alpha <= 0:
                msg = f"alpha must be a finite positive number when dynamic_lambda=True, got {self.alpha}."
                raise ValueError(msg)
            if self.batch_prune_threshold is not None and not (
                0 <= self.batch_prune_threshold <= 1
            ):
                msg = f"batch_prune_threshold must be in [0, 1] or None, got {self.batch_prune_threshold}."
                raise ValueError(msg)

        # Process batch keys
        self.batch_codes, self.n_levels = _get_batch_codes(batch_df, batch_key)
        self.n_batches = int(self.n_levels.sum())
        self.n_covariates = int(self.n_levels.size)

    def fit(self, x: np.ndarray) -> np.ndarray:
        """Run Harmony.

        Returns
        -------
        z_corr
            Batch-corrected embedding matrix (n_cells x d).
        """
        # Ensure input is C-contiguous float array (infer dtype from x)
        x = np.ascontiguousarray(x)
        n_cells = x.shape[0]

        # Normalize input for clustering
        z_norm = _normalize_rows_l2(x)

        # Compute batch proportions
        n_b = np.bincount(self.batch_codes.ravel(), minlength=self.n_batches).astype(
            x.dtype
        )
        pr_b = (n_b / n_cells).reshape(-1, 1)

        theta_arr = _get_theta_array(self.theta, self.n_levels, x.dtype)

        # Set default n_clusters (needed before tau discounting)
        if self.n_clusters is None:
            n_clusters = int(min(100, n_cells / 30))
            n_clusters = max(n_clusters, 2)
        else:
            n_clusters = self.n_clusters

        # Apply tau discounting to theta
        if self.tau > 0:
            theta_arr = theta_arr * (1 - np.exp(-n_b / (n_clusters * self.tau)) ** 2)

        # Initialize centroids and state arrays
        r, e, o, obj_init = _initialize_centroids(
            z_norm,
            self.batch_codes,
            self.n_batches,
            pr_b,
            n_clusters=n_clusters,
            sigma=self.sigma,
            theta=theta_arr,
            rng=self._rng,
            stabilized_penalty=self.stabilized_penalty,
        )

        # Main Harmony loop
        objectives_harmony = [obj_init]
        with tqdm(
            range(self.max_iter_harmony), disable=settings.verbosity < Verbosity.info
        ) as bar:
            for i in bar:
                r, e, o, obj = self._cluster(
                    z_norm, pr_b, r=r, e=e, o=o, theta=theta_arr
                )
                if obj is not None:
                    objectives_harmony.append(obj)

                # Compute per-(k,b) ridge regularization
                lambda_kb = _compute_lambda_kb(
                    e,
                    o=o,
                    n_b=n_b,
                    alpha=self.alpha,
                    threshold=self.batch_prune_threshold,
                    ridge_lambda=self.ridge_lambda,
                    dynamic_lambda=self.dynamic_lambda,
                )

                z_hat = self._correct(x, r, o, lambda_kb=lambda_kb)
                z_norm = _normalize_rows_l2(z_hat)
                if self._is_convergent(objectives_harmony, self.tol_harmony):
                    log.info(f"Harmony converged in {i + 1} iterations")
                    break
            else:
                log.info(
                    f"Harmony did not converge after {self.max_iter_harmony} iterations."
                )

        return z_hat

    @staticmethod
    def _is_convergent(objectives: list[float], tol: float) -> bool:
        """Check Harmony convergence."""
        if len(objectives) < 2:
            return False
        obj_old = objectives[-2]
        obj_new = objectives[-1]
        return (obj_old - obj_new) < tol * abs(obj_old)

    def _cluster(
        self,
        z_norm: np.ndarray,
        pr_b: np.ndarray,
        *,
        r: np.ndarray,
        e: np.ndarray,
        o: np.ndarray,
        theta: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float | None]:
        """Perform clustering step. Modifies r, e, o in-place."""
        return _clustering(
            z_norm,
            self.batch_codes,
            self.n_batches,
            pr_b,
            r=r,
            e=e,
            o=o,
            theta=theta,
            sigma=self.sigma,
            max_iter=self.max_iter_clustering,
            tol=self.tol_clustering,
            block_proportion=self.block_proportion,
            stabilized_penalty=self.stabilized_penalty,
            rng=self._rng,
        )

    def _correct(
        self,
        x: np.ndarray,
        r: np.ndarray,
        o: np.ndarray,
        *,
        lambda_kb: np.ndarray,
    ) -> np.ndarray:
        """Perform correction step."""
        if self.n_covariates > 1:
            return _correction_multi(
                x,
                self.batch_codes,
                self.n_batches,
                r,
                lambda_kb=lambda_kb,
            )

        batch_codes = self.batch_codes[:, 0]
        return _correction_fast(
            x,
            batch_codes,
            self.n_batches,
            r,
            o,
            lambda_kb=lambda_kb,
        )


def _get_batch_codes(
    batch_df: pd.DataFrame,
    batch_key: str | Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
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


def _get_theta_array(
    theta: float | Sequence[float],
    n_levels: np.ndarray,
    dtype: np.dtype,
) -> np.ndarray:
    """Normalize scalar, per-variable, or per-category theta values."""
    levels = np.atleast_1d(n_levels).astype(np.int64, copy=False)
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


def _normalize_rows_l1(r: np.ndarray) -> None:
    """L1 normalize each row of r in-place (rows sum to 1)."""
    row_sums = r.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-12)
    r /= row_sums


def _initialize_centroids(
    z_norm: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    pr_b: np.ndarray,
    *,
    n_clusters: int,
    sigma: float,
    theta: np.ndarray,
    rng: np.random.Generator,
    stabilized_penalty: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Initialize cluster centroids using K-means."""
    kmeans = KMeans(
        n_clusters=n_clusters,
        random_state=_legacy_random_state(rng, always_state=True),
        max_iter=25,
    )
    kmeans.fit(z_norm)

    # Centroids
    y = kmeans.cluster_centers_.copy()
    y_norm = _normalize_rows_l2(y)

    # Compute soft cluster assignments r
    term = -2.0 / sigma
    r = _compute_r(z_norm, y_norm, term)
    _normalize_rows_l1(r)

    # Initialize e (expected) and o (observed)
    r_sum = r.sum(axis=0)
    e = pr_b @ r_sum.reshape(1, -1)
    # o[b, k] = sum of r[i, k] for cells in marginal category b
    o = np.zeros((n_batches, n_clusters), dtype=z_norm.dtype)
    for codes in batch_codes.T:
        np.add.at(o, codes, r)

    # Compute initial objective
    obj = _compute_objective(
        y_norm,
        z_norm,
        r,
        theta=theta,
        sigma=sigma,
        o=o,
        e=e,
        stabilized_penalty=stabilized_penalty,
    )

    return r, e, o, obj


def _compute_r(
    z: np.ndarray,
    y: np.ndarray,
    term: float,
) -> np.ndarray:
    """Compute soft cluster assignments using NumPy dot."""
    dots = z @ y.T
    return np.exp(term * (1.0 - dots))


def _clustering(  # noqa: PLR0913
    z_norm: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    pr_b: np.ndarray,
    *,
    r: np.ndarray,
    e: np.ndarray,
    o: np.ndarray,
    theta: np.ndarray,
    sigma: float,
    max_iter: int,
    tol: float,
    block_proportion: float,
    rng: np.random.Generator,
    stabilized_penalty: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float | None]:
    """Run clustering iterations (modifies r, e, o in-place)."""
    n_cells = z_norm.shape[0]
    k = r.shape[1]
    n_blocks = min(n_cells, 1 // block_proportion)
    term = -2.0 / sigma

    objectives_clustering = []

    # Pre-allocate work arrays
    y = np.empty((k, z_norm.shape[1]), dtype=z_norm.dtype)
    y_norm = np.empty_like(y)

    for _ in range(max_iter):
        # Compute cluster centroids: y = r.T @ z_norm, then normalize
        np.dot(r.T, z_norm, out=y)
        norms = np.linalg.norm(y, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-12)
        np.divide(y, norms, out=y_norm)

        # Randomly shuffle cell indices
        idx_list = rng.permutation(n_cells)

        # Process blocks. Every cell contributes once to each batch variable,
        # while the expected counts depend only on its cluster assignment.
        for cell_idx in np.array_split(idx_list, n_blocks):
            cell_codes = batch_codes[cell_idx]

            # Remove old r contribution from all marginal categories and from e.
            r_old = r[cell_idx]
            r_old_sum = r_old.sum(axis=0)
            for codes in cell_codes.T:
                np.add.at(o, codes, -r_old)
            e -= pr_b * r_old_sum

            # Compute new r values
            dots = z_norm[cell_idx, :] @ y_norm.T
            r_new = np.empty_like(dots)

            # Apply the product of the marginal penalties in log space. This
            # both follows the Harmony formulation and avoids overflow when
            # several variables have large finite penalty factors.
            if stabilized_penalty:
                # Harmony2: denominator is (O + E + 1)
                log_penalty = theta.T * (np.log(e + 1.0) - np.log(o + e + 1.0))
            else:
                # Harmony1: denominator is (O + 1)
                log_penalty = theta.T * (np.log(e + 1.0) - np.log(o + 1.0))
            log_r_new = term * (1.0 - dots)
            log_r_new += log_penalty[cell_codes].sum(axis=1)
            log_r_new -= log_r_new.max(axis=1, keepdims=True)
            np.exp(log_r_new, out=r_new)

            # Normalize rows to sum to 1
            row_sums = r_new.sum(axis=1, keepdims=True)
            row_sums = np.maximum(row_sums, 1e-12)
            r_new /= row_sums

            # Store back
            r[cell_idx, :] = r_new

            # Add new r contribution to o and e
            r_new_sum = r_new.sum(axis=0)
            for codes in cell_codes.T:
                np.add.at(o, codes, r_new)
            e += pr_b * r_new_sum

        # Compute objective
        obj = _compute_objective(
            y_norm,
            z_norm,
            r,
            theta=theta,
            sigma=sigma,
            o=o,
            e=e,
            stabilized_penalty=stabilized_penalty,
        )
        objectives_clustering.append(obj)

        # Check convergence
        if _is_convergent_clustering(objectives_clustering, tol):
            obj = objectives_clustering[-1]
            break
    else:
        obj = None

    return r, e, o, obj


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


def _correction_multi(
    x: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    r: np.ndarray,
    *,
    lambda_kb: np.ndarray,
) -> np.ndarray:
    """Apply the exact ridge correction for a multi-covariate design."""
    _, d = x.shape
    z = x.copy()
    sentinel = x.dtype.type(_SUPPRESS_PENALTY)

    for k_idx, r_k in enumerate(r.T):
        marginal_r = np.zeros(n_batches, dtype=x.dtype)
        rhs = np.zeros((n_batches + 1, d), dtype=x.dtype)
        weighted_x = r_k[:, np.newaxis] * x
        rhs[0] = weighted_x.sum(axis=0)
        for codes in batch_codes.T:
            np.add.at(marginal_r, codes, r_k)
            np.add.at(rhs, codes + 1, weighted_x)

        gram = np.zeros((n_batches + 1, n_batches + 1), dtype=x.dtype)
        gram[0, 0] = r_k.sum()
        gram[0, 1:] = marginal_r
        gram[1:, 0] = marginal_r
        gram[np.arange(1, n_batches + 1), np.arange(1, n_batches + 1)] = marginal_r
        for left in range(batch_codes.shape[1] - 1):
            left_codes = batch_codes[:, left] + 1
            for right in range(left + 1, batch_codes.shape[1]):
                right_codes = batch_codes[:, right] + 1
                np.add.at(gram, (left_codes, right_codes), r_k)
                np.add.at(gram, (right_codes, left_codes), r_k)

        active = lambda_kb[:, k_idx] < sentinel
        retained = np.concatenate(([True], active))
        retained_idx = np.flatnonzero(retained)
        gram[retained_idx[1:], retained_idx[1:]] += lambda_kb[active, k_idx]

        gram_reduced = gram[np.ix_(retained_idx, retained_idx)]
        rhs_reduced = rhs[retained]
        try:
            if np.linalg.matrix_rank(gram_reduced) < gram_reduced.shape[0]:
                raise np.linalg.LinAlgError
            w_reduced = np.linalg.solve(gram_reduced, rhs_reduced)
        except np.linalg.LinAlgError:
            w_reduced = np.linalg.lstsq(gram_reduced, rhs_reduced, rcond=None)[0]

        w = np.zeros((n_batches + 1, d), dtype=x.dtype)
        w[retained] = w_reduced
        w[0] = 0
        correction = np.zeros_like(x)
        for codes in batch_codes.T:
            correction += w[codes + 1]
        z -= r_k[:, np.newaxis] * correction

    return z


def _correction_fast(
    x: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    r: np.ndarray,
    o: np.ndarray,
    *,
    lambda_kb: np.ndarray,
) -> np.ndarray:
    """Fast correction method using precomputed factors."""
    _, d = x.shape

    z = x.copy()
    dtype = x.dtype
    p = np.eye(n_batches + 1, dtype=dtype)

    # Compute all cluster/category right-hand sides together. This avoids
    # rebuilding the same batch masks once per cluster and lets BLAS process
    # the full responsibility matrix for each category.
    phi_t_x = np.empty((r.shape[1], n_batches + 1, d), dtype=dtype)
    phi_t_x[:, 0] = r.T @ x
    for batch in range(n_batches):
        mask = batch_codes == batch
        phi_t_x[:, batch + 1] = r[mask].T @ x[mask]

    for k_idx, (o_k, r_k) in enumerate(zip(o.T, r.T, strict=True)):
        lam_k = lambda_kb[:, k_idx]

        factor = (1.0 / (o_k + lam_k)).astype(dtype)
        c = np.sum(o_k) + np.sum(-factor * o_k**2)
        c_inv = dtype.type(1.0 / c)

        p[0, 1:] = -factor * o_k

        p_t_b_inv = np.zeros((n_batches + 1, n_batches + 1), dtype=dtype)
        p_t_b_inv[0, 0] = c_inv
        p_t_b_inv[1:, 1:] = np.diag(factor)
        p_t_b_inv[1:, 0] = p[0, 1:] * c_inv

        inv_mat = p_t_b_inv @ p

        w = inv_mat @ phi_t_x[k_idx]
        w[0, :] = 0

        w_batch = w[batch_codes + 1, :]
        z -= r_k[:, np.newaxis] * w_batch

    return z


def _compute_objective(
    y_norm: np.ndarray,
    z_norm: np.ndarray,
    r: np.ndarray,
    *,
    theta: np.ndarray,
    sigma: float,
    o: np.ndarray,
    e: np.ndarray,
    stabilized_penalty: bool = True,
) -> float:
    """Compute Harmony objective function."""
    zy = z_norm @ y_norm.T
    kmeans_error = np.sum(r * 2.0 * (1.0 - zy))

    r_row_sums = r.sum(axis=1, keepdims=True)
    r_normalized = r / np.clip(r_row_sums, 1e-12, None)
    entropy = sigma * np.sum(r_normalized * np.log(r_normalized + 1e-12))

    if stabilized_penalty:
        # Harmony2: numerator is (O + E + 1)
        log_ratio = np.log((o + e + 1) / (e + 1))
    else:
        # Harmony1: numerator is (O + 1)
        log_ratio = np.log((o + 1) / (e + 1))
    diversity_penalty = sigma * np.sum(theta @ (o * log_ratio))

    return kmeans_error + entropy + diversity_penalty


def _is_convergent_clustering(
    objectives: list,
    tol: float,
    window_size: int = 3,
) -> bool:
    """Check clustering convergence using window."""
    if len(objectives) < window_size + 1:
        return False

    obj_old = sum(objectives[-window_size - 1 : -1])
    obj_new = sum(objectives[-window_size:])

    return (obj_old - obj_new) < tol * abs(obj_old)
