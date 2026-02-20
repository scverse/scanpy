from __future__ import annotations

from dataclasses import KW_ONLY, InitVar, dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix  # noqa: TID251
from sklearn.cluster import KMeans
from tqdm.auto import tqdm

from ... import logging as log
from ..._settings import settings
from ..._settings.verbosity import Verbosity

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    import pandas as pd

    from ..._compat import CSBase


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
    theta: float | Sequence[float] | None
    sigma: float
    n_clusters: int | None
    max_iter_harmony: int
    max_iter_clustering: int
    tol_harmony: float
    tol_clustering: float
    ridge_lambda: float
    correction_method: Literal["fast", "original"]
    block_proportion: float
    random_state: int | None
    sparse: bool

    batch_codes: np.ndarray = field(init=False)
    n_batches: int = field(init=False)

    def __post_init__(
        self, batch_df: pd.DataFrame, batch_key: str | Sequence[str]
    ) -> None:
        if self.max_iter_harmony < 1:
            msg = "max_iter_harmony must be >= 1"
            raise ValueError(msg)

        # Process batch keys
        self.batch_codes, self.n_batches = _get_batch_codes(batch_df, batch_key)

    def fit(self, x: np.ndarray) -> np.ndarray:
        """Run Harmony.

        Returns
        -------
        z_corr
            Batch-corrected embedding matrix (n_cells x d).
        """
        if self.random_state is not None:
            np.random.seed(self.random_state)

        # Ensure input is C-contiguous float array (infer dtype from x)
        x = np.ascontiguousarray(x)
        n_cells = x.shape[0]

        # Normalize input for clustering
        z_norm = _normalize_rows_l2(x)

        # Build phi matrix (one-hot encoding of batches)
        if self.sparse:
            phi = _one_hot_encode_sparse(self.batch_codes, self.n_batches, x.dtype)
            n_b = np.asarray(phi.sum(axis=0)).ravel()
        else:
            phi = _one_hot_encode(self.batch_codes, self.n_batches, x.dtype)
            n_b = phi.sum(axis=0)
        pr_b = (n_b / n_cells).reshape(-1, 1)

        # Set default theta
        if self.theta is None:
            theta_arr = np.ones(self.n_batches, dtype=x.dtype) * 2.0
        elif isinstance(self.theta, (int, float)):
            theta_arr = np.ones(self.n_batches, dtype=x.dtype) * float(self.theta)
        else:
            theta_arr = np.array(self.theta, dtype=x.dtype)
        theta_arr = theta_arr.reshape(1, -1)

        # Set default n_clusters
        if self.n_clusters is None:
            n_clusters = int(min(100, n_cells / 30))
            n_clusters = max(n_clusters, 2)
        else:
            n_clusters = self.n_clusters

        # Initialize centroids and state arrays
        r, e, o, obj_init = _initialize_centroids(
            z_norm,
            phi,
            pr_b,
            n_clusters=n_clusters,
            sigma=self.sigma,
            theta=theta_arr,
            random_state=self.random_state,
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
                z_hat = self._correct(x, r, o)
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
        """Perform clustering step."""
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
        )

    def _correct(self, x: np.ndarray, r: np.ndarray, o: np.ndarray) -> np.ndarray:
        """Perform correction step."""
        if self.correction_method == "fast":
            return _correction_fast(
                x,
                self.batch_codes,
                self.n_batches,
                r,
                o,
                ridge_lambda=self.ridge_lambda,
            )
        else:
            return _correction_original(
                x,
                self.batch_codes,
                self.n_batches,
                r,
                ridge_lambda=self.ridge_lambda,
            )


def _get_batch_codes(
    batch_df: pd.DataFrame,
    batch_key: str | Sequence[str],
) -> tuple[np.ndarray, int]:
    """Get batch codes from DataFrame."""
    if isinstance(batch_key, str):
        batch_vec = batch_df[batch_key]
    elif len(batch_key) == 1:
        batch_vec = batch_df[batch_key[0]]
    else:
        df = batch_df[list(batch_key)].astype("str")
        batch_vec = df.apply(",".join, axis=1)

    batch_cat = batch_vec.astype("category")
    codes = batch_cat.cat.codes.to_numpy(copy=True)
    n_batches = len(batch_cat.cat.categories)

    return codes.astype(np.int32), n_batches


def _one_hot_encode(
    codes: np.ndarray,
    n_categories: int,
    dtype: np.dtype,
) -> np.ndarray:
    """One-hot encode category codes."""
    n = len(codes)
    phi = np.zeros((n, n_categories), dtype=dtype)
    phi[np.arange(n), codes] = 1.0
    return phi


def _one_hot_encode_sparse(
    codes: np.ndarray,
    n_categories: int,
    dtype: np.dtype,
):
    """One-hot encode category codes as sparse CSR matrix."""
    n = len(codes)
    data = np.ones(n, dtype=dtype)
    indices = codes.astype(np.int32)
    indptr = np.arange(n + 1, dtype=np.int32)
    return csr_matrix((data, indices, indptr), shape=(n, n_categories))


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
    phi: np.ndarray | CSBase,
    pr_b: np.ndarray,
    *,
    n_clusters: int,
    sigma: float,
    theta: np.ndarray,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Initialize cluster centroids using K-means."""
    kmeans = KMeans(
        n_clusters=n_clusters, random_state=random_state, n_init=10, max_iter=25
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
    o = phi.T @ r

    # Compute initial objective
    obj = _compute_objective(y_norm, z_norm, r, theta=theta, sigma=sigma, o=o, e=e)

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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float | None]:
    """Run clustering iterations (modifies r, e, o in-place)."""
    n_cells = z_norm.shape[0]
    k = r.shape[1]
    block_size = max(1, int(n_cells * block_proportion))
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
        idx_list = np.random.permutation(n_cells)

        # Process blocks
        pos = 0
        while pos < n_cells:
            end_pos = min(pos + block_size, n_cells)
            block_idx = idx_list[pos:end_pos]

            for b in range(n_batches):
                mask = batch_codes[block_idx] == b
                if not np.any(mask):
                    continue

                cell_idx = block_idx[mask]

                # Remove old r contribution from o and e
                r_old = r[cell_idx, :]
                r_old_sum = r_old.sum(axis=0)
                o[b, :] -= r_old_sum
                e -= pr_b * r_old_sum

                # Compute new r values
                dots = z_norm[cell_idx, :] @ y_norm.T
                r_new = np.exp(term * (1.0 - dots))

                # Apply penalty
                penalty = ((e[b, :] + 1.0) / (o[b, :] + 1.0)) ** theta[0, b]
                r_new *= penalty

                # Normalize rows to sum to 1
                row_sums = r_new.sum(axis=1, keepdims=True)
                row_sums = np.maximum(row_sums, 1e-12)
                r_new /= row_sums

                # Store back
                r[cell_idx, :] = r_new

                # Add new r contribution to o and e
                r_new_sum = r_new.sum(axis=0)
                o[b, :] += r_new_sum
                e += pr_b * r_new_sum

            pos = end_pos

        # Compute objective
        obj = _compute_objective(y_norm, z_norm, r, theta=theta, sigma=sigma, o=o, e=e)
        objectives_clustering.append(obj)

        # Check convergence
        if _is_convergent_clustering(objectives_clustering, tol):
            obj = objectives_clustering[-1]
            break
    else:
        obj = None

    return r, e, o, obj


def _correction_original(
    x: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    r: np.ndarray,
    *,
    ridge_lambda: float,
) -> np.ndarray:
    """Original correction method - per-cluster ridge regression."""
    _, d = x.shape
    k = r.shape[1]

    # Ridge regularization matrix (don't penalize intercept)
    id_mat = np.eye(n_batches + 1)
    id_mat[0, 0] = 0
    lambda_mat = ridge_lambda * id_mat

    z = x.copy()

    for k_idx in range(k):
        r_k = r[:, k_idx]

        r_sum_total = r_k.sum()
        r_sum_per_batch = np.zeros(n_batches, dtype=x.dtype)
        for b in range(n_batches):
            r_sum_per_batch[b] = r_k[batch_codes == b].sum()

        phi_t_phi = np.zeros((n_batches + 1, n_batches + 1), dtype=x.dtype)
        phi_t_phi[0, 0] = r_sum_total
        phi_t_phi[0, 1:] = r_sum_per_batch
        phi_t_phi[1:, 0] = r_sum_per_batch
        phi_t_phi[1:, 1:] = np.diag(r_sum_per_batch)
        phi_t_phi += lambda_mat

        phi_t_x = np.zeros((n_batches + 1, d), dtype=x.dtype)
        phi_t_x[0, :] = r_k @ x
        for b in range(n_batches):
            mask = batch_codes == b
            phi_t_x[b + 1, :] = r_k[mask] @ x[mask]

        try:
            w = np.linalg.solve(phi_t_phi, phi_t_x)
        except np.linalg.LinAlgError:
            w = np.linalg.lstsq(phi_t_phi, phi_t_x, rcond=None)[0]

        w[0, :] = 0
        w_batch = w[batch_codes + 1, :]
        z -= r_k[:, np.newaxis] * w_batch

    return z


def _correction_fast(
    x: np.ndarray,
    batch_codes: np.ndarray,
    n_batches: int,
    r: np.ndarray,
    o: np.ndarray,
    *,
    ridge_lambda: float,
) -> np.ndarray:
    """Fast correction method using precomputed factors."""
    _, d = x.shape
    k = r.shape[1]

    z = x.copy()
    p = np.eye(n_batches + 1)

    for k_idx in range(k):
        o_k = o[:, k_idx]
        n_k = np.sum(o_k)

        factor = 1.0 / (o_k + ridge_lambda)
        c = n_k + np.sum(-factor * o_k**2)
        c_inv = 1.0 / c

        p[0, 1:] = -factor * o_k

        p_t_b_inv = np.zeros((n_batches + 1, n_batches + 1))
        p_t_b_inv[0, 0] = c_inv
        p_t_b_inv[1:, 1:] = np.diag(factor)
        p_t_b_inv[1:, 0] = p[0, 1:] * c_inv

        inv_mat = p_t_b_inv @ p

        r_k = r[:, k_idx]
        phi_t_x = np.zeros((n_batches + 1, d), dtype=x.dtype)
        phi_t_x[0, :] = r_k @ x
        for b in range(n_batches):
            mask = batch_codes == b
            phi_t_x[b + 1, :] = r_k[mask] @ x[mask]

        w = inv_mat @ phi_t_x
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
) -> float:
    """Compute Harmony objective function."""
    zy = z_norm @ y_norm.T
    kmeans_error = np.sum(r * 2.0 * (1.0 - zy))

    r_row_sums = r.sum(axis=1, keepdims=True)
    r_normalized = r / np.clip(r_row_sums, 1e-12, None)
    entropy = sigma * np.sum(r_normalized * np.log(r_normalized + 1e-12))

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
