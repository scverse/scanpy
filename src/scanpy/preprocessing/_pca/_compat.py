from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from fast_array_utils.stats import mean_var
from packaging.version import Version
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.utils import check_array
from sklearn.utils.extmath import svd_flip

from ..._compat import pkg_version
from ..._utils.random import _accepts_legacy_random_state, _legacy_random_state

if TYPE_CHECKING:
    from typing import Any, Literal

    from numpy.typing import NDArray
    from sklearn.decomposition import PCA

    from ..._compat import CSBase
    from ..._utils.random import RNGLike, SeedLike


SCIPY_1_15 = pkg_version("scipy") >= Version("1.15.0rc1")


@_accepts_legacy_random_state(None)
def _pca_compat_sparse(
    x: CSBase,
    n_pcs: int,
    *,
    solver: Literal["arpack", "lobpcg"],
    mu: NDArray[np.floating] | None = None,
    rng: SeedLike | RNGLike | None = None,
) -> tuple[NDArray[np.floating], PCA]:
    """Sparse PCA for scikit-learn <1.4."""
    rng = np.random.default_rng(rng)
    # this exists only to be stored in our PCA container object
    random_state_meta = _legacy_random_state(rng)

    x = check_array(x, accept_sparse=["csr", "csc"])

    if mu is None:
        mu = np.asarray(x.mean(0)).flatten()[None, :]
    ones = np.ones(x.shape[0])[None, :].dot

    def mat_op(v: NDArray[np.floating]):
        return (x @ v) - (mu @ v)

    def rmat_op(v: NDArray[np.floating]):
        return (x.T.conj() @ v) - (mu.T @ ones(v))

    linop = LinearOperator(
        dtype=x.dtype,
        shape=x.shape,
        matvec=mat_op,
        matmat=mat_op,
        rmatvec=rmat_op,
        rmatmat=rmat_op,
    )

    random_init = rng.uniform(size=np.min(x.shape))
    kw: Any = (
        dict(rng=rng) if SCIPY_1_15 else dict(random_state=_legacy_random_state(rng))
    )
    u, s, v = svds(linop, solver=solver, k=n_pcs, v0=random_init, **kw)
    # u_based_decision was changed in https://github.com/scikit-learn/scikit-learn/pull/27491
    u, v = svd_flip(u, v, u_based_decision=not SCIPY_1_15)
    idx = np.argsort(-s)
    v = v[idx, :]

    x_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (x.shape[0] - 1)

    total_var = mean_var(x, correction=1, axis=0)[1].sum()
    ev_ratio = ev / total_var

    from sklearn.decomposition import PCA

    pca = PCA(n_components=n_pcs, svd_solver=solver, random_state=random_state_meta)
    pca.explained_variance_ = ev
    pca.explained_variance_ratio_ = ev_ratio
    pca.components_ = v
    return x_pca, pca


@_accepts_legacy_random_state(None)
def _pca_compat_sparse_row_center(
    x: CSBase,
    row_center: NDArray[np.floating],
    n_pcs: int,
    *,
    solver: Literal["arpack", "lobpcg"],
    rng: SeedLike | RNGLike | None = None,
) -> tuple[NDArray[np.floating], PCA]:
    """Sparse PCA for matrices represented as ``x - row_center[:, None]``."""
    rng = np.random.default_rng(rng)
    random_state_meta = _legacy_random_state(rng)

    x = check_array(x, accept_sparse=["csr", "csc"])
    row_center = np.asarray(row_center, dtype=np.float64).ravel()
    if row_center.shape[0] != x.shape[0]:
        msg = (
            f"`row_center` length must match the number of observations, "
            f"got {row_center.shape[0]} for {x.shape[0]} observations."
        )
        raise ValueError(msg)

    n_obs, _ = x.shape
    row_center_sum = row_center.sum()
    mu = (np.asarray(x.sum(axis=0)).ravel() - row_center_sum) / n_obs
    mu = mu[None, :]
    ones = np.ones(n_obs)[None, :].dot

    def mat_op(v: NDArray[np.floating]):
        row_scale = np.sum(v, axis=0)
        row_offset = (
            row_center * row_scale
            if np.ndim(row_scale) == 0
            else row_center[:, None] * row_scale
        )
        return (x @ v) - row_offset - (mu @ v)

    def rmat_op(v: NDArray[np.floating]):
        row_offset = row_center @ v
        if np.ndim(row_offset) > 0:
            row_offset = row_offset[None, :]
        return (x.T.conj() @ v) - row_offset - (mu.T @ ones(v))

    linop = LinearOperator(
        dtype=x.dtype,
        shape=x.shape,
        matvec=mat_op,
        matmat=mat_op,
        rmatvec=rmat_op,
        rmatmat=rmat_op,
    )

    random_init = rng.uniform(size=np.min(x.shape))
    kw: Any = (
        dict(rng=rng) if SCIPY_1_15 else dict(random_state=_legacy_random_state(rng))
    )
    u, s, v = svds(linop, solver=solver, k=n_pcs, v0=random_init, **kw)
    u, v = svd_flip(u, v, u_based_decision=not SCIPY_1_15)
    idx = np.argsort(-s)
    v = v[idx, :]

    x_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (n_obs - 1)

    x_sum = np.asarray(x.sum(axis=0)).ravel()
    x_sq_sum = np.asarray(x.multiply(x).sum(axis=0)).ravel()
    x_row_center_sum = np.asarray(x.T @ row_center).ravel()
    y_sum = x_sum - row_center_sum
    y_sq_sum = x_sq_sum - 2.0 * x_row_center_sum + np.sum(row_center**2)
    total_var = np.sum((y_sq_sum - y_sum**2 / n_obs) / (n_obs - 1))
    ev_ratio = ev / total_var

    from sklearn.decomposition import PCA

    pca = PCA(n_components=n_pcs, svd_solver=solver, random_state=random_state_meta)
    pca.explained_variance_ = ev
    pca.explained_variance_ratio_ = ev_ratio
    pca.components_ = v
    return x_pca, pca
