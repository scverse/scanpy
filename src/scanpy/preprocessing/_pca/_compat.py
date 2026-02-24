from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from fast_array_utils.stats import mean_var
from packaging.version import Version
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.utils import check_array
from sklearn.utils.extmath import svd_flip

from ..._compat import pkg_version
from ..._utils.random import accepts_legacy_random_state, legacy_random_state

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray
    from sklearn.decomposition import PCA

    from ..._compat import CSBase
    from ..._utils.random import RNGLike, SeedLike


@accepts_legacy_random_state(None)
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
    random_init = rng.uniform(size=np.min(x.shape))
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

    u, s, v = svds(linop, solver=solver, k=n_pcs, v0=random_init)
    # u_based_decision was changed in https://github.com/scikit-learn/scikit-learn/pull/27491
    u, v = svd_flip(
        u, v, u_based_decision=pkg_version("scikit-learn") < Version("1.5.0rc1")
    )
    idx = np.argsort(-s)
    v = v[idx, :]

    x_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (x.shape[0] - 1)

    total_var = mean_var(x, correction=1, axis=0)[1].sum()
    ev_ratio = ev / total_var

    from sklearn.decomposition import PCA

    pca = PCA(
        n_components=n_pcs, svd_solver=solver, random_state=legacy_random_state(rng)
    )
    pca.explained_variance_ = ev
    pca.explained_variance_ratio_ = ev_ratio
    pca.components_ = v
    return x_pca, pca
