from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from packaging.version import Version
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.utils import check_array, check_random_state
from sklearn.utils.extmath import svd_flip

from ..._compat import pkg_version
from .._utils import _get_mean_var

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray
    from sklearn.decomposition import PCA

    from ..._compat import CSBase
    from ..._utils.random import _LegacyRandom


def _pca_compat_sparse(
    x: CSBase,
    n_pcs: int,
    *,
    solver: Literal["arpack", "lobpcg"],
    mu: NDArray[np.floating] | None = None,
    random_state: _LegacyRandom = None,
) -> tuple[NDArray[np.floating], PCA]:
    """Sparse PCA for scikit-learn <1.4."""
    random_state = check_random_state(random_state)
    np.random.set_state(random_state.get_state())
    random_init = np.random.rand(np.min(x.shape))
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

    X_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (x.shape[0] - 1)

    total_var = _get_mean_var(x)[1].sum()
    ev_ratio = ev / total_var

    from sklearn.decomposition import PCA

    pca = PCA(n_components=n_pcs, svd_solver=solver, random_state=random_state)
    pca.explained_variance_ = ev
    pca.explained_variance_ratio_ = ev_ratio
    pca.components_ = v
    return X_pca, pca
