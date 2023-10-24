from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD

from ..._utils import AnyRandom
from .utils import sparse_var, sparse_multiply, sparse_zscore

if TYPE_CHECKING:
    from .core import Scrublet


def mean_center(self: Scrublet) -> None:
    gene_means = self._E_obs_norm.mean(0)
    self._E_obs_norm = self._E_obs_norm - gene_means
    if self._E_sim_norm is not None:
        self._E_sim_norm = self._E_sim_norm - gene_means


def normalize_variance(self: Scrublet) -> None:
    gene_stdevs = np.sqrt(sparse_var(self._E_obs_norm))
    self._E_obs_norm = sparse_multiply(self._E_obs_norm.T, 1 / gene_stdevs).T
    if self._E_sim_norm is not None:
        self._E_sim_norm = sparse_multiply(self._E_sim_norm.T, 1 / gene_stdevs).T


def zscore(self: Scrublet) -> None:
    gene_means = self._E_obs_norm.mean(0)
    gene_stdevs = np.sqrt(sparse_var(self._E_obs_norm))
    self._E_obs_norm = np.array(
        sparse_zscore(self._E_obs_norm, gene_mean=gene_means, gene_stdev=gene_stdevs)
    )
    if self._E_sim_norm is not None:
        self._E_sim_norm = np.array(
            sparse_zscore(
                self._E_sim_norm, gene_mean=gene_means, gene_stdev=gene_stdevs
            )
        )


def truncated_svd(
    self: Scrublet,
    n_prin_comps: int = 30,
    *,
    random_state: AnyRandom = 0,
    algorithm: Literal["arpack", "randomized"] = "arpack",
) -> None:
    svd = TruncatedSVD(
        n_components=n_prin_comps, random_state=random_state, algorithm=algorithm
    ).fit(self._E_obs_norm)
    self.set_manifold(svd.transform(self._E_obs_norm), svd.transform(self._E_sim_norm))


def pca(
    self: Scrublet,
    n_prin_comps: int = 50,
    *,
    random_state: AnyRandom = 0,
    svd_solver: Literal["auto", "full", "arpack", "randomized"] = "arpack",
) -> None:
    if sparse.issparse(self._E_obs_norm):
        X_obs = self._E_obs_norm.toarray()
    else:
        X_obs = self._E_obs_norm
    if sparse.issparse(self._E_sim_norm):
        X_sim = self._E_sim_norm.toarray()
    else:
        X_sim = self._E_sim_norm

    pca = PCA(
        n_components=n_prin_comps, random_state=random_state, svd_solver=svd_solver
    ).fit(X_obs)
    self.set_manifold(pca.transform(X_obs), pca.transform(X_sim))
