from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy import sparse

from .sparse_utils import sparse_multiply, sparse_var, sparse_zscore

if TYPE_CHECKING:
    from ..._utils import AnyRandom
    from .core import Scrublet


def mean_center(self: Scrublet) -> None:
    gene_means = self._counts_obs_norm.mean(0)
    self._counts_obs_norm = sparse.csc_matrix(self._counts_obs_norm - gene_means)
    if self._counts_sim_norm is not None:
        self._counts_sim_norm = sparse.csc_matrix(self._counts_sim_norm - gene_means)


def normalize_variance(self: Scrublet) -> None:
    gene_stdevs = np.sqrt(sparse_var(self._counts_obs_norm, axis=0))
    self._counts_obs_norm = sparse_multiply(self._counts_obs_norm.T, 1 / gene_stdevs).T
    if self._counts_sim_norm is not None:
        self._counts_sim_norm = sparse_multiply(
            self._counts_sim_norm.T, 1 / gene_stdevs
        ).T


def zscore(self: Scrublet) -> None:
    gene_means = self._counts_obs_norm.mean(0)
    gene_stdevs = np.sqrt(sparse_var(self._counts_obs_norm, axis=0))
    self._counts_obs_norm = sparse_zscore(
        self._counts_obs_norm, gene_mean=gene_means, gene_stdev=gene_stdevs
    )
    if self._counts_sim_norm is not None:
        self._counts_sim_norm = sparse_zscore(
            self._counts_sim_norm, gene_mean=gene_means, gene_stdev=gene_stdevs
        )


def truncated_svd(
    self: Scrublet,
    n_prin_comps: int = 30,
    *,
    random_state: AnyRandom = 0,
    algorithm: Literal["arpack", "randomized"] = "arpack",
) -> None:
    if self._counts_sim_norm is None:
        raise RuntimeError("_counts_sim_norm is not set")
    from sklearn.decomposition import TruncatedSVD

    svd = TruncatedSVD(
        n_components=n_prin_comps, random_state=random_state, algorithm=algorithm
    ).fit(self._counts_obs_norm)
    self.set_manifold(
        svd.transform(self._counts_obs_norm), svd.transform(self._counts_sim_norm)
    )


def pca(
    self: Scrublet,
    n_prin_comps: int = 50,
    *,
    random_state: AnyRandom = 0,
    svd_solver: Literal["auto", "full", "arpack", "randomized"] = "arpack",
) -> None:
    if self._counts_sim_norm is None:
        raise RuntimeError("_counts_sim_norm is not set")
    from sklearn.decomposition import PCA

    X_obs = self._counts_obs_norm.toarray()
    X_sim = self._counts_sim_norm.toarray()

    pca = PCA(
        n_components=n_prin_comps, random_state=random_state, svd_solver=svd_solver
    ).fit(X_obs)
    self.set_manifold(pca.transform(X_obs), pca.transform(X_sim))
