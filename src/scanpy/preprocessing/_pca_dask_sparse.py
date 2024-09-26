from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, cast, overload

import numpy as np
import scipy.linalg
from numpy.typing import NDArray

from ._utils import _get_mean_var

if TYPE_CHECKING:
    from typing import Literal

    from scipy import sparse

    from .._compat import DaskArray

    CSMatrix = sparse.csr_matrix | sparse.csc_matrix


class PCASparseFit(Protocol):
    _n_components: int

    n_components_: int
    n_samples_: int
    n_features_in_: int
    dtype_: np.dtype
    mean_: NDArray[np.floating]
    components_: NDArray[np.floating]
    explained_variance_: NDArray[np.floating]
    explained_variance_ratio_: NDArray[np.floating]
    noise_variance_: NDArray[np.floating]

    def fit(self, x: DaskArray) -> PCASparseFit: ...
    def transform(self, X: DaskArray) -> DaskArray: ...


@dataclass
class PCASparseDask:
    _n_components: int | None = None

    def fit(self, x: DaskArray) -> PCASparseFit:
        self = cast(PCASparseFit, self)  # this makes `self` into the fitted version
        assert isinstance(x.shape, tuple)
        self.n_components_ = (
            min(x.shape[:2]) if self._n_components is None else self._n_components
        )
        self.n_samples_ = x.shape[0]
        self.n_features_in_ = x.shape[1] if x.ndim == 2 else 1
        self.dtype_ = x.dtype
        covariance, self.mean_ = _cov_sparse_dask(x)
        self.explained_variance_, self.components_ = scipy.linalg.eigh(
            covariance, lower=False
        )
        # NOTE: We reverse the eigen vector and eigen values here
        # because cupy provides them in ascending order. Make a copy otherwise
        # it is not C_CONTIGUOUS anymore and would error when converting to
        # CumlArray
        self.explained_variance_ = self.explained_variance_[::-1]

        self.components_ = np.flip(self.components_, axis=1)

        self.components_ = self.components_.T[: self.n_components_, :]

        self.explained_variance_ratio_ = self.explained_variance_ / np.sum(
            self.explained_variance_
        )
        if self.n_components_ < min(self.n_samples_, self.n_features_in_):
            self.noise_variance_ = self.explained_variance_[self.n_components_ :].mean()
        else:
            self.noise_variance_ = np.array([0.0])
        self.explained_variance_ = self.explained_variance_[: self.n_components_]

        self.explained_variance_ratio_ = self.explained_variance_ratio_[
            : self.n_components_
        ]
        return self

    def transform(self: PCASparseFit, x: DaskArray) -> DaskArray:
        def _transform(X_part, mean_, components_):
            pre_mean = mean_ @ components_.T
            mean_impact = np.ones((X_part.shape[0], 1)) @ pre_mean.reshape(1, -1)
            X_transformed = X_part.dot(components_.T) - mean_impact
            return X_transformed

        X_pca = x.map_blocks(
            _transform,
            mean_=self.mean_,
            components_=self.components_,
            dtype=x.dtype,
            chunks=(x.chunks[0], self.n_components_),
            meta=np.zeros([0], dtype=x.dtype),
        )

        return X_pca

    def fit_transform(self, x: DaskArray, y: DaskArray | None = None) -> DaskArray:
        if y is None:
            y = x
        return self.fit(x).transform(y)


@overload
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: Literal[False] = False
) -> tuple[NDArray[np.floating], NDArray[np.floating]]: ...
@overload
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: Literal[True]
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]: ...
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: bool = False
) -> (
    tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
    | tuple[NDArray[np.floating], NDArray[np.floating]]
):
    """
    Computes the mean and the covariance of matrix X of
    the form Cov(X, X) = E(XX) - E(X)E(X)

    This is a temporary fix for
    cuml issue #5475 and cupy issue #7699,
    where the operation `x.T.dot(x)` did not work for
    larger sparse matrices.

    Parameters
    ----------

    x : cupyx.scipy.sparse of size (m, n)
    return_gram : boolean (default = False)
        If True, gram matrix of the form (1 / n) * X.T.dot(X)
        will be returned.
        When True, a copy will be created
        to store the results of the covariance.
        When False, the local gram matrix result
        will be overwritten
    return_mean: boolean (default = False)
        If True, the Maximum Likelihood Estimate used to
        calculate the mean of X and X will be returned,
        of the form (1 / n) * mean(X) and (1 / n) * mean(X)

    Returns
    -------

    result : cov(X, X) when return_gram and return_mean are False
            cov(X, X), gram(X, X) when return_gram is True,
            return_mean is False
            cov(X, X), mean(X), mean(X) when return_gram is False,
            return_mean is True
            cov(X, X), gram(X, X), mean(X), mean(X)
            when return_gram is True and return_mean is True
    """
    import dask

    from ._kernels._pca_sparse_kernel import (
        _copy_kernel,
        _cov_kernel,
        _gramm_kernel_csr,
    )

    compute_mean_cov = _gramm_kernel_csr(x.dtype)
    compute_mean_cov.compile()
    n_cols = x.shape[1]

    def gram_block(x_part: CSMatrix):
        gram_matrix: CSMatrix = x_part.T @ x_part
        return gram_matrix.toarray()[None, ...]  # need new axis for summing

    n_blocks = len(x.to_delayed().ravel())
    gram_matrix = x.map_blocks(
        gram_block,
        new_axis=(1,),
        chunks=((1,) * n_blocks, (x.shape[1],), (x.shape[1],)),
        meta=np.array([]),
        dtype=x.dtype,
    ).sum(axis=0)
    mean_x, _ = _get_mean_var(x)
    gram_matrix, mean_x = cast(
        tuple[NDArray[np.floating], NDArray[np.floating]],
        dask.compute(gram_matrix, mean_x),
    )
    mean_x = mean_x.astype(x.dtype)
    copy_gram = _copy_kernel(x.dtype)
    block = (32, 32)
    grid = (math.ceil(n_cols / block[0]), math.ceil(n_cols / block[1]))
    copy_gram(
        grid,
        block,
        (gram_matrix, n_cols),
    )

    gram_matrix *= 1 / x.shape[0]

    if return_gram:
        cov_result = np.zeros(
            (gram_matrix.shape[0], gram_matrix.shape[0]),
            dtype=gram_matrix.dtype,
        )
    else:
        cov_result = gram_matrix

    compute_cov = _cov_kernel(gram_matrix.dtype)

    block_size = (32, 32)
    grid_size = (math.ceil(gram_matrix.shape[0] / 8),) * 2
    compute_cov(
        grid_size,
        block_size,
        (cov_result, gram_matrix, mean_x, mean_x, gram_matrix.shape[0]),
    )

    if return_gram:
        return cov_result, gram_matrix, mean_x
    return cov_result, mean_x
