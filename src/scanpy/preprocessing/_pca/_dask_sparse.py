from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, cast, overload

import numpy as np
import scipy.linalg
from numpy.typing import NDArray

from .._utils import _get_mean_var

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike
    from scipy import sparse

    from ..._compat import DaskArray

    CSMatrix = sparse.csr_matrix | sparse.csc_matrix


@dataclass
class PCASparseDask:
    n_components: int | None = None

    def fit(self, x: DaskArray) -> PCASparseFit:
        # this method makes `self` into the fitted version
        self.__class__ = PCASparseFit
        self = cast(PCASparseFit, self)

        assert isinstance(x.shape, tuple)
        self.n_components_ = (
            min(x.shape) if self.n_components is None else self.n_components
        )
        self.n_samples_ = x.shape[0]
        self.n_features_in_ = x.shape[1] if x.ndim == 2 else 1
        self.dtype_ = x.dtype
        covariance, self.mean_ = _cov_sparse_dask(x)
        self.explained_variance_, self.components_ = scipy.linalg.eigh(
            covariance, lower=False
        )

        # Arrange eigenvectors and eigenvalues in descending order
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

    def fit_transform(self, x: DaskArray, y: DaskArray | None = None) -> DaskArray:
        if y is None:
            y = x
        return self.fit(x).transform(y)


@dataclass
class PCASparseFit(PCASparseDask):
    n_components_: int = field(init=False)
    n_samples_: int = field(init=False)
    n_features_in_: int = field(init=False)
    dtype_: np.dtype = field(init=False)
    mean_: NDArray[np.floating] = field(init=False)
    components_: NDArray[np.floating] = field(init=False)
    explained_variance_: NDArray[np.floating] = field(init=False)
    explained_variance_ratio_: NDArray[np.floating] = field(init=False)
    noise_variance_: NDArray[np.floating] = field(init=False)

    def transform(self, x: DaskArray) -> DaskArray:
        if TYPE_CHECKING:
            # The type checker does not understand imports from dask.array
            import dask.array.core as da
        else:
            import dask.array as da

        def transform_block(
            x_part: CSMatrix,
            mean_: NDArray[np.floating],
            components_: NDArray[np.floating],
        ):
            pre_mean = mean_ @ components_.T
            mean_impact = np.ones((x_part.shape[0], 1)) @ pre_mean.reshape(1, -1)
            return (x_part @ components_.T) - mean_impact

        return da.map_blocks(
            transform_block,
            x,
            mean_=self.mean_,
            components_=self.components_,
            chunks=(x.chunks[0], self.n_components_),
            meta=np.zeros([0], dtype=x.dtype),
            dtype=x.dtype,
        )


@overload
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: Literal[False] = False, dtype: DTypeLike | None = None
) -> tuple[NDArray[np.floating], NDArray[np.floating]]: ...
@overload
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: Literal[True], dtype: DTypeLike | None = None
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]: ...
def _cov_sparse_dask(
    x: DaskArray, *, return_gram: bool = False, dtype: DTypeLike | None = None
) -> (
    tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
    | tuple[NDArray[np.floating], NDArray[np.floating]]
):
    """\
    Computes the covariance matrix and row/col means of matrix `x`.

    Parameters
    ----------

    x
        A sparse matrix
    return_gram
        If `True`, the gram matrix will be returned and a copy will be created
        to store the results of the covariance,
        while if `False`, the local gram matrix result will be overwritten.
        (only used for unit testing at the moment)
    dtype
        The data type of the result (excluding the means)

    Returns
    -------

    :math:`\\cov(X, X)`
        The covariance matrix of `x` in the form :math:`\\cov(X, X) = \\E(XX) - \\E(X)\\E(X)`.
    :math:`\\gram(X, X)`
        When return_gram is `True`, the gram matrix of `x` in the form :math:`\\frac{1}{n} X.T \\dot X`.
    :math:`\\mean(X)`
        The row means of `x`.
    """
    if TYPE_CHECKING:
        import dask.array.core as da
        import dask.base as dask
    else:
        import dask
        import dask.array as da

    if dtype is None:
        dtype = np.float64 if np.issubdtype(x.dtype, np.integer) else x.dtype
    else:
        dtype = np.dtype(dtype)

    def gram_block(x_part: CSMatrix):
        gram_matrix: CSMatrix = x_part.T @ x_part
        return gram_matrix.toarray()[None, ...]  # need new axis for summing

    gram_matrix_dask: DaskArray = da.map_blocks(
        gram_block,
        x,
        new_axis=(1,),
        chunks=((1,) * x.blocks.size, (x.shape[1],), (x.shape[1],)),
        meta=np.array([], dtype=x.dtype),
        dtype=x.dtype,
    ).sum(axis=0)
    mean_x_dask, _ = _get_mean_var(x)
    gram_matrix, mean_x = cast(
        tuple[NDArray, NDArray[np.float64]],
        dask.compute(gram_matrix_dask, mean_x_dask),
    )
    gram_matrix = gram_matrix.astype(dtype)
    gram_matrix /= x.shape[0]

    cov_result = gram_matrix.copy() if return_gram else gram_matrix
    cov_result -= mean_x[:, None] @ mean_x[None, :]

    if return_gram:
        return cov_result, gram_matrix, mean_x
    return cov_result, mean_x
