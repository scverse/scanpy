from __future__ import annotations

from typing import Any, Literal

import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix
from sklearn.base import check_is_fitted
from sklearn.neighbors import KNeighborsTransformer

from ... import settings


class RapidsKNNTransformer(KNeighborsTransformer):
    """Compute nearest neighbors using RAPIDS cuml.

    See :class:`cuml.neighbors.NearestNeighbors`.
    """

    def __init__(self, *args, handle=None, **kwargs) -> None:
        from cuml.neighbors import NearestNeighbors

        super().__init__(*args, **kwargs)
        self.nn = NearestNeighbors(
            n_neighbors=self.n_neighbors,
            # https://docs.rapids.ai/api/cuml/nightly/api/#verbosity-levels
            verbose=settings.verbosity + 2,
            handle=handle,
            metric=self.metric,
            p=self.p,
            # algo_params=???,
            metric_params=self.metric_params,
            n_jobs=self.n_job,
        )

    def fit(self, X: np.ndarray, y: Any = None) -> RapidsKNNTransformer:
        """
        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            The data to compute nearest neighbors for.
        y
            Not used, present for API consistency by convention.
        """

        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        self.nn.fit(X_contiguous)

    def transform(self, X: np.ndarray) -> csr_matrix:
        check_is_fitted(self)
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        return self.nn.kneighbors_graph(X_contiguous)

    def fit_transform(self, X: np.ndarray, y: Any = None) -> csr_matrix:
        return self.fit(X, y).transform(X)

    def kneighbors(
        self,
        X: np.ndarray | None = None,
        n_neighbors: int | None = None,
        return_distance: bool = True,
    ) -> ndarray | tuple[ndarray, ndarray]:
        return self.nn.kneighbors(X, n_neighbors, return_distance)

    def kneighbors_graph(
        self,
        X: np.ndarray | None = None,
        n_neighbors: int | None = None,
        mode: Literal['distance', 'connectivity'] = 'connectivity',
    ):
        return self.nn.kneighbors_graph(X, n_neighbors, mode)
