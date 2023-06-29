from __future__ import annotations

from typing import Any, Literal, Mapping

import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn_ann.utils import TransformerChecksMixin

from ... import settings


_Algorithm = Literal['auto', 'rbc', 'brute', 'ivfflat', 'ivfpq']
_Metric = Literal[
    'l1',
    'cityblock',
    'taxicab',
    'manhattan',
    'euclidean',
    'l2',
    'braycurtis',
    'canberra',
    'minkowski',
    'chebyshev',
    'jensenshannon',
    'cosine',
    'correlation',
]


class RapidsKNNTransformer(TransformerChecksMixin, TransformerMixin, BaseEstimator):
    """Compute nearest neighbors using RAPIDS cuml.

    See :class:`cuml.neighbors.NearestNeighbors`.
    """

    def __init__(
        self,
        *,
        handle=None,
        algorithm: _Algorithm = 'auto',
        n_neighbors: int,
        metric: _Metric = "euclidean",
        p: int = 2,
        algo_params: Mapping[str, Any] | None = None,
        metric_params: Mapping[str, Any] | None = None,
    ) -> None:
        from cuml.neighbors import NearestNeighbors

        self.n_neighbors = n_neighbors
        self.metric = metric
        self.p = p
        self.nn = NearestNeighbors(
            n_neighbors=n_neighbors,
            # https://docs.rapids.ai/api/cuml/nightly/api/#verbosity-levels
            verbose=settings.verbosity + 2,
            handle=handle,
            algorithm=algorithm,
            metric=metric,
            p=p,
            algo_params=algo_params,
            metric_params=metric_params,
            # TODO: output_type=...,
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
        return self

    def transform(self, X: np.ndarray) -> csr_matrix:
        self._transform_checks(X)
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        return self.nn.kneighbors_graph(X_contiguous)

    def fit_transform(self, X: np.ndarray, y: Any = None) -> csr_matrix:
        return self.fit(X, y).transform(X)

    def _more_tags(self):
        return {
            "requires_y": False,
            # TODO: are these correct?
            "preserves_dtype": [np.float32],
            "non_deterministic": True,
        }
