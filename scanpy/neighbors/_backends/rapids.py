from __future__ import annotations

from typing import Any, Literal, Mapping

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.base import BaseEstimator, TransformerMixin, check_is_fitted
from sklearn.exceptions import NotFittedError

from ... import settings
from ._common import TransformerChecksMixin, mappings

_Backend = Literal['rapids']
_Algorithm = Literal['rbc', 'brute', 'ivfflat', 'ivfpq']
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

_backends: dict[_Backend, set[_Algorithm]] = {
    'rapids': {'rbc', 'brute', 'ivfflat', 'ivfpq'},
}

BACKENDS, ALGORITHMS = mappings(_backends)


CudaArrayLike = np.ndarray  # TODO


class RapidsKNNTransformer(TransformerChecksMixin, TransformerMixin, BaseEstimator):
    """Compute nearest neighbors using RAPIDS cuml.

    See :class:`cuml.neighbors.NearestNeighbors`.
    """

    def __init__(
        self,
        *,
        handle=None,
        algorithm: _Algorithm | Literal['auto'] = 'auto',
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

    def __sklearn_is_fitted__(self) -> bool:
        try:
            check_is_fitted(self.nn)
        except NotFittedError:
            return False
        else:
            return True

    def fit(self, X: CudaArrayLike, y: Any = None) -> RapidsKNNTransformer:
        """Index data for knn search."""
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        self.nn.fit(X_contiguous)
        return self

    def transform(self, X: CudaArrayLike) -> csr_matrix:
        """Perform knn search on the index."""
        self._transform_checks(X)
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        return self.nn.kneighbors_graph(X_contiguous)

    def _more_tags(self):
        return {
            "requires_y": False,
            # TODO: are these correct?
            "preserves_dtype": [np.float32],
            "non_deterministic": True,
        }
