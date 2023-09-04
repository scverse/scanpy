from __future__ import annotations

from typing import Any, Literal
from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import csr_matrix
from sklearn.base import BaseEstimator, TransformerMixin, check_is_fitted
from sklearn.exceptions import NotFittedError

from ... import settings
from ._common import TransformerChecksMixin

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
            output_type='input',  # could also be None to respect global setting
        )

    def __sklearn_is_fitted__(self) -> bool:
        try:
            check_is_fitted(self.nn)
        except NotFittedError:
            return False
        else:
            return True

    def fit(self, X: ArrayLike, y: Any = None) -> RapidsKNNTransformer:
        """Index data for knn search."""
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        self.nn.fit(X_contiguous)
        return self

    def transform(self, X: ArrayLike) -> csr_matrix:
        """Perform knn search on the index."""
        self._transform_checks(X)
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        return self.nn.kneighbors_graph(X_contiguous)

    def _more_tags(self) -> dict[str, Any]:
        """See :label:`sklearn:estimator_tags`"""
        return {
            "requires_y": False,
            "preserves_dtype": [np.float32],
            "non_deterministic": True,
        }
