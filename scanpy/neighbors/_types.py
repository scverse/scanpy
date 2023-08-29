from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, Protocol, Union as _U
from collections.abc import Callable

import numpy as np
from scipy.sparse import spmatrix

if TYPE_CHECKING:
    from typing import Self


_Method = Literal['umap', 'gauss']

_KnownTransformer = Literal['pynndescent', 'rapids']

_MetricFn = Callable[[np.ndarray, np.ndarray], float]
# from sklearn.metrics.pairwise_distances.__doc__:
_MetricSparseCapable = Literal[
    'cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'
]
_MetricScipySpatial = Literal[
    'braycurtis',
    'canberra',
    'chebyshev',
    'correlation',
    'dice',
    'hamming',
    'jaccard',
    'kulsinski',
    'mahalanobis',
    'minkowski',
    'rogerstanimoto',
    'russellrao',
    'seuclidean',
    'sokalmichener',
    'sokalsneath',
    'sqeuclidean',
    'yule',
]
_Metric = _U[_MetricSparseCapable, _MetricScipySpatial]


class KnnTransformerLike(Protocol):
    """See :class:`~sklearn.neighbors.KNeighborsTransformer`."""

    def fit(self, X, y: None = None):
        ...

    def transform(self, X) -> spmatrix:
        ...

    # from TransformerMixin
    def fit_transform(self, X, y: None = None) -> spmatrix:
        ...

    # from BaseEstimator
    def get_params(self, deep: bool = True) -> dict[str, Any]:
        ...

    def set_params(self, **params: Any) -> Self:
        ...
