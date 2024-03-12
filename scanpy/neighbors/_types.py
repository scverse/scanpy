from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, Protocol
from typing import Callable as _C
from typing import Union as _U

import numpy as np

if TYPE_CHECKING:
    from typing import Self

    from scipy.sparse import spmatrix


_Method = Literal["umap", "gauss"]

_KnownTransformer = Literal["pynndescent", "sklearn", "rapids"]

_MetricFn = _C[[np.ndarray, np.ndarray], float]
# from sklearn.metrics.pairwise_distances.__doc__:
_MetricSparseCapable = Literal[
    "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan"
]
_MetricScipySpatial = Literal[
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
]
_Metric = _U[_MetricSparseCapable, _MetricScipySpatial]


class KnnTransformerLike(Protocol):
    """See :class:`~sklearn.neighbors.KNeighborsTransformer`."""

    def fit(self, X, y: None = None): ...

    def transform(self, X) -> spmatrix: ...

    # from TransformerMixin
    def fit_transform(self, X, y: None = None) -> spmatrix: ...

    # from BaseEstimator
    def get_params(self, deep: bool = True) -> dict[str, Any]: ...

    def set_params(self, **params: Any) -> Self: ...
