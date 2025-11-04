from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Literal, Protocol

import numpy as np

if TYPE_CHECKING:
    from typing import Any, Self

    from .._compat import CSRBase

__all__ = [
    "KnnTransformerLike",
    "_KnownTransformer",
    "_Method",
    "_Metric",
    "_MetricFn",
    "_MetricScipySpatial",
    "_MetricSparseCapable",
]


type _Method = Literal["umap", "gauss"]
type _KnownTransformer = Literal["pynndescent", "sklearn", "rapids"]

type _MetricFn = Callable[[np.ndarray, np.ndarray], float]
# from sklearn.metrics.pairwise_distances.__doc__:
type _MetricSparseCapable = Literal[
    "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan"
]
type _MetricScipySpatial = Literal[
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
type _Metric = _MetricSparseCapable | _MetricScipySpatial


class KnnTransformerLike(Protocol):
    """See :class:`~sklearn.neighbors.KNeighborsTransformer`."""

    def fit(self, x, /, y: None = None): ...
    def transform(self, x, /) -> CSRBase: ...

    # from TransformerMixin
    def fit_transform(self, x, /, y: None = None) -> CSRBase: ...

    # from BaseEstimator
    def get_params(self, *, deep: bool = True) -> dict[str, Any]: ...
    def set_params(self, **params: Any) -> Self: ...
