from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Literal, Protocol

import numpy as np

if TYPE_CHECKING:
    from typing import Any, Self

    from .._compat import CSRBase


# These two are used with get_literal_vals elsewhere
_Method = Literal["umap", "gauss"]
_KnownTransformer = Literal["pynndescent", "sklearn", "rapids"]

# sphinx-autodoc-typehints can’t transitively import types from if TYPE_CHECKING blocks,
# so these four needs to be importable

_MetricFn = Callable[[np.ndarray, np.ndarray], float]
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
_Metric = _MetricSparseCapable | _MetricScipySpatial


class KnnTransformerLike(Protocol):
    """See :class:`~sklearn.neighbors.KNeighborsTransformer`."""

    def fit(self, X, y: None = None): ...
    def transform(self, X) -> CSRBase: ...

    # from TransformerMixin
    def fit_transform(self, X, y: None = None) -> CSRBase: ...

    # from BaseEstimator
    def get_params(self, *, deep: bool = True) -> dict[str, Any]: ...
    def set_params(self, **params: Any) -> Self: ...
