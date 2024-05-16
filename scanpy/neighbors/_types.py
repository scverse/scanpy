from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Self, Union

    import numpy as np
    from scipy.sparse import spmatrix

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
    _Metric = Union[_MetricSparseCapable, _MetricScipySpatial]

# These two are used with get_Args elsewhere
_Method = Literal["umap", "gauss"]
_KnownTransformer = Literal["pynndescent", "sklearn", "rapids"]


class KnnTransformerLike(Protocol):
    """See :class:`~sklearn.neighbors.KNeighborsTransformer`."""

    def fit(self, X, y: None = None): ...

    def transform(self, X) -> spmatrix: ...

    # from TransformerMixin
    def fit_transform(self, X, y: None = None) -> spmatrix: ...

    # from BaseEstimator
    def get_params(self, deep: bool = True) -> dict[str, Any]: ...

    def set_params(self, **params: Any) -> Self: ...
