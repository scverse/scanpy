from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Literal, Protocol, TypedDict

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, NotRequired, Self, TypeAlias

    from .._compat import CSRBase
    from .._utils.random import _LegacyRandom

    # TODO: make `type` when https://github.com/sphinx-doc/sphinx/pull/13508 is released
    RPForestDict: TypeAlias = Mapping[str, Mapping[str, np.ndarray]]  # noqa: UP040

__all__ = [
    "KnnTransformerLike",
    "KwdsForTransformer",
    "NeighborsDict",
    "NeighborsParams",
    "_KnownTransformer",
    "_Method",
    "_Metric",
    "_MetricFn",
    "_MetricScipySpatial",
    "_MetricSparseCapable",
]

type _Method = Literal["umap", "gauss", "jaccard"]
type _KnownTransformer = Literal["pynndescent", "sklearn"]

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


class KwdsForTransformer(TypedDict):
    """Keyword arguments passed to a _KnownTransformer.

    IMPORTANT: when changing the parameters set here,
    update the “*ignored*” part in the parameter docs!
    """

    n_neighbors: int
    metric: _Metric | _MetricFn
    metric_params: Mapping[str, Any]
    rng: NotRequired[np.random.Generator]


class NeighborsDict(TypedDict):
    connectivities_key: str
    distances_key: str
    params: NeighborsParams
    rp_forest: NotRequired[RPForestDict]


class NeighborsParams(TypedDict):
    n_neighbors: int
    method: _Method
    metric: _Metric | _MetricFn | None
    random_state: NotRequired[_LegacyRandom]
    metric_kwds: NotRequired[Mapping[str, Any]]
    use_rep: NotRequired[str | list[str]]  # see `scanpy.get.get._rep_to_json`
    n_pcs: NotRequired[int]
