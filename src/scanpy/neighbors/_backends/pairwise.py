from __future__ import annotations

from dataclasses import KW_ONLY, dataclass
from typing import TYPE_CHECKING

from sklearn.base import TransformerMixin

from .._common import (
    _get_indices_distances_from_dense_matrix,
    _get_sparse_matrix_from_indices_distances,
)

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Literal, Self

    from numpy.typing import NDArray

    from ..._utils import _CSMatrix

    _Metric = Literal["cityblock", "cosine", "euclidean", "l1", "l2", "manhattan"]
    _MatrixLike = NDArray | _CSMatrix


@dataclass
class PairwiseDistancesTransformer(TransformerMixin):
    _: KW_ONLY
    algorithm: Literal["brute"]
    n_jobs: int
    n_neighbors: int
    metric: _Metric
    metric_params: Mapping[str, object]

    def fit(self, x: _MatrixLike) -> Self:
        self.x_ = x
        return self

    def transform(self, y: _MatrixLike | None) -> _CSMatrix:
        from sklearn.metrics import pairwise_distances

        d_arr = pairwise_distances(self.x_, y, metric=self.metric, **self.metric_params)
        ind, dist = _get_indices_distances_from_dense_matrix(d_arr, self.n_neighbors)
        return _get_sparse_matrix_from_indices_distances(ind, dist, keep_self=True)
