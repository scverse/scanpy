from __future__ import annotations

from typing import Literal, Union
from collections.abc import Callable

import numpy as np


_Method = Literal['umap', 'gauss', 'rapids']
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
_Metric = Union[_MetricSparseCapable, _MetricScipySpatial]
