from typing import Union

import numpy as np

from ..enums import _Metric, _MetricFn


def compute_neighbors_rapids(
    X: np.ndarray, n_neighbors: int, metric: Union[_Metric, _MetricFn] = 'euclidean'
):
    """Compute nearest neighbors using RAPIDS cuml.
    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        The data to compute nearest neighbors for.
    n_neighbors
        The number of neighbors to use.
    metric
        The metric to use to compute distances in high dimensional space.
        This string must match a valid predefined metric in RAPIDS cuml.
        Returns
    -------
    **knn_indices**, **knn_dists** : np.arrays of shape (n_observations, n_neighbors)
    """
    from cuml.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
    X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
    nn.fit(X_contiguous)
    knn_dist, knn_indices = nn.kneighbors(X_contiguous)
    return knn_indices, knn_dist
