from __future__ import annotations

import warnings
from types import MappingProxyType
from typing import Any, Union, Mapping, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, coo_matrix
from sklearn.utils import check_random_state
from sklearn.base import BaseEstimator, TransformerMixin
from pynndescent import NNDescent

from ... import settings
from ..._utils import AnyRandom
from .._enums import _Metric, _MetricFn
from .._common import _get_sparse_matrix_from_indices_distances
from ._common import TransformerChecksMixin, mappings


_Backend = _Algorithm = Literal['umap']

_backends: dict[_Backend, set[_Algorithm]] = {
    'umap': {'umap'},
}

BACKENDS, ALGORITHMS = mappings(_backends)


class UMAPKNNTransformer(TransformerChecksMixin, TransformerMixin, BaseEstimator):
    """This is from umap.fuzzy_simplicial_set [McInnes18]_.
    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        The data to be modelled as a fuzzy simplicial set.
    n_neighbors
        The number of neighbors to use to approximate geodesic distance.
        Larger numbers induce more global estimates of the manifold that can
        miss finer detail, while smaller values will focus on fine manifold
        structure to the detriment of the larger picture.
    random_state
        A state capable being used as a numpy random state.
    metric
        The metric to use to compute distances in high dimensional space.
        If a string is passed it must match a valid predefined metric. If
        a general metric is required a function that takes two 1d arrays and
        returns a float can be provided. For performance purposes it is
        required that this be a numba jit'd function. Valid string metrics
        include:
            * euclidean
            * manhattan
            * chebyshev
            * minkowski
            * canberra
            * braycurtis
            * mahalanobis
            * wminkowski
            * seuclidean
            * cosine
            * correlation
            * haversine
            * hamming
            * jaccard
            * dice
            * russelrao
            * kulsinski
            * rogerstanimoto
            * sokalmichener
            * sokalsneath
            * yule
        Metrics that take arguments (such as minkowski, mahalanobis etc.)
        can have arguments passed via the metric_kwds dictionary. At this
        time care must be taken and dictionary elements must be ordered
        appropriately; this will hopefully be fixed in the future.
    metric_kwds
        Arguments to pass on to the metric, such as the ``p`` value for
        Minkowski distance.
    angular
        Whether to use angular/cosine distance for the random projection
        forest for seeding NN-descent to determine approximate nearest
        neighbors.
    verbose
        Whether to report information on the current progress of the algorithm.
    Returns
    -------
    **knn_indices**, **knn_dists** : np.arrays of shape (n_observations, n_neighbors)
    """

    # fit attributes are supposed to only exist once fitted
    fit_X_obs_: int
    knn_indices_: NDArray[np.int32]
    knn_dists_: NDArray[np.float32]
    forest_: NNDescent | None

    def __init__(
        self,
        n_neighbors: int,
        random_state: AnyRandom = None,
        metric: Union[_Metric, _MetricFn] = 'euclidean',
        metric_kwds: Mapping[str, Any] = MappingProxyType({}),
        angular: bool = False,
        verbose: bool = False,
    ) -> None:
        self.n_neighbors = n_neighbors
        self.random_state = check_random_state(random_state)
        self.metric = metric
        self.metric_kwds = metric_kwds
        self.angular = angular
        self.verbose = verbose

    def fit(self, X: np.ndarray | csr_matrix, y: Any = None) -> UMAPKNNTransformer:
        """
        Calculate approximate knn index using PyNNDescent.

        Since the PyNNDescent index is very close to a knn graph,
        X’s knn graph already gets cached at this point.
        """
        with warnings.catch_warnings():
            # umap 0.5.0
            warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
            from umap.umap_ import nearest_neighbors

        self.fit_X_obs_ = X.shape[0]
        self.knn_indices_, self.knn_dists_, self.forest_ = nearest_neighbors(
            X,
            self.n_neighbors,
            random_state=self.random_state,
            metric=self.metric,
            metric_kwds=self.metric_kwds,
            angular=self.angular,
            verbose=self.verbose,
            n_jobs=settings.n_jobs,
        )

        return self

    def transform(self, X: np.ndarray | csr_matrix | None = None) -> csr_matrix:
        """Perform a knn search on the index.

        Parameters
        ----------
        X
            Data to search knn for.
            If None, simply return the knn graph created as part of indexing.
        """
        self._transform_checks(X=X)
        if X is None:
            return _get_sparse_matrix_from_indices_distances(
                self.knn_indices_,
                self.knn_dists_,
                n_obs=self.fit_X_obs_,
                n_neighbors=self.n_neighbors,
            )
        else:
            raise NotImplementedError("")  # TODO

    def fit_transform(self, X: np.ndarray | csr_matrix, y: Any = None) -> csr_matrix:
        """Create a knn graph for X."""
        return self.fit(X, y).transform()

    def _more_tags(self):
        return {
            "requires_y": False,
            # TODO: is this correct?
            "preserves_dtype": [np.float32],
        }


def compute_connectivities(
    knn_indices: NDArray[np.int32],
    knn_dists: NDArray[np.float32],
    *,
    n_obs: int,
    n_neighbors: int,
    set_op_mix_ratio: float = 1.0,
    local_connectivity: float = 1.0,
) -> csr_matrix:
    """\
    This is from umap.fuzzy_simplicial_set [McInnes18]_.

    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        from umap.umap_ import fuzzy_simplicial_set

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(
        X,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if isinstance(connectivities, tuple):
        # In umap-learn 0.4, this returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    return connectivities.tocsr()
