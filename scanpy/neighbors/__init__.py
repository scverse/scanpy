import numpy as np
import scipy
from scipy.sparse import issparse
from scipy.sparse import coo_matrix
from sklearn.metrics import pairwise_distances
from .. import settings
from .. import logging as logg
from .. utils import doc_params
from .. import utils
from ..tools._utils import choose_representation, doc_use_rep, doc_n_pcs

N_DCS = 15  # default number of diffusion components
N_PCS = 50  # default number of PCs


@doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def neighbors(
        adata,
        n_neighbors=15,
        n_pcs=None,
        use_rep=None,
        knn=True,
        random_state=0,
        method='umap',
        metric='euclidean',
        metric_kwds={},
        copy=False):
    """\
    Compute a neighborhood graph of observations [McInnes18]_.

    The neighbor search efficiency of this heavily relies on UMAP [McInnes18]_,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='diffmap'`,
    connectivities are computed according to [Coifman05]_, in the adaption of
    [Haghverdi16]_.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    n_neighbors : `int`, optional (default: 15)
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    {n_pcs}
    {use_rep}
    knn : `bool`, optional (default: `True`)
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    method : {{'umap', 'gauss', `None`}}  (default: `'umap'`)
        Use 'umap' [McInnes18]_ or 'gauss' (Gauss kernel following [Coifman05]_
        with adaptive width [Haghverdi16]_) for computing connectivities.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:
    connectivities : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    distances : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """
    logg.info('computing neighbors', r=True)
    adata = adata.copy() if copy else adata
    if adata.isview:  # we shouldn't need this here...
        adata._init_as_actual(adata.copy())
    neighbors = Neighbors(adata)
    neighbors.compute_neighbors(
        n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs, use_rep=use_rep,
        method=method, metric=metric, metric_kwds=metric_kwds,
        random_state=random_state)
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': method}
    adata.uns['neighbors']['distances'] = neighbors.distances
    adata.uns['neighbors']['connectivities'] = neighbors.connectivities
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns[\'neighbors\']`\n'
        '    \'distances\', weighted adjacency matrix\n'
        '    \'connectivities\', weighted adjacency matrix')
    return adata if copy else None


def compute_euclidean_distances_using_matrix_mult(X, Y):
    """Compute euclidean distance matrix for data arrays X and Y.

    Parameters
    ----------
    X : np.ndarray
        Data array (rows store observations, columns store variables).
    Y : np.ndarray
        Data array (rows store observations, columns store variables).

    Returns
    -------
    distances : np.ndarray
        Distance matrix.
    """
    XX = np.einsum('ij,ij->i', X, X)[:, None]
    if X is Y:
        YY = XX
    else:
        YY = np.einsum('ij,ij->i', Y, Y)[:, None]
    distances = np.dot(X, Y.T)
    distances *= -2
    distances += XX
    distances += YY.T
    np.maximum(distances, 0, out=distances)
    if X is Y:
        distances.flat[::distances.shape[0] + 1] = 0.
    # print(distances)
    # print(distances[distances < 0])
    # distances[distances < 0] = 0  # set to 0
    distances = np.sqrt(distances)
    return distances


def compute_neighbors_numpy_chunk(X, Y, n_neighbors):
    """Compute distance matrix in euclidean norm for chunk.
    """
    D = compute_euclidean_distances_using_matrix_mult(X, Y)
    chunk_range = np.arange(D.shape[0])[:, None]
    indices_chunk = np.argpartition(D, n_neighbors-1, axis=1)[:, :n_neighbors]
    indices_chunk = indices_chunk[chunk_range,
                                  np.argsort(D[chunk_range, indices_chunk])]
    distances_chunk = D[chunk_range, indices_chunk]
    # we know that a point has zero-distance to itself:
    #     set the first distance to zero
    #     if we don't this, for large data (not treated by tests)
    #     we might introduce spurious small non-zero values
    #     which affects backwards compat
    distances_chunk[:, 0] = 0
    return indices_chunk, distances_chunk


def compute_neighbors_numpy(X, n_neighbors, knn=True):
    """Compute distance matrix in uared euclidean norm.
    """
    if not knn:
        D = compute_euclidean_distances_using_matrix_mult(X, X)
        indices, distances = get_indices_distances_from_dense_matrix(D, n_neighbors)
        return D, indices, distances
    # assume we can fit at max 20000 data points into memory
    len_chunk = np.ceil(min(20000, X.shape[0])).astype(int)
    n_chunks = np.ceil(X.shape[0] / len_chunk).astype(int)
    chunks = [np.arange(start, min(start + len_chunk, X.shape[0]))
             for start in range(0, n_chunks * len_chunk, len_chunk)]
    indices = np.zeros((X.shape[0], n_neighbors), dtype=int)
    distances = np.zeros((X.shape[0], n_neighbors), dtype=np.float32)
    for i_chunk, chunk in enumerate(chunks):
        indices_chunk, distances_chunk = compute_neighbors_numpy_chunk(
            X[chunk], X, n_neighbors)
        indices[chunk] = indices_chunk
        distances[chunk] = distances_chunk
    D = get_sparse_matrix_from_indices_distances_numpy(
        indices, distances, X.shape[0], n_neighbors)
    return D, indices, distances


def compute_neighbors_umap(
        X, n_neighbors, random_state=None,
        metric='euclidean', metric_kwds={}, angular=False,
        verbose=False):
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
    n_neighbors: int
        The number of neighbors to use to approximate geodesic distance.
        Larger numbers induce more global estimates of the manifold that can
        miss finer detail, while smaller values will focus on fine manifold
        structure to the detriment of the larger picture.
    random_state: numpy RandomState or equivalent
        A state capable being used as a numpy random state.
    metric: string or function (optional, default 'euclidean')
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
    metric_kwds: dict (optional, default {})
        Arguments to pass on to the metric, such as the ``p`` value for
        Minkowski distance.
    angular: bool (optional, default False)
        Whether to use angular/cosine distance for the random projection
        forest for seeding NN-descent to determine approximate nearest
        neighbors.
    verbose: bool (optional, default False)
        Whether to report information on the current progress of the algorithm.

    Returns
    -------
    knn_indices, knn_dists : np.arrays of shape (n_observations, n_neighbors)
    """
    from .umap import sparse
    from .umap.umap_ import rptree_leaf_array, make_nn_descent
    from .umap import distances as dist
    from .umap import sparse
    import scipy
    from sklearn.utils import check_random_state

    INT32_MIN = np.iinfo(np.int32).min + 1
    INT32_MAX = np.iinfo(np.int32).max - 1

    random_state = check_random_state(random_state)

    if metric == 'precomputed':
        # Note that this does not support sparse distance matrices yet ...
        # Compute indices of n nearest neighbors
        knn_indices = np.argsort(X)[:, :n_neighbors]
        # Compute the nearest neighbor distances
        #   (equivalent to np.sort(X)[:, :n_neighbors])
        knn_dists = X[np.arange(X.shape[0])[:, None], knn_indices].copy()
    else:
        if callable(metric):
            distance_func = metric
        elif metric in dist.named_distances:
            distance_func = dist.named_distances[metric]
        else:
            raise ValueError('Metric is neither callable, ' +
                             'nor a recognised string')

        if metric in ('cosine', 'correlation', 'dice', 'jaccard'):
            angular = True

        rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

        if scipy.sparse.isspmatrix_csr(X):
            if metric in sparse.sparse_named_distances:
                distance_func = sparse.sparse_named_distances[metric]
                if metric in sparse.sparse_need_n_features:
                    metric_kwds['n_features'] = X.shape[1]
            else:
                raise ValueError('Metric {} not supported for sparse ' +
                                'data'.format(metric))
            metric_nn_descent = sparse.make_sparse_nn_descent(
                distance_func, tuple(metric_kwds.values()))
            leaf_array = rptree_leaf_array(X, n_neighbors,
                                           rng_state, n_trees=10,
                                           angular=angular)
            knn_indices, knn_dists = metric_nn_descent(X.indices,
                                                       X.indptr,
                                                       X.data,
                                                       X.shape[0],
                                                       n_neighbors,
                                                       rng_state,
                                                       max_candidates=60,
                                                       rp_tree_init=True,
                                                       leaf_array=leaf_array,
                                                       verbose=verbose)
        else:
            metric_nn_descent = make_nn_descent(distance_func,
                                                tuple(metric_kwds.values()))
            # TODO: Hacked values for now
            n_trees = 5 + int(round((X.shape[0]) ** 0.5 / 20.0))
            n_iters = max(5, int(round(np.log2(X.shape[0]))))

            leaf_array = rptree_leaf_array(X, n_neighbors,
                                           rng_state, n_trees=n_trees,
                                           angular=angular)
            knn_indices, knn_dists = metric_nn_descent(X,
                                                       n_neighbors,
                                                       rng_state,
                                                       max_candidates=60,
                                                       rp_tree_init=True,
                                                       leaf_array=leaf_array,
                                                       n_iters=n_iters,
                                                       verbose=verbose)

        if np.any(knn_indices < 0):
            logg.warn('Failed to correctly find n_neighbors for some samples. '
                 'Results may be less than ideal. Try re-running with '
                 'different parameters.')

    return knn_indices, knn_dists


def compute_connectivities_umap(knn_indices, knn_dists,
        n_obs, n_neighbors, set_op_mix_ratio=1.0,
        local_connectivity=1.0, bandwidth=1.0):
    """This is from umap.fuzzy_simplicial_set [McInnes18]_.

    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    from .umap.umap_ import smooth_knn_dist

    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    sims = np.zeros((n_obs * n_neighbors), dtype=np.float64)
    dists = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    sigmas, rhos = smooth_knn_dist(knn_dists, n_neighbors,
                                   local_connectivity=local_connectivity)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                sim = 0.0
                dist = 0.0
            elif knn_dists[i, j] - rhos[i] <= 0.0:
                sim = 1.0
                dist = knn_dists[i, j]
            else:
                sim = np.exp(-((knn_dists[i, j] - rhos[i]) / (sigmas[i] *
                                                              bandwidth)))
                dist = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            sims[i * n_neighbors + j] = sim
            dists[i * n_neighbors + j] = dist

    connectivities = coo_matrix((sims, (rows, cols)),
                               shape=(n_obs, n_obs))
    connectivities.eliminate_zeros()

    distances = coo_matrix((dists, (rows, cols)),
                           shape=(n_obs, n_obs))
    distances.eliminate_zeros()

    transpose = connectivities.transpose()

    prod_matrix = connectivities.multiply(transpose)

    connectivities = set_op_mix_ratio * (connectivities + transpose - prod_matrix) + \
             (1.0 - set_op_mix_ratio) * prod_matrix

    connectivities.eliminate_zeros()
    return distances.tocsr(), connectivities.tocsr()


def get_sparse_matrix_from_indices_distances_umap(knn_indices, knn_dists, n_obs, n_neighbors):
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)),
                                      shape=(n_obs, n_obs))
    return result.tocsr()


def get_sparse_matrix_from_indices_distances_numpy(indices, distances, n_obs, n_neighbors):
    n_nonzero = n_obs * n_neighbors
    indptr = np.arange(0, n_nonzero + 1, n_neighbors)
    D = scipy.sparse.csr_matrix((distances.copy().ravel(),  # copy the data, otherwise strange behavior here
                                indices.copy().ravel(),
                                indptr),
                                shape=(n_obs, n_obs))
    D.eliminate_zeros()
    return D


def get_indices_distances_from_sparse_matrix(D, n_neighbors):
    indices = np.zeros((D.shape[0], n_neighbors), dtype=int)
    distances = np.zeros((D.shape[0], n_neighbors), dtype=D.dtype)
    n_neighbors_m1 = n_neighbors - 1
    for i in range(indices.shape[0]):
        neighbors = D[i].nonzero()  # 'true' and 'spurious' zeros
        indices[i, 0] = i
        distances[i, 0] = 0
        # account for the fact that there might be more than n_neighbors
        # due to an approximate search
        # [the point itself was not detected as its own neighbor during the search]
        if len(neighbors[1]) > n_neighbors_m1:
            sorted_indices = np.argsort(D[i][neighbors].A1)[:n_neighbors_m1]
            indices[i, 1:] = neighbors[1][sorted_indices]
            distances[i, 1:] = D[i][
                neighbors[0][sorted_indices], neighbors[1][sorted_indices]]
        else:
            indices[i, 1:] = neighbors[1]
            distances[i, 1:] = D[i][neighbors]
    return indices, distances


def get_indices_distances_from_dense_matrix(D, n_neighbors):
    sample_range = np.arange(D.shape[0])[:, None]
    indices = np.argpartition(D, n_neighbors-1, axis=1)[:, :n_neighbors]
    indices = indices[sample_range, np.argsort(D[sample_range, indices])]
    distances = D[sample_range, indices]
    return indices, distances


def _backwards_compat_get_full_X_diffmap(adata):
    if 'X_diffmap0' in adata.obs:
        return np.c_[adata.obs['X_diffmap0'].values[:, None],
                     adata.obsm['X_diffmap']]
    else:
        return adata.obsm['X_diffmap']


def _backwards_compat_get_full_eval(adata):
    if 'X_diffmap0' in adata.obs:
        return np.r_[1, adata.uns['diffmap_evals']]
    else:
        return adata.uns['diffmap_evals']


class OnFlySymMatrix():
    """Emulate a matrix where elements are calculated on the fly.
    """
    def __init__(self, get_row, shape, DC_start=0, DC_end=-1, rows=None, restrict_array=None):
        self.get_row = get_row
        self.shape = shape
        self.DC_start = DC_start
        self.DC_end = DC_end
        self.rows = {} if rows is None else rows
        self.restrict_array = restrict_array  # restrict the array to a subset

    def __getitem__(self, index):
        if isinstance(index, int) or isinstance(index, np.integer):
            if self.restrict_array is None:
                glob_index = index
            else:
                # map the index back to the global index
                glob_index = self.restrict_array[index]
            if glob_index not in self.rows:
                self.rows[glob_index] = self.get_row(glob_index)
            row = self.rows[glob_index]
            if self.restrict_array is None:
                return row
            else:
                return row[self.restrict_array]
        else:
            if self.restrict_array is None:
                glob_index_0, glob_index_1 = index
            else:
                glob_index_0 = self.restrict_array[index[0]]
                glob_index_1 = self.restrict_array[index[1]]
            if glob_index_0 not in self.rows:
                self.rows[glob_index_0] = self.get_row(glob_index_0)
            return self.rows[glob_index_0][glob_index_1]

    def restrict(self, index_array):
        """Generate a view restricted to a subset of indices.
        """
        new_shape = index_array.shape[0], index_array.shape[0]
        return OnFlySymMatrix(self.get_row, new_shape, DC_start=self.DC_start,
                              DC_end=self.DC_end,
                              rows=self.rows, restrict_array=index_array)


class Neighbors():
    """Data represented as graph of nearest neighbors.

    Represent a data matrix as a graph of nearest neighbor relations (edges)
    among data points (nodes).

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        An annotated data matrix.

    Attributes
    ----------
    distances
    connectivities
    transitions
    transitions_sym
    eigen_values
    eigen_basis
    laplacian
    distances_dpt

    Methods
    -------
    compute_neighbors
    compute_transitions
    compute_eigen
    to_igraph
    """

    def __init__(self, adata):
        self._adata = adata
        self._init_iroot()
        # use the graph in adata
        info_str = ''
        self.knn = None
        self._distances = None
        self._connectivities = None
        self._number_connected_components = None
        if 'neighbors' in adata.uns:
            if 'distances' in adata.uns['neighbors']:
                self.knn = issparse(adata.uns['neighbors']['distances'])
                self._distances = adata.uns['neighbors']['distances']
            if 'connectivities' in adata.uns['neighbors']:
                self.knn = issparse(adata.uns['neighbors']['connectivities'])
                self._connectivities = adata.uns['neighbors']['connectivities']
            if 'params' in adata.uns['neighbors']:
                self.n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
            else:
                # estimating n_neighbors
                if self._connectivities is None:
                    self.n_neighbors = int(self._distances.count_nonzero() / self._distances.shape[0])
                else:
                    self.n_neighbors = int(self._connectivities.count_nonzero() / self._connectivities.shape[0] / 2)
            info_str += '`.distances` `.connectivities` '
            self._number_connected_components = 1
            if issparse(self._connectivities):
                from scipy.sparse.csgraph import connected_components
                self._connected_components = connected_components(self._connectivities)
        if 'X_diffmap' in adata.obsm_keys():
            self._eigen_values = _backwards_compat_get_full_eval(adata)
            self._eigen_basis = _backwards_compat_get_full_X_diffmap(adata)
            self.n_dcs = len(self._eigen_values)
            info_str += '`.eigen_values` `.eigen_basis` `.distances_dpt`'
        else:
            self._eigen_values = None
            self._eigen_basis = None
            self.n_dcs = None
        if info_str != '':
            logg.msg('    initialized {}'.format(info_str), v=4)

    @property
    def distances(self):
        """Distances between data points (sparse matrix).
        """
        return self._distances

    @property
    def connectivities(self):
        """Connectivities between data points (sparse matrix).
        """
        return self._connectivities

    @property
    def transitions(self):
        """Transition matrix (sparse matrix).

        Note: this has not been tested, in contrast to `transitions_sym`.

        Is conjugate to the symmetrized transition matrix via::

            self.transitions = self.Z *  self.transitions_sym / self.Z

        where ``self.Z`` is the diagonal matrix storing the normalization of the
        underlying kernel matrix.
        """
        if issparse(self.Z):
            Zinv = self.Z.power(-1)
        else:
            Zinv = np.diag(1./np.diag(self.Z))
        return self.Z.dot(self.transitions_sym).dot(Zinv)

    @property
    def transitions_sym(self):
        """Symmetrized transition matrix (sparse matrix).

        Is conjugate to the transition matrix via::

            self.transitions_sym = self.Z /  self.transitions * self.Z

        where ``self.Z`` is the diagonal matrix storing the normalization of the
        underlying kernel matrix.
        """
        return self._transitions_sym

    @property
    def eigen_values(self):
        """Eigen values of transition matrix (numpy array).
        """
        return self._eigen_values

    @property
    def eigen_basis(self):
        """Eigen basis of transition matrix (numpy array).
        """
        return self._eigen_basis

    @property
    def laplacian(self):
        """Graph laplacian (sparse matrix).
        """
        return self._laplacian

    @property
    def distances_dpt(self):
        """DPT distances (on-fly matrix).

        This is yields [Haghverdi16]_, Eq. 15 from the supplement with the
        extensions of [Wolf17i]_, supplement on random-walk based distance
        measures.
        """
        return OnFlySymMatrix(self._get_dpt_row, shape=self._adata.shape)

    def to_igraph(self):
        """Generate igraph from connectiviies.
        """
        return utils.get_igraph_from_adjacency(self.connectivities)

    @doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
    def compute_neighbors(
            self,
            n_neighbors=30,
            knn=True,
            n_pcs=None,
            use_rep=None,
            method='umap',
            random_state=0,
            precompute_metric=None,
            metric='euclidean',
            metric_kwds={}):
        """\
        Compute distances and connectivities of neighbors.

        Parameters
        ----------
        n_neighbors : `int`, optional (default: 30)
             Use this number of nearest neighbors.
        knn : `bool`, optional (default: `True`)
             Restrict result to `n_neighbors` nearest neighbors.
        {n_pcs}
        {use_rep}

        Returns
        -------
        Writes attributes `.distances` and `.connectivities`.
        """
        if n_neighbors > self._adata.shape[0]:  # very small datasets
            n_neighbors = 1 + int(0.5*self._adata.shape[0])
            logg.warn('n_obs too small: adjusting to `n_neighbors = {}`'
                      .format(n_neighbors))
        if method == 'umap' and not knn:
            raise ValueError('`method = \'umap\' only with `knn = True`.')
        if method not in {'umap', 'gauss'}:
            raise ValueError('`method` needs to be \'umap\' or \'gauss\'.')
        if self._adata.shape[0] >= 10000 and not knn:
            logg.warn(
                'Using high n_obs without `knn=True` takes a lot of memory...')
        self.n_neighbors = n_neighbors
        self.knn = knn
        X = choose_representation(self._adata, use_rep=use_rep, n_pcs=n_pcs)
        # neighbor search
        use_dense_distances = (metric == 'euclidean' and X.shape[0] < 8192) or knn == False
        if use_dense_distances:
            # standard eulcidean case for relatively small matrices
            self._distances, knn_indices, knn_distances = compute_neighbors_numpy(
                X, n_neighbors, knn=knn)
        else:
            # non-euclidean case and approx nearest neighbors
            if X.shape[0] < 4096:
                X = pairwise_distances(X, metric=metric, **metric_kwds)
                metric = 'precomputed'
            knn_indices, knn_distances = compute_neighbors_umap(
                X, n_neighbors, random_state, metric=metric, metric_kwds=metric_kwds)
        logg.msg('computed neighbors', t=True, v=4)
        if not use_dense_distances or method == 'umap':
            # we need self._distances also for method == 'gauss' if we didn't
            # use dense distances
            self._distances, self._connectivities = compute_connectivities_umap(
                knn_indices, knn_distances, self._adata.shape[0], self.n_neighbors)
        # overwrite the umap connectivities if method is 'gauss'
        # self._distances is unaffected by this
        if method == 'gauss':
            self._compute_connectivities_diffmap()
        logg.msg('computed connectivities', t=True, v=4)
        self._number_connected_components = 1
        if issparse(self._connectivities):
            from scipy.sparse.csgraph import connected_components
            self._connected_components = connected_components(self._connectivities)

    def _compute_connectivities_diffmap(self, density_normalize=True):
        # init distances
        if self.knn:
            Dsq = self._distances.power(2)
            indices, distances_sq = get_indices_distances_from_sparse_matrix(
                Dsq, self.n_neighbors)
        else:
            Dsq = np.power(self._distances, 2)
            indices, distances_sq = get_indices_distances_from_dense_matrix(
                Dsq, self.n_neighbors)

        # exclude the first point, the 0th neighbor
        indices = indices[:, 1:]
        distances_sq = distances_sq[:, 1:]

        # choose sigma, the heuristic here doesn't seem to make much of a difference,
        # but is used to reproduce the figures of Haghverdi et al. (2016)
        if self.knn:
            # as the distances are not sorted
            # we have decay within the n_neighbors first neighbors
            sigmas_sq = np.median(distances_sq, axis=1)
        else:
            # the last item is already in its sorted position through argpartition
            # we have decay beyond the n_neighbors neighbors
            sigmas_sq = distances_sq[:, -1]/4
        sigmas = np.sqrt(sigmas_sq)

        # compute the symmetric weight matrix
        if not issparse(self._distances):
            Num = 2 * np.multiply.outer(sigmas, sigmas)
            Den = np.add.outer(sigmas_sq, sigmas_sq)
            W = np.sqrt(Num/Den) * np.exp(-Dsq/Den)
            # make the weight matrix sparse
            if not self.knn:
                mask = W > 1e-14
                W[mask == False] = 0
            else:
                # restrict number of neighbors to ~k
                # build a symmetric mask
                mask = np.zeros(Dsq.shape, dtype=bool)
                for i, row in enumerate(indices):
                    mask[i, row] = True
                    for j in row:
                        if i not in set(indices[j]):
                            W[j, i] = W[i, j]
                            mask[j, i] = True
                # set all entries that are not nearest neighbors to zero
                W[mask == False] = 0
        else:
            W = Dsq.copy()  # need to copy the distance matrix here; what follows is inplace
            for i in range(len(Dsq.indptr[:-1])):
                row = Dsq.indices[Dsq.indptr[i]:Dsq.indptr[i+1]]
                num = 2 * sigmas[i] * sigmas[row]
                den = sigmas_sq[i] + sigmas_sq[row]
                W.data[Dsq.indptr[i]:Dsq.indptr[i+1]] = np.sqrt(num/den) * np.exp(
                    -Dsq.data[Dsq.indptr[i]: Dsq.indptr[i+1]] / den)
            W = W.tolil()
            for i, row in enumerate(indices):
                for j in row:
                    if i not in set(indices[j]):
                        W[j, i] = W[i, j]
            W = W.tocsr()

        self._connectivities = W

    def compute_transitions(self, density_normalize=True):
        """Compute transition matrix.

        Parameters
        ----------
        density_normalize : `bool`
            The density rescaling of Coifman and Lafon (2006): Then only the
            geometry of the data matters, not the sampled density.

        Returns
        -------
        Makes attributes `.transitions_sym` and `.transitions` available.
        """
        W = self._connectivities
        # density normalization as of Coifman et al. (2005)
        # ensures that kernel matrix is independent of sampling density
        if density_normalize:
            # q[i] is an estimate for the sampling density at point i
            # it's also the degree of the underlying graph
            q = np.asarray(W.sum(axis=0))
            if not issparse(W):
                Q = np.diag(1.0/q)
            else:
                Q = scipy.sparse.spdiags(1.0/q, 0, W.shape[0], W.shape[0])
            K = Q.dot(W).dot(Q)
        else:
            K = W

        # z[i] is the square root of the row sum of K
        z = np.sqrt(np.asarray(K.sum(axis=0)))
        if not issparse(K):
            self.Z = np.diag(1.0/z)
        else:
            self.Z = scipy.sparse.spdiags(1.0/z, 0, K.shape[0], K.shape[0])
        self._transitions_sym = self.Z.dot(K).dot(self.Z)
        logg.msg('computed transitions', v=4, time=True)

    def compute_eigen(self, n_comps=15, sym=None, sort='decrease'):
        """Compute eigen decomposition of transition matrix.

        Parameters
        ----------
        n_comps : `int`
            Number of eigenvalues/vectors to be computed, set `n_comps = 0` if
            you need all eigenvectors.
        sym : `bool`
            Instead of computing the eigendecomposition of the assymetric
            transition matrix, computed the eigendecomposition of the symmetric
            Ktilde matrix.
        matrix : sparse matrix, np.ndarray, optional (default: `.connectivities`)
            Matrix to diagonalize. Merely for testing and comparison purposes.

        Returns
        -------
        Writes the following attributes.
        eigen_values : np.ndarray
            Eigenvalues of transition matrix.
        eigen_basis : np.ndarray
             Matrix of eigenvectors (stored in columns).  `.eigen_basis` is
             projection of data matrix on right eigenvectors, that is, the
             projection on the diffusion components.  these are simply the
             components of the right eigenvectors and can directly be used for
             plotting.
        """
        np.set_printoptions(precision=10)
        if self._transitions_sym is None:
            raise ValueError('Run `.compute_transitions` first.')
        matrix = self._transitions_sym
        # compute the spectrum
        if n_comps == 0:
            evals, evecs = scipy.linalg.eigh(matrix)
        else:
            n_comps = min(matrix.shape[0]-1, n_comps)
            # ncv = max(2 * n_comps + 1, int(np.sqrt(matrix.shape[0])))
            ncv = None
            which = 'LM' if sort == 'decrease' else 'SM'
            # it pays off to increase the stability with a bit more precision
            matrix = matrix.astype(np.float64)
            evals, evecs = scipy.sparse.linalg.eigsh(matrix, k=n_comps,
                                                  which=which, ncv=ncv)
            evals, evecs = evals.astype(np.float32), evecs.astype(np.float32)
        if sort == 'decrease':
            evals = evals[::-1]
            evecs = evecs[:, ::-1]
        logg.info('    eigenvalues of transition matrix\n'
                  '    {}'.format(str(evals).replace('\n', '\n    ')))
        if self._connected_components[0] > len(evals)/2:
            logg.warn('Transition matrix has many disconnected components!')
        self._eigen_values = evals
        self._eigen_basis = evecs

    def _init_iroot(self):
        self.iroot = None
        # set iroot directly
        if 'iroot' in self._adata.uns:
            if self._adata.uns['iroot'] >= self._adata.n_obs:
                logg.warn('Root cell index {} does not exist for {} samples. '
                          'Is ignored.'
                          .format(self._adata.uns['iroot'], self._adata.n_obs))
            else:
                self.iroot = self._adata.uns['iroot']
            return
        # set iroot via xroot
        xroot = None
        if 'xroot' in self._adata.uns: xroot = self._adata.uns['xroot']
        elif 'xroot' in self._adata.var: xroot = self._adata.var['xroot']
        # see whether we can set self.iroot using the full data matrix
        if xroot is not None and xroot.size == self._adata.shape[1]:
            self._set_iroot_via_xroot(xroot)

    def _get_dpt_row(self, i):
        use_mask = False
        if self._connected_components[0] > 1:
            use_mask = True
            label = self._connected_components[1][i]
            mask = self._connected_components[1] == label
        row = sum([(self.eigen_values[l]/(1-self.eigen_values[l])
                     * (self.eigen_basis[i, l] - self.eigen_basis[:, l]))**2
                   # account for float32 precision
                    for l in range(0, self.eigen_values.size) if self.eigen_values[l] < 0.9994])
        row += sum([(self.eigen_basis[i, l] - self.eigen_basis[:, l])**2
                    for l in range(0, self.eigen_values.size) if self.eigen_values[l] >= 0.9994])
        if not use_mask:
            return np.sqrt(row)
        else:
            row[~mask] = np.inf
            return np.sqrt(row)

    def _compute_Lp_matrix(self):
        """See Fouss et al. (2006) and von Luxburg et al. (2007).

        See Proposition 6 in von Luxburg (2007) and the inline equations
        right in the text above.
        """
        self.Lp = sum([1/self.eigen_values[i]
                      * np.outer(self.eigen_basis[:, i], self.eigen_basis[:, i])
                      for i in range(1, self.eigen_values.size)])

    def _compute_C_matrix(self):
        """See Fouss et al. (2006) and von Luxburg et al. (2007).

        This is the commute-time matrix. It's a squared-euclidian distance
        matrix in \mathbb{R}^n.
        """
        self.C = np.repeat(np.diag(self.Lp)[:, np.newaxis],
                           self.Lp.shape[0], axis=1)
        self.C += np.repeat(np.diag(self.Lp)[np.newaxis, :],
                            self.Lp.shape[0], axis=0)
        self.C -= 2*self.Lp
        # the following is much slower
        # self.C = np.zeros(self.Lp.shape)
        # for i in range(self.Lp.shape[0]):
        #     for j in range(self.Lp.shape[1]):
        #         self.C[i, j] = self.Lp[i, i] + self.Lp[j, j] - 2*self.Lp[i, j]
        volG = np.sum(self.z)
        self.C *= volG
        settings.mt(0, 'computed commute distance matrix')
        self.distances_dpt = self.C

    def _compute_MFP_matrix(self):
        """See Fouss et al. (2006).

        This is the mean-first passage time matrix. It's not a distance.

        Mfp[i, k] := m(k|i) in the notation of Fouss et al. (2006). This
        corresponds to the standard notation for transition matrices (left index
        initial state, right index final state, i.e. a right-stochastic
        matrix, with each row summing to one).
        """
        self.MFP = np.zeros(self.Lp.shape)
        for i in range(self.Lp.shape[0]):
            for k in range(self.Lp.shape[1]):
                for j in range(self.Lp.shape[1]):
                    self.MFP[i, k] += (self.Lp[i, j] - self.Lp[i, k]
                                       - self.Lp[k, j] + self.Lp[k, k]) * self.z[j]
        settings.mt(0, 'computed mean first passage time matrix')
        self.distances_dpt = self.MFP

    def _set_pseudotime(self):
        """Return pseudotime with respect to root point.
        """
        self.pseudotime = self.distances_dpt[self.iroot].copy()
        self.pseudotime /= np.max(self.pseudotime[self.pseudotime < np.inf])

    def _set_iroot_via_xroot(self, xroot):
        """Determine the index of the root cell.

        Given an expression vector, find the observation index that is closest
        to this vector.

        Parameters
        ----------
        xroot : np.ndarray
            Vector that marks the root cell, the vector storing the initial
            condition, only relevant for computing pseudotime.
        """
        if self._adata.shape[1] != xroot.size:
            raise ValueError(
                'The root vector you provided does not have the '
                'correct dimension.')
        # this is the squared distance
        dsqroot = 1e10
        iroot = 0
        for i in range(self._adata.shape[0]):
            diff = self._adata.X[i, :] - xroot
            dsq = diff.dot(diff)
            if dsq < dsqroot:
                dsqroot = dsq
                iroot = i
                if np.sqrt(dsqroot) < 1e-10: break
        logg.msg('setting root index to', iroot, v=4)
        if self.iroot is not None and iroot != self.iroot:
            logg.warn('Changing index of iroot from {} to {}.'.format(self.iroot, iroot))
        self.iroot = iroot

    def _test_embed(self):
        """
        Checks and tests for embed.
        """
        # pl.semilogy(w,'x',label=r'$ \widetilde K$')
        # pl.show()
        if settings.verbosity > 2:
            # output of spectrum of K for comparison
            w, v = np.linalg.eigh(self.K)
            logg.msg('spectrum of K (kernel)')
        if settings.verbosity > 3:
            # direct computation of spectrum of T
            w, vl, vr = scipy.linalg.eig(self.T, left=True)
            logg.msg('spectrum of transition matrix (should be same as of Ktilde)')
