from __future__ import annotations

import sys
from types import MappingProxyType
from typing import (
    TYPE_CHECKING,
    TypedDict,
    Union,
    Optional,
    Any,
    NamedTuple,
    Literal,
    get_args,
)
from collections.abc import Mapping, MutableMapping, Callable
from warnings import warn

import numpy as np
import scipy
from anndata import AnnData
from scipy.sparse import issparse, csr_matrix
from sklearn.utils import check_random_state
from pynndescent import NNDescent, PyNNDescentTransformer

if TYPE_CHECKING:
    from igraph import Graph
    from ._types import KnnTransformerLike

from . import _connectivity
from ._types import _Metric, _MetricFn, _Method, _KnownTransformer
from ._common import (
    _get_indices_distances_from_sparse_matrix,
    _get_sparse_matrix_from_indices_distances,
)
from .. import logging as logg
from .. import _utils, settings
from .._utils import _doc_params, AnyRandom, NeighborsView
from ..tools._utils import _choose_representation, doc_use_rep, doc_n_pcs

if sys.version_info >= (3, 9):
    RPForestDict = Mapping[str, Mapping[str, np.ndarray]]
else:
    RPForestDict = Mapping


N_DCS = 15  # default number of diffusion components
# Backwards compat, constants should be defined in only one place.
N_PCS = settings.N_PCS


class KwdsForTransformer(TypedDict):
    """Keyword arguments passed to a _KnownTransformer.

    IMPORTANT: when changing the parameters set here,
    update the “*ignored*” part in the parameter docs!
    """

    n_neighbors: int
    metric: Union[_Metric, _MetricFn]
    metric_params: Mapping[str, Any]
    random_state: AnyRandom


@_doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: Optional[int] = None,
    *,
    use_rep: Optional[str] = None,
    knn: bool = True,
    method: _Method = 'umap',
    transformer: KnnTransformerLike | _KnownTransformer | None = None,
    metric: Union[_Metric, _MetricFn] = 'euclidean',
    metric_kwds: Mapping[str, Any] = MappingProxyType({}),
    random_state: AnyRandom = 0,
    key_added: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Compute a neighborhood graph of observations [McInnes18]_.

    The neighbor search efficiency of this heavily relies on UMAP [McInnes18]_,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='gauss'`,
    connectivities are computed according to [Coifman05]_, in the adaption of
    [Haghverdi16]_.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.

        *ignored if ``transformer`` is an instance.*
    {n_pcs}
    {use_rep}
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    method
        Use 'umap' [McInnes18]_ or 'gauss' (Gauss kernel following [Coifman05]_
        with adaptive width [Haghverdi16]_) for computing connectivities.
    transformer
        Approximate kNN search implementation following the API of
        :class:`~sklearn.neighbors.KNeighborsTransformer`.
        Also accepts the following known options:

        `None` (the default)
            Behavior depends on data size.
            For small data, we will calculate exact kNN, otherwise we use
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'pynndescent'`
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'rapids'`
            A transformer based on :class:`cuml.neighbors.NearestNeighbors`.
    metric
        A known metric’s name or a callable that returns a distance.

        *ignored if ``transformer`` is an instance.*
    metric_kwds
        Options for the metric.

        *ignored if ``transformer`` is an instance.*
    random_state
        A numpy random seed.

        *ignored if ``transformer`` is an instance.*
    key_added
        If not specified, the neighbors data is stored in .uns['neighbors'],
        distances and connectivities are stored in .obsp['distances'] and
        .obsp['connectivities'] respectively.
        If specified, the neighbors data is added to .uns[key_added],
        distances are stored in .obsp[key_added+'_distances'] and
        connectivities in .obsp[key_added+'_connectivities'].
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:

    See `key_added` parameter description for the storage path of
    connectivities and distances.

    **connectivities** : sparse matrix of dtype `float32`.
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    **distances** : sparse matrix of dtype `float64`.
        Instead of decaying weights, this stores distances for each pair of
        neighbors.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> # Basic usage
    >>> sc.pp.neighbors(adata, 20, metric='cosine')
    >>> # Provide your own transformer for more control and flexibility
    >>> from sklearn.neighbors import KNeighborsTransformer
    >>> transformer = KNeighborsTransformer(n_neighbors=10, metric='cosine', algorithm='ball_tree')
    >>> sc.pp.neighbors(adata, transformer=transformer)
    >>> # now you can e.g. access the index: `transformer._tree`
    """
    start = logg.info('computing neighbors')
    adata = adata.copy() if copy else adata
    if adata.is_view:  # we shouldn't need this here...
        adata._init_as_actual(adata.copy())
    neighbors = Neighbors(adata)
    neighbors.compute_neighbors(
        n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        knn=knn,
        method=method,
        transformer=transformer,
        metric=metric,
        metric_kwds=metric_kwds,
        random_state=random_state,
    )

    if key_added is None:
        key_added = 'neighbors'
        conns_key = 'connectivities'
        dists_key = 'distances'
    else:
        conns_key = key_added + '_connectivities'
        dists_key = key_added + '_distances'

    adata.uns[key_added] = {}

    neighbors_dict = adata.uns[key_added]

    neighbors_dict['connectivities_key'] = conns_key
    neighbors_dict['distances_key'] = dists_key

    neighbors_dict['params'] = dict(
        n_neighbors=neighbors.n_neighbors,
        method=method,
        random_state=random_state,
        metric=metric,
    )
    if metric_kwds:
        neighbors_dict['params']['metric_kwds'] = metric_kwds
    if use_rep is not None:
        neighbors_dict['params']['use_rep'] = use_rep
    if n_pcs is not None:
        neighbors_dict['params']['n_pcs'] = n_pcs

    adata.obsp[dists_key] = neighbors.distances
    adata.obsp[conns_key] = neighbors.connectivities

    if neighbors.rp_forest is not None:
        neighbors_dict['rp_forest'] = neighbors.rp_forest
    logg.info(
        '    finished',
        time=start,
        deep=(
            f'added to `.uns[{key_added!r}]`\n'
            f'    `.obsp[{dists_key!r}]`, distances for each pair of neighbors\n'
            f'    `.obsp[{conns_key!r}]`, weighted adjacency matrix'
        ),
    )
    return adata if copy else None


class FlatTree(NamedTuple):
    hyperplanes: None
    offsets: None
    children: None
    indices: None


def _backwards_compat_get_full_X_diffmap(adata: AnnData) -> np.ndarray:
    if 'X_diffmap0' in adata.obs:
        return np.c_[adata.obs['X_diffmap0'].values[:, None], adata.obsm['X_diffmap']]
    else:
        return adata.obsm['X_diffmap']


def _backwards_compat_get_full_eval(adata: AnnData):
    if 'X_diffmap0' in adata.obs:
        return np.r_[1, adata.uns['diffmap_evals']]
    else:
        return adata.uns['diffmap_evals']


def _make_forest_dict(forest):
    d = {}
    props = ('hyperplanes', 'offsets', 'children', 'indices')
    for prop in props:
        d[prop] = {}
        sizes = np.fromiter(
            (getattr(tree, prop).shape[0] for tree in forest), dtype=int
        )
        d[prop]['start'] = np.zeros_like(sizes)
        if prop == 'offsets':
            dims = sizes.sum()
        else:
            dims = (sizes.sum(), getattr(forest[0], prop).shape[1])
        dtype = getattr(forest[0], prop).dtype
        dat = np.empty(dims, dtype=dtype)
        start = 0
        for i, size in enumerate(sizes):
            d[prop]['start'][i] = start
            end = start + size
            dat[start:end] = getattr(forest[i], prop)
            start = end
        d[prop]['data'] = dat
    return d


class OnFlySymMatrix:
    """Emulate a matrix where elements are calculated on the fly."""

    def __init__(
        self,
        get_row: Callable[[Any], np.ndarray],
        shape: tuple[int, int],
        DC_start: int = 0,
        DC_end: int = -1,
        rows: Optional[MutableMapping[Any, np.ndarray]] = None,
        restrict_array: Optional[np.ndarray] = None,
    ):
        self.get_row = get_row
        self.shape = shape
        self.DC_start = DC_start
        self.DC_end = DC_end
        self.rows = {} if rows is None else rows
        self.restrict_array = restrict_array  # restrict the array to a subset

    def __getitem__(self, index):
        if isinstance(index, (int, np.integer)):
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
        """Generate a view restricted to a subset of indices."""
        new_shape = index_array.shape[0], index_array.shape[0]
        return OnFlySymMatrix(
            self.get_row,
            new_shape,
            DC_start=self.DC_start,
            DC_end=self.DC_end,
            rows=self.rows,
            restrict_array=index_array,
        )


class Neighbors:
    """\
    Data represented as graph of nearest neighbors.

    Represent a data matrix as a graph of nearest neighbor relations (edges)
    among data points (nodes).

    Parameters
    ----------
    adata
        Annotated data object.
    n_dcs
        Number of diffusion components to use.
    neighbors_key
        Where to look in .uns and .obsp for neighbors data
    """

    def __init__(
        self,
        adata: AnnData,
        n_dcs: Optional[int] = None,
        neighbors_key: Optional[str] = None,
    ):
        self._adata = adata
        self._init_iroot()
        # use the graph in adata
        info_str = ''
        self.knn: Optional[bool] = None
        self._distances: Union[np.ndarray, csr_matrix, None] = None
        self._connectivities: Union[np.ndarray, csr_matrix, None] = None
        self._transitions_sym: Union[np.ndarray, csr_matrix, None] = None
        self._number_connected_components: Optional[int] = None
        self._rp_forest: Optional[RPForestDict] = None
        if neighbors_key is None:
            neighbors_key = 'neighbors'
        if neighbors_key in adata.uns:
            neighbors = NeighborsView(adata, neighbors_key)
            if 'distances' in neighbors:
                self.knn = issparse(neighbors['distances'])
                self._distances = neighbors['distances']
            if 'connectivities' in neighbors:
                self.knn = issparse(neighbors['connectivities'])
                self._connectivities = neighbors['connectivities']
            if 'rp_forest' in neighbors:
                self._rp_forest = neighbors['rp_forest']
            if 'params' in neighbors:
                self.n_neighbors = neighbors['params']['n_neighbors']
            else:

                def count_nonzero(a: Union[np.ndarray, csr_matrix]) -> int:
                    return a.count_nonzero() if issparse(a) else np.count_nonzero(a)

                # estimating n_neighbors
                if self._connectivities is None:
                    self.n_neighbors = int(
                        count_nonzero(self._distances) / self._distances.shape[0]
                    )
                else:
                    self.n_neighbors = int(
                        count_nonzero(self._connectivities)
                        / self._connectivities.shape[0]
                        / 2
                    )
            info_str += '`.distances` `.connectivities` '
            self._number_connected_components = 1
            if issparse(self._connectivities):
                from scipy.sparse.csgraph import connected_components

                self._connected_components = connected_components(self._connectivities)
                self._number_connected_components = self._connected_components[0]
        if 'X_diffmap' in adata.obsm_keys():
            self._eigen_values = _backwards_compat_get_full_eval(adata)
            self._eigen_basis = _backwards_compat_get_full_X_diffmap(adata)
            if n_dcs is not None:
                if n_dcs > len(self._eigen_values):
                    raise ValueError(
                        'Cannot instantiate using `n_dcs`={}. '
                        'Compute diffmap/spectrum with more components first.'.format(
                            n_dcs
                        )
                    )
                self._eigen_values = self._eigen_values[:n_dcs]
                self._eigen_basis = self._eigen_basis[:, :n_dcs]
            self.n_dcs = len(self._eigen_values)
            info_str += '`.eigen_values` `.eigen_basis` `.distances_dpt`'
        else:
            self._eigen_values = None
            self._eigen_basis = None
            self.n_dcs = None
        if info_str != '':
            logg.debug(f'    initialized {info_str}')

    @property
    def rp_forest(self) -> Optional[RPForestDict]:
        return self._rp_forest

    @property
    def distances(self) -> Union[np.ndarray, csr_matrix, None]:
        """Distances between data points (sparse matrix)."""
        return self._distances

    @property
    def connectivities(self) -> Union[np.ndarray, csr_matrix, None]:
        """Connectivities between data points (sparse matrix)."""
        return self._connectivities

    @property
    def transitions(self) -> Union[np.ndarray, csr_matrix]:
        """Transition matrix (sparse matrix).

        Is conjugate to the symmetrized transition matrix via::

            self.transitions = self.Z *  self.transitions_sym / self.Z

        where ``self.Z`` is the diagonal matrix storing the normalization of the
        underlying kernel matrix.

        Notes
        -----
        This has not been tested, in contrast to `transitions_sym`.
        """
        if issparse(self.Z):
            Zinv = self.Z.power(-1)
        else:
            Zinv = np.diag(1.0 / np.diag(self.Z))
        return self.Z @ self.transitions_sym @ Zinv

    @property
    def transitions_sym(self) -> Union[np.ndarray, csr_matrix, None]:
        """Symmetrized transition matrix (sparse matrix).

        Is conjugate to the transition matrix via::

            self.transitions_sym = self.Z /  self.transitions * self.Z

        where ``self.Z`` is the diagonal matrix storing the normalization of the
        underlying kernel matrix.
        """
        return self._transitions_sym

    @property
    def eigen_values(self) -> np.ndarray:
        """Eigen values of transition matrix."""
        return self._eigen_values

    @property
    def eigen_basis(self) -> np.ndarray:
        """Eigen basis of transition matrix."""
        return self._eigen_basis

    @property
    def distances_dpt(self) -> OnFlySymMatrix:
        """DPT distances.

        This is yields [Haghverdi16]_, Eq. 15 from the supplement with the
        extensions of [Wolf19]_, supplement on random-walk based distance
        measures.
        """
        return OnFlySymMatrix(self._get_dpt_row, shape=self._adata.shape)

    def to_igraph(self) -> Graph:
        """Generate igraph from connectiviies."""
        return _utils.get_igraph_from_adjacency(self.connectivities)

    @_doc_params(n_pcs=doc_n_pcs, use_rep=doc_use_rep)
    def compute_neighbors(
        self,
        n_neighbors: int = 30,
        n_pcs: Optional[int] = None,
        *,
        use_rep: Optional[str] = None,
        knn: bool = True,
        method: _Method = 'umap',
        transformer: KnnTransformerLike | _KnownTransformer | None = None,
        metric: Union[_Metric, _MetricFn] = 'euclidean',
        metric_kwds: Mapping[str, Any] = MappingProxyType({}),
        random_state: AnyRandom = 0,
    ) -> None:
        """\
        Compute distances and connectivities of neighbors.

        Parameters
        ----------
        n_neighbors
            Use this number of nearest neighbors.
        {n_pcs}
        {use_rep}
        knn
            Restrict result to `n_neighbors` nearest neighbors.

        Returns
        -------
        Writes sparse graph attributes `.distances` and `.connectivities`.
        """
        start_neighbors = logg.debug('computing neighbors')
        if n_neighbors > self._adata.shape[0]:  # very small datasets
            n_neighbors = 1 + int(0.5 * self._adata.shape[0])
            logg.warning(f'n_obs too small: adjusting to `n_neighbors = {n_neighbors}`')

        # default keyword arguments when `transformer` is not an instance
        transformer_kwds_default = KwdsForTransformer(
            n_neighbors=n_neighbors,
            metric=metric,
            metric_params=metric_kwds,  # most use _params, not _kwds
            random_state=random_state,
        )
        method, transformer, shortcut = self._handle_transformer(
            method, transformer, knn=knn, kwds=transformer_kwds_default
        )

        if self._adata.shape[0] >= 10000 and not knn:
            logg.warning('Using high n_obs without `knn=True` takes a lot of memory...')
        # do not use the cached rp_forest
        self._rp_forest = None
        self.n_neighbors = n_neighbors
        self.knn = knn
        X = _choose_representation(self._adata, use_rep=use_rep, n_pcs=n_pcs)
        self._distances = transformer.fit_transform(X)
        knn_indices, knn_distances = _get_indices_distances_from_sparse_matrix(
            self._distances, n_neighbors
        )
        if shortcut:
            # self._distances is a sparse matrix with a diag of 1, fix that
            self._distances[np.diag_indices_from(self.distances)] = 0
            if knn:  # remove too far away entries and keep as sparse
                self._distances = _get_sparse_matrix_from_indices_distances(
                    knn_indices, knn_distances, self._adata.n_obs, n_neighbors
                )
            else:  # convert to dense
                self._distances = self._distances.toarray()
        if (index := getattr(transformer, "index_", None)) and isinstance(
            index, NNDescent
        ):
            # very cautious here
            try:
                self._rp_forest = _make_forest_dict(index)
            except Exception:  # TODO catch the correct exception
                pass

        start_connect = logg.debug('computed neighbors', time=start_neighbors)
        if method == 'umap':
            self._connectivities = _connectivity.umap(
                knn_indices,
                knn_distances,
                n_obs=self._adata.shape[0],
                n_neighbors=self.n_neighbors,
            )
        elif method == 'gauss':
            self._connectivities = _connectivity.gauss(
                self._distances, self.n_neighbors, knn=self.knn
            )
        else:
            msg = f'{method!r} should have been coerced in _handle_transform_args'
            raise AssertionError(msg)
        logg.debug('computed connectivities', time=start_connect)
        self._number_connected_components = 1
        if issparse(self._connectivities):
            from scipy.sparse.csgraph import connected_components

            self._connected_components = connected_components(self._connectivities)
            self._number_connected_components = self._connected_components[0]

    def _handle_transformer(
        self,
        method: _Method | Literal['gauss'],
        transformer: KnnTransformerLike | _KnownTransformer | None,
        *,
        knn: bool,
        kwds: KwdsForTransformer,
    ) -> tuple[_Method, KnnTransformerLike, bool]:
        """Return effective `method` and transformer.

        `method` will be coerced to `'gauss'` or `'umap'`.
        `transformer` is coerced from a str or instance to an instance class.

        If `transformer` is `None` and there are few data points,
        `transformer` will be set to a brute force
        :class:`~sklearn.neighbors.KNeighborsTransformer`.

        If `transformer` is `None` and there are many data points,
        `transformer` will be set like `umap` does (i.e. to a
        ~`pynndescent.PyNNDescentTransformer` with custom `n_trees` and `n_iter`).
        """
        # legacy logic
        use_dense_distances = (
            kwds['metric'] == 'euclidean' and self._adata.n_obs < 8192
        ) or not knn
        shortcut = transformer is None and (
            use_dense_distances or self._adata.n_obs < 4096
        )

        # Coerce `method` to 'gauss' or 'umap'
        if method == 'rapids':
            if transformer is not None:
                msg = "Can’t specify both `method = 'rapids'` and `transformer`."
                raise ValueError(msg)
            msg = "method = 'rapids' is deprecated. Use transformer = 'rapids'."
            warn(msg, FutureWarning)
            method = 'umap'
            transformer = 'rapids'
        elif method not in (methods := set(get_args(_Method))):
            msg = f'`method` needs to be one of {methods}.'
            raise ValueError(msg)

        # Validate `knn`
        conn_method = 'gauss' if method == 'gauss' else 'umap'
        if not knn and not (conn_method == 'gauss' and transformer is None):
            # “knn=False” seems to be only intended for method “gauss”
            msg = f'`method = {method!r} only with `knn = True`.'
            raise ValueError(msg)

        # Coerce `transformer` to an instance
        if shortcut:
            from sklearn.neighbors import KNeighborsTransformer

            assert transformer is None
            transformer = KNeighborsTransformer(
                algorithm='brute',
                n_jobs=settings.n_jobs,
                n_neighbors=self._adata.n_obs - 1,  # ignore n_neighbors
                metric=kwds['metric'],
                metric_params=dict(kwds['metric_params']),  # needs dict
                # no random_state
            )
        elif transformer is None or transformer == 'pynndescent':
            kwds = kwds.copy()
            kwds['metric_kwds'] = kwds.pop('metric_params')
            if transformer is None:
                # Use defaults from UMAP’s `nearest_neighbors` function
                kwds.update(
                    n_jobs=settings.n_jobs,
                    n_trees=min(64, 5 + int(round((self._adata.n_obs) ** 0.5 / 20.0))),
                    n_iters=max(5, int(round(np.log2(self._adata.n_obs)))),
                )
            transformer = PyNNDescentTransformer(**kwds)
        elif transformer == 'rapids':
            from scanpy.neighbors._backends.rapids import RapidsKNNTransformer

            transformer = RapidsKNNTransformer(**kwds)
        elif isinstance(transformer, str):
            msg = (
                f'Unknown transformer: {transformer}. '
                f'Try passing a class or one of {set(get_args(_KnownTransformer))}'
            )
            raise ValueError(msg)
        # else `transformer` is probably an instance
        return conn_method, transformer, shortcut

    def compute_transitions(self, density_normalize: bool = True):
        """\
        Compute transition matrix.

        Parameters
        ----------
        density_normalize
            The density rescaling of Coifman and Lafon (2006): Then only the
            geometry of the data matters, not the sampled density.

        Returns
        -------
        Makes attributes `.transitions_sym` and `.transitions` available.
        """
        start = logg.info('computing transitions')
        W = self._connectivities
        # density normalization as of Coifman et al. (2005)
        # ensures that kernel matrix is independent of sampling density
        if density_normalize:
            # q[i] is an estimate for the sampling density at point i
            # it's also the degree of the underlying graph
            q = np.asarray(W.sum(axis=0))
            if not issparse(W):
                Q = np.diag(1.0 / q)
            else:
                Q = scipy.sparse.spdiags(1.0 / q, 0, W.shape[0], W.shape[0])
            K = Q @ W @ Q
        else:
            K = W

        # z[i] is the square root of the row sum of K
        z = np.sqrt(np.asarray(K.sum(axis=0)))
        if not issparse(K):
            self.Z = np.diag(1.0 / z)
        else:
            self.Z = scipy.sparse.spdiags(1.0 / z, 0, K.shape[0], K.shape[0])
        self._transitions_sym = self.Z @ K @ self.Z
        logg.info('    finished', time=start)

    def compute_eigen(
        self,
        n_comps: int = 15,
        sym: Optional[bool] = None,
        sort: Literal['decrease', 'increase'] = 'decrease',
        random_state: AnyRandom = 0,
    ):
        """\
        Compute eigen decomposition of transition matrix.

        Parameters
        ----------
        n_comps
            Number of eigenvalues/vectors to be computed, set `n_comps = 0` if
            you need all eigenvectors.
        sym
            Instead of computing the eigendecomposition of the assymetric
            transition matrix, computed the eigendecomposition of the symmetric
            Ktilde matrix.
        random_state
            A numpy random seed

        Returns
        -------
        Writes the following attributes.

        eigen_values : numpy.ndarray
            Eigenvalues of transition matrix.
        eigen_basis : numpy.ndarray
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
            n_comps = min(matrix.shape[0] - 1, n_comps)
            # ncv = max(2 * n_comps + 1, int(np.sqrt(matrix.shape[0])))
            ncv = None
            which = 'LM' if sort == 'decrease' else 'SM'
            # it pays off to increase the stability with a bit more precision
            matrix = matrix.astype(np.float64)

            # Setting the random initial vector
            random_state = check_random_state(random_state)
            v0 = random_state.standard_normal(matrix.shape[0])
            evals, evecs = scipy.sparse.linalg.eigsh(
                matrix, k=n_comps, which=which, ncv=ncv, v0=v0
            )
            evals, evecs = evals.astype(np.float32), evecs.astype(np.float32)
        if sort == 'decrease':
            evals = evals[::-1]
            evecs = evecs[:, ::-1]
        logg.info(
            '    eigenvalues of transition matrix\n'
            '    {}'.format(str(evals).replace('\n', '\n    '))
        )
        if self._number_connected_components > len(evals) / 2:
            logg.warning('Transition matrix has many disconnected components!')
        self._eigen_values = evals
        self._eigen_basis = evecs

    def _init_iroot(self):
        self.iroot = None
        # set iroot directly
        if 'iroot' in self._adata.uns:
            if self._adata.uns['iroot'] >= self._adata.n_obs:
                logg.warning(
                    f'Root cell index {self._adata.uns["iroot"]} does not '
                    f'exist for {self._adata.n_obs} samples. It’s ignored.'
                )
            else:
                self.iroot = self._adata.uns['iroot']
            return
        # set iroot via xroot
        xroot = None
        if 'xroot' in self._adata.uns:
            xroot = self._adata.uns['xroot']
        elif 'xroot' in self._adata.var:
            xroot = self._adata.var['xroot']
        # see whether we can set self.iroot using the full data matrix
        if xroot is not None and xroot.size == self._adata.shape[1]:
            self._set_iroot_via_xroot(xroot)

    def _get_dpt_row(self, i: int) -> np.ndarray:
        mask = None
        if self._number_connected_components > 1:
            label = self._connected_components[1][i]
            mask = self._connected_components[1] == label
        row = sum(
            (
                self.eigen_values[j]
                / (1 - self.eigen_values[j])
                * (self.eigen_basis[i, j] - self.eigen_basis[:, j])
            )
            ** 2
            # account for float32 precision
            for j in range(0, self.eigen_values.size)
            if self.eigen_values[j] < 0.9994
        )
        # thanks to Marius Lange for pointing Alex to this:
        # we will likely remove the contributions from the stationary state below when making
        # backwards compat breaking changes, they originate from an early implementation in 2015
        # they never seem to have deteriorated results, but also other distance measures (see e.g.
        # PAGA paper) don't have it, which makes sense
        row += sum(
            (self.eigen_basis[i, k] - self.eigen_basis[:, k]) ** 2
            for k in range(0, self.eigen_values.size)
            if self.eigen_values[k] >= 0.9994
        )
        if mask is not None:
            row[~mask] = np.inf
        return np.sqrt(row)

    def _set_pseudotime(self):
        """Return pseudotime with respect to root point."""
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
                'The root vector you provided does not have the ' 'correct dimension.'
            )
        # this is the squared distance
        dsqroot = 1e10
        iroot = 0
        for i in range(self._adata.shape[0]):
            diff = self._adata.X[i, :] - xroot
            dsq = diff @ diff
            if dsq < dsqroot:
                dsqroot = dsq
                iroot = i
                if np.sqrt(dsqroot) < 1e-10:
                    break
        logg.debug(f'setting root index to {iroot}')
        if self.iroot is not None and iroot != self.iroot:
            logg.warning(f'Changing index of iroot from {self.iroot} to {iroot}.')
        self.iroot = iroot
