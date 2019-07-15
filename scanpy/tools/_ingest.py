import pandas as pd
import numpy as np
from sklearn.utils import check_random_state
from scipy.sparse import issparse
from anndata.core.alignedmapping import AxisArrays

from ..preprocessing._simple import N_PCS
from ..neighbors import _rp_forest_generate


class Ingest:

    def _init_umap(self, adata):
        from umap import UMAP

        self._umap = UMAP(
            metric = adata.uns['neighbors']['params']['metric']
        )

        self._umap.embedding_ = adata.obsm['X_umap']
        self._umap._raw_data = self._rep
        self._umap._sparse_data = issparse(self._rep)
        self._umap._small_data = self._rep.shape[0] < 4096
        self._umap._metric_kwds = adata.uns['neighbors']['params'].get('metric_kwds', {})
        self._umap._n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
        self._umap._initial_alpha = self._umap.learning_rate

        self._umap._random_init = self._random_init
        self._umap._tree_init = self._tree_init
        self._umap._search = self._search

        self._umap._rp_forest = self._rp_forest

        self._umap._search_graph = self._search_graph

        self._umap._a = adata.uns['umap']['params']['a']
        self._umap._b = adata.uns['umap']['params']['b']

        self._umap._input_hash = None

    def _init_neighbors(self, adata):
        from umap.distances import named_distances
        from umap.nndescent import make_initialisations, make_initialized_nnd_search

        if 'use_rep' in adata.uns['neighbors']['params']:
            self._use_rep = adata.uns['neighbors']['params']['use_rep']
            self._rep = adata.X if self._use_rep == 'X' else adata.obsm[use_rep]
        elif 'n_pcs' in adata.uns['neighbors']['params']:
            self._use_rep = 'X_pca'
            self._n_pcs = adata.uns['neighbors']['params']['n_pcs']
            self._rep = adata.obsm['X_pca'][:, :self._n_pcs]
        elif adata.n_vars > N_PCS and 'X_pca' in adata.obsm.keys():
            self._use_rep = 'X_pca'
            self._rep = adata.obsm['X_pca'][:, :N_PCS]
            self._n_pcs = self._rep.shape[1]

        if 'metric_kwds' in adata.uns['neighbors']['params']:
            dist_args = tuple(adata.uns['neighbors']['params']['metric_kwds'].values())
        else:
            dist_args = ()
        dist_func = named_distances[adata.uns['neighbors']['params']['metric']]
        self._random_init, self._tree_init = make_initialisations(dist_func, dist_args)
        self._search = make_initialized_nnd_search(dist_func, dist_args)

        search_graph = adata.uns['neighbors']['distances'].copy()
        search_graph.data = (search_graph.data > 0).astype(np.int8)
        self._search_graph = search_graph.maximum(search_graph.transpose())

        if 'rp_forest' in adata.uns['neighbors']:
            self._rp_forest = _rp_forest_generate(adata.uns['neighbors']['rp_forest'])
        else:
            self._rp_forest = None

    def __init__(self, adata):
        #assume rep is X if all initializations fail to identify it
        self._rep = adata.X
        self._use_rep = 'X'

        self._n_pcs = None

        self._adata_ref = adata
        self._adata_new = None

        if 'PCs' in adata.varm:
            self._pca_basis = adata.varm['PCs']

        if 'neighbors' in adata.uns:
            self._init_neighbors(adata)

        if 'X_umap' in adata.obsm:
            self._init_umap(adata)

        self._obsm = None
        self._obs = None

        self._indices = None
        self._distances = None

    def _pca(self, n_pcs=None):
        #todo - efficient implementation for sparse matrices
        #todo - check if should be centered
        #todo - check if used highly_variable
        X = self._adata_new.X
        X = X.toarray() if issparse(X) else X.copy()
        X -= X.mean(axis=0)
        X_pca = np.dot(X, self._pca_basis[:, :n_pcs])
        return X_pca

    def _same_rep(self):
        adata = self._adata_new
        if self._n_pcs is not None:
            return self._pca(self._n_pcs)
        if self._use_rep == 'X':
            return adata.X
        if self._use_rep in adata.obsm.keys():
            return adata.obsm[self._use_rep]
        return adata.X

    def transform(self, adata_new):
        self._obs = pd.DataFrame(index=adata_new.obs.index)
        #not sure if it should be AxisArrays
        self._obsm = AxisArrays(adata_new, 0)

        self._adata_new = adata_new
        self._obsm['rep'] = self._same_rep()

    def neighbors(self, k=10, queue_size=5, random_state=0):
        from umap.nndescent import initialise_search
        from umap.utils import deheap_sort
        from umap.umap_ import INT32_MAX, INT32_MIN

        random_state = check_random_state(random_state)
        rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

        train = self._rep
        test = self._obsm['rep']

        init = initialise_search(self._rp_forest, train, test, int(k * queue_size),
                                 self._random_init, self._tree_init, rng_state)

        result = self._search(train, self._search_graph.indptr, self._search_graph.indices, init, test)
        indices, dists = deheap_sort(result)
        self._indices, self._distances = indices[:, :k], dists[:, :k]

    def _umap_transform(self):
        return self._umap.transform(self._obsm['rep'])

    def map_embedding(self, method):
        if method == 'umap':
            self._obsm['X_umap'] = self._umap_transform()
        else:
            raise NotImplementedError('Ingest supports only umap embeddings for now.')

    def map_labels(self, labels, k=10, **kwargs):
        cat_array = self._adata_ref.obs[labels]

        values = [cat_array[inds].mode()[0] for inds in self._indices]
        self._obs[labels] = pd.Categorical(values=values, categories=cat_array.cat.categories)

    def to_adata(self, inplace=False):
        adata = self._adata_new if inplace else self._adata_new.copy()
        adata.obsm.update(self._obsm)

        adata.obs.update(self._obs)
        new_cols = ~self._obs.columns.isin(adata.obs.columns)
        adata.obs = pd.concat([adata.obs, self._obs.loc[:, new_cols]], axis=1, copy=False)

        if not inplace:
            return adata
