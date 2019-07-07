import pandas as pd
import numpy as np
from umap import UMAP
from umap.distances import named_distances
from umap.nndescent import make_initialisations, make_initialized_nnd_search
from umap.umap_ import INT32_MAX, INT32_MIN
from umap.utils import deheap_sort
from sklearn.utils import check_random_state
from scipy.sparse import issparse

from ..preprocessing._simple import N_PCS
from ..neighbors import _rp_forest_generate


class Ingest:
    def __init__(self, adata):
        #assume rep is X if all initializations fail to identify it
        self._rep = adata.X
        self._use_rep = 'X'

        self._n_pcs = None

        #maybe don't need it, rather initialize Ingest class with all needed cluster keys
        self._adata = adata

        if 'PCs' in adata.varm:
            self._pca_basis = adata.varm['PCs']

        if 'neighbors' not in adata.uns:
            return

        #need to split all these into separate init functions
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

        if 'X_umap' not in adata.obsm:
            return

        self._umap = UMAP(
            metric = adata.uns['neighbors']['params']['metric']
        )

        self._umap.embedding_ = adata.obsm['X_umap']
        self._umap._raw_data = self._rep
        self._umap._sparse_data = issparse(self._rep)
        self._umap._small_data = self._rep.shape[0] < 4096
        self._umap._metric_kwds = adata.uns['neighbors']['params']['metric_kwds']
        self._umap._n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
        self._umap._initial_alpha = self._umap.learning_rate

        self._umap._random_init = self._random_init
        self._umap._tree_init = self._tree_init
        self._umap._search = self._search

        self._umap._rp_forest = self._rp_forest

        self._umap._search_graph = self._search_graph

        self._umap._a = adata.uns['umap']['params']['a']
        self._umap._b = adata.uns['umap']['params']['b']

    def pca(self, adata_small, n_pcs=None, inplace=True):
        #todo - efficient implementation for sparse matrices
        rep = adata_small.X
        rep = rep.toarray() if issparse(rep) else rep.copy()
        rep -= rep.mean(axis=0)
        X_pca = np.dot(rep, self._pca_basis[:, :n_pcs])
        if inplace:
            adata_small.obsm['X_pca'] = X_pca
        else:
            return X_pca

    def same_rep(self, adata_small):
        if self._n_pcs is not None:
            return self.pca(adata_small, self._n_pcs, inplace=False)
        if self._use_rep == 'X':
            return adata_small.X
        if self._use_rep in adata.obsm.keys():
            return adata.obsm[self._use_rep]
        return adata_small.X

    def neighbors(self, adata_small, k=10, queue_size=5, random_state=0):
        random_state = check_random_state(random_state)
        rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

        train = self._rep
        test = self.same_rep(adata_small)

        init = initialise_search(rp_forest, train, test, int(k * queue_size),
                                 self._random_init, self._tree_init, rng_state)

        result = self._search(train, self._search_graph.indptr, self._search_graph.indices, init, test)
        indices, dists = deheap_sort(result)
        return indices[:, :k], dists[:, :k]

    def umap(self, adata_small):
        rep = self.same_rep(adata_small)
        adata_small.obsm['X_umap'] = self._umap.transform(rep)

    def knn_classify(self, adata_small, classes_key, k=10, **kwargs):
        #i.e. ingest.knn_classify(adata_small, 'louvain')
        cat_array = self._adata.obs[classes_key]
        neighbors, _ = self.neighbors(adata_small, k, **kwargs)

        values = [cat_array[inds].mode()[0] for inds in neighbors]
        adata_small.obs[classes_key] = pd.Categorical(values=values, categories=cat_array.cat.categories)
