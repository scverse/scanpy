import numpy as np
from umap import UMAP
from umap.distances import named_distances
from umap.nndescent import make_initialisations, make_initialized_nnd_search
from ..neighbors import _rp_forest_generate
from scipy.sparse import issparse

class Ingest:

    def __init__(self, adata):
        #need to take care of representation
        rep = adata.X

        #maybe don't need it, rather initialize Ingest class with all needed cluster keys
        self._adata = adata

        if 'PCs' in adata.varm:
            self._pca_basis = adata.varm['PCs']
        if 'X_umap' in adata.obsm:
            self._umap = UMAP(
                metric = adata.uns['neighbors']['params']['metric']
            )

            self._umap.embedding_ = adata.obsm['X_umap']
            self._umap._raw_data = rep
            self._umap._sparse_data = issparse(rep)
            self._umap._small_data = rep.shape[0] < 4096
            self._umap._metric_kwds = adata.uns['neighbors']['params']['metric_kwds']
            self._umap._n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
            self._umap._initial_alpha = self._umap.learning_rate

            dist_func = named_distances[self._umap.metric]
            dist_args = tuple(self._umap._metric_kwds.values())

            self._umap._random_init, self._umap._tree_init = make_initialisations(
                dist_func,
                dist_args
            )
            self._umap._search = make_initialized_nnd_search(dist_func, dist_args)

            if 'rp_forest' in adata.uns['neighbors']:
                self._umap._rp_forest = _rp_forest_generate(adata.uns['neighbors']['rp_forest'])
            else:
                self._umap._rp_forest = None

            search_graph = adata.uns['neighbors']['distances'].copy()
            search_graph.data = (search_graph.data > 0).astype(np.int8)
            search_graph = search_graph.maximum(search_graph.transpose())
            self._umap._search_graph = search_graph

            self._umap._a = adata.uns['umap']['params']['a']
            self._umap._b = adata.uns['umap']['params']['b']

    def umap(self, adata_small):
        #need to take care of representation
        rep = adata_small.X
        adata_small.obsm['X_umap'] = self._umap.transform(rep)

    def pca(self, adata_small):
        #todo - efficient implementation for sparse matrices
        rep = adata_small.X
        rep = rep.toarray() if issparse(rep) else rep.copy()
        rep -= rep.mean(axis=0)
        adata_small.obsm['X_pca'] = np.dot(rep, self._pca_basis)

    def knn_classify(self, adata_small, classes_key, k):
        #i.e. ingest.knn_classify(adata_small, 'louvain')
        cat_array = self._adata.obs[classes_key]
        #todo
