from collections.abc import MutableMapping
from typing import Iterable, Union, Optional

import pandas as pd
import numpy as np
from packaging import version
from sklearn.utils import check_random_state
from scipy.sparse import issparse
from anndata import AnnData

from .. import settings
from .. import logging as logg
from ..neighbors import _rp_forest_generate
from .._utils import NeighborsView
from .._compat import pkg_version

ANNDATA_MIN_VERSION = version.parse("0.7rc1")


def ingest(
    adata: AnnData,
    adata_ref: AnnData,
    obs: Optional[Union[str, Iterable[str]]] = None,
    embedding_method: Union[str, Iterable[str]] = ('umap', 'pca'),
    labeling_method: str = 'knn',
    neighbors_key: Optional[str] = None,
    inplace: bool = True,
    **kwargs,
):
    """\
    Map labels and embeddings from reference data to new data.

    :tutorial:`integrating-data-using-ingest`

    Integrates embeddings and annotations of an `adata` with a reference dataset
    `adata_ref` through projecting on a PCA (or alternate
    model) that has been fitted on the reference data. The function uses a knn
    classifier for mapping labels and the UMAP package [McInnes18]_ for mapping
    the embeddings.

    .. note::

        We refer to this *asymmetric* dataset integration as *ingesting*
        annotations from reference data to new data. This is different from
        learning a joint representation that integrates both datasets in an
        unbiased way, as CCA (e.g. in Seurat) or a conditional VAE (e.g. in
        scVI) would do.

    You need to run :func:`~scanpy.pp.neighbors` on `adata_ref` before
    passing it.

    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes. This is the dataset without labels and
        embeddings.
    adata_ref
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
        Variables (`n_vars` and `var_names`) of `adata_ref` should be the same
        as in `adata`.
        This is the dataset with labels and embeddings
        which need to be mapped to `adata`.
    obs
        Labels' keys in `adata_ref.obs` which need to be mapped to `adata.obs`
        (inferred for observation of `adata`).
    embedding_method
        Embeddings in `adata_ref` which need to be mapped to `adata`.
        The only supported values are 'umap' and 'pca'.
    labeling_method
        The method to map labels in `adata_ref.obs` to `adata.obs`.
        The only supported value is 'knn'.
    neighbors_key
        If not specified, ingest looks adata_ref.uns['neighbors']
        for neighbors settings and adata_ref.obsp['distances'] for
        distances (default storage places for pp.neighbors).
        If specified, ingest looks adata_ref.uns[neighbors_key] for
        neighbors settings and
        adata_ref.obsp[adata_ref.uns[neighbors_key]['distances_key']] for distances.
    inplace
        Only works if `return_joint=False`.
        Add labels and embeddings to the passed `adata` (if `True`)
        or return a copy of `adata` with mapped embeddings and labels.

    Returns
    -------
    * if `inplace=False` returns a copy of `adata`
      with mapped embeddings and labels in `obsm` and `obs` correspondingly
    * if `inplace=True` returns `None` and updates `adata.obsm` and `adata.obs`
      with mapped embeddings and labels

    Example
    -------
    Call sequence:

    >>> import scanpy as sc
    >>> sc.pp.neighbors(adata_ref)
    >>> sc.tl.umap(adata_ref)
    >>> sc.tl.ingest(adata, adata_ref, obs='cell_type')

    .. _ingest PBMC tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/integrating-pbmcs-using-ingest.html
    .. _ingest Pancreas tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/integrating-pancreas-using-ingest.html
    """
    # anndata version check
    anndata_version = pkg_version("anndata")
    if anndata_version < ANNDATA_MIN_VERSION:
        raise ValueError(
            f'ingest only works correctly with anndata>={ANNDATA_MIN_VERSION} '
            f'(you have {anndata_version}) as prior to {ANNDATA_MIN_VERSION}, '
            '`AnnData.concatenate` did not concatenate `.obsm`.'
        )

    start = logg.info('running ingest')
    obs = [obs] if isinstance(obs, str) else obs
    embedding_method = (
        [embedding_method] if isinstance(embedding_method, str) else embedding_method
    )
    labeling_method = (
        [labeling_method] if isinstance(labeling_method, str) else labeling_method
    )

    if len(labeling_method) == 1 and len(obs or []) > 1:
        labeling_method = labeling_method * len(obs)

    ing = Ingest(adata_ref, neighbors_key)
    ing.fit(adata)

    for method in embedding_method:
        ing.map_embedding(method)

    if obs is not None:
        ing.neighbors(**kwargs)
        for i, col in enumerate(obs):
            ing.map_labels(col, labeling_method[i])

    logg.info('    finished', time=start)
    return ing.to_adata(inplace)


class _DimDict(MutableMapping):
    def __init__(self, dim, axis=0, vals=None):
        self._data = {}
        self._dim = dim
        self._axis = axis
        if vals is not None:
            self.update(vals)

    def __setitem__(self, key, value):
        if value.shape[self._axis] != self._dim:
            raise ValueError(
                f"Value passed for key '{key}' is of incorrect shape. "
                f"Value has shape {value.shape[self._axis]} "
                f"for dimension {self._axis} while "
                f"it should have {self._dim}."
            )
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return f"{type(self).__name__}({self._data})"


class Ingest:
    """\
    Class to map labels and embeddings from existing data to new data.

    You need to run :func:`~scanpy.pp.neighbors` on `adata` before
    initializing Ingest with it.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix of shape `n_obs` × `n_vars`
        with embeddings and labels.
    """

    def _init_umap(self, adata):
        import umap as u

        if not self._use_pynndescent:
            u.umap_._HAVE_PYNNDESCENT = False

        self._umap = u.UMAP(
            metric=self._metric,
            random_state=adata.uns['umap']['params'].get('random_state', 0),
        )

        self._umap._initial_alpha = self._umap.learning_rate
        self._umap._raw_data = self._rep
        self._umap.knn_dists = None

        self._umap._validate_parameters()

        self._umap.embedding_ = adata.obsm['X_umap']
        self._umap._sparse_data = issparse(self._rep)
        self._umap._small_data = self._rep.shape[0] < 4096
        self._umap._metric_kwds = self._metric_kwds

        self._umap._n_neighbors = self._n_neighbors
        self._umap.n_neighbors = self._n_neighbors

        if not self._use_pynndescent:
            if self._random_init is not None or self._tree_init is not None:
                self._umap._random_init = self._random_init
                self._umap._tree_init = self._tree_init
                self._umap._search = self._search

            if self._dist_func is not None:
                self._umap._input_distance_func = self._dist_func

            self._umap._rp_forest = self._rp_forest

            self._umap._search_graph = self._search_graph
        else:
            self._umap._knn_search_index = self._nnd_idx

        self._umap._a = adata.uns['umap']['params']['a']
        self._umap._b = adata.uns['umap']['params']['b']

        self._umap._input_hash = None

    def _init_dist_search(self, dist_args):
        from functools import partial
        from umap.nndescent import initialise_search
        from umap.distances import named_distances

        self._random_init = None
        self._tree_init = None

        self._initialise_search = None
        self._search = None

        self._dist_func = None

        dist_func = named_distances[self._metric]

        if pkg_version('umap-learn') < version.parse("0.4.0"):
            from umap.nndescent import (
                make_initialisations,
                make_initialized_nnd_search,
            )

            self._random_init, self._tree_init = make_initialisations(
                dist_func, dist_args
            )
            _initialise_search = partial(
                initialise_search,
                init_from_random=self._random_init,
                init_from_tree=self._tree_init,
            )
            _search = make_initialized_nnd_search(dist_func, dist_args)

        else:
            from numba import njit
            from umap.nndescent import initialized_nnd_search

            @njit
            def partial_dist_func(x, y):
                return dist_func(x, y, *dist_args)

            _initialise_search = partial(initialise_search, dist=partial_dist_func)
            _search = partial(initialized_nnd_search, dist=partial_dist_func)

            self._dist_func = partial_dist_func

        self._initialise_search = _initialise_search
        self._search = _search

    def _init_pynndescent(self, distances):
        from pynndescent import NNDescent

        self._use_pynndescent = True

        first_col = np.arange(distances.shape[0])[:, None]
        init_indices = np.hstack((first_col, np.stack(distances.tolil().rows)))

        self._nnd_idx = NNDescent(
            data=self._rep,
            metric=self._metric,
            metric_kwds=self._metric_kwds,
            n_neighbors=self._n_neighbors,
            init_graph=init_indices,
            random_state=self._neigh_random_state,
        )

        # temporary hack for the broken forest storage
        from pynndescent.rp_trees import make_forest

        current_random_state = check_random_state(self._nnd_idx.random_state)
        self._nnd_idx._rp_forest = make_forest(
            self._nnd_idx._raw_data,
            self._nnd_idx.n_neighbors,
            self._nnd_idx.n_search_trees,
            self._nnd_idx.leaf_size,
            self._nnd_idx.rng_state,
            current_random_state,
            self._nnd_idx.n_jobs,
            self._nnd_idx._angular_trees,
        )

    def _init_neighbors(self, adata, neighbors_key):
        neighbors = NeighborsView(adata, neighbors_key)

        self._n_neighbors = neighbors['params']['n_neighbors']

        if 'use_rep' in neighbors['params']:
            self._use_rep = neighbors['params']['use_rep']
            self._rep = adata.X if self._use_rep == 'X' else adata.obsm[self._use_rep]
        elif 'n_pcs' in neighbors['params']:
            self._use_rep = 'X_pca'
            self._n_pcs = neighbors['params']['n_pcs']
            self._rep = adata.obsm['X_pca'][:, : self._n_pcs]
        elif adata.n_vars > settings.N_PCS and 'X_pca' in adata.obsm.keys():
            self._use_rep = 'X_pca'
            self._rep = adata.obsm['X_pca'][:, : settings.N_PCS]
            self._n_pcs = self._rep.shape[1]

        if 'metric_kwds' in neighbors['params']:
            self._metric_kwds = neighbors['params']['metric_kwds']
            dist_args = tuple(self._metric_kwds.values())
        else:
            self._metric_kwds = {}
            dist_args = ()

        self._metric = neighbors['params']['metric']

        if pkg_version('umap-learn') < version.parse("0.5.0"):
            self._init_dist_search(dist_args)

            search_graph = neighbors['distances'].copy()
            search_graph.data = (search_graph.data > 0).astype(np.int8)
            self._search_graph = search_graph.maximum(search_graph.transpose())

            if 'rp_forest' in neighbors:
                self._rp_forest = _rp_forest_generate(neighbors['rp_forest'])
            else:
                self._rp_forest = None
        else:
            self._neigh_random_state = neighbors['params'].get('random_state', 0)
            self._init_pynndescent(neighbors['distances'])

    def _init_pca(self, adata):
        self._pca_centered = adata.uns['pca']['params']['zero_center']
        self._pca_use_hvg = adata.uns['pca']['params']['use_highly_variable']

        if self._pca_use_hvg and 'highly_variable' not in adata.var.keys():
            raise ValueError('Did not find adata.var[\'highly_variable\'].')

        if self._pca_use_hvg:
            self._pca_basis = adata.varm['PCs'][adata.var['highly_variable']]
        else:
            self._pca_basis = adata.varm['PCs']

    def __init__(self, adata, neighbors_key=None):
        # assume rep is X if all initializations fail to identify it
        self._rep = adata.X
        self._use_rep = 'X'

        self._n_pcs = None

        self._adata_ref = adata
        self._adata_new = None

        # only use with umap > 0.5.0
        self._use_pynndescent = False

        if 'pca' in adata.uns:
            self._init_pca(adata)

        if neighbors_key is None:
            neighbors_key = 'neighbors'

        if neighbors_key in adata.uns:
            self._init_neighbors(adata, neighbors_key)
        else:
            raise ValueError(
                f'There is no neighbors data in `adata.uns["{neighbors_key}"]`.\n'
                'Please run pp.neighbors.'
            )

        if 'X_umap' in adata.obsm:
            self._init_umap(adata)

        self._obsm = None
        self._obs = None
        self._labels = None

        self._indices = None
        self._distances = None

    def _pca(self, n_pcs=None):
        X = self._adata_new.X
        X = X.toarray() if issparse(X) else X.copy()
        if self._pca_use_hvg:
            X = X[:, self._adata_ref.var['highly_variable']]
        if self._pca_centered:
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

    def fit(self, adata_new):
        """\
        Map `adata_new` to the same representation as `adata`.

        This function identifies the representation which was used to
        calculate neighbors in 'adata' and maps `adata_new` to
        this representation.
        Variables (`n_vars` and `var_names`) of `adata_new` should be the same
        as in `adata`.

        `adata` refers to the :class:`~anndata.AnnData` object
        that is passed during the initialization of an Ingest instance.
        """
        ref_var_names = self._adata_ref.var_names.str.upper()
        new_var_names = adata_new.var_names.str.upper()

        if not ref_var_names.equals(new_var_names):
            raise ValueError(
                'Variables in the new adata are different '
                'from variables in the reference adata'
            )

        self._obs = pd.DataFrame(index=adata_new.obs.index)
        self._obsm = _DimDict(adata_new.n_obs, axis=0)

        self._adata_new = adata_new
        self._obsm['rep'] = self._same_rep()

    def neighbors(self, k=None, queue_size=5, epsilon=0.1, random_state=0):
        """\
        Calculate neighbors of `adata_new` observations in `adata`.

        This function calculates `k` neighbors in `adata` for
        each observation of `adata_new`.
        """
        from umap.umap_ import INT32_MAX, INT32_MIN

        random_state = check_random_state(random_state)
        rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

        train = self._rep
        test = self._obsm['rep']

        if k is None:
            k = self._n_neighbors

        if self._use_pynndescent:
            self._nnd_idx.search_rng_state = rng_state

            self._indices, self._distances = self._nnd_idx.query(test, k, epsilon)

        else:
            from umap.utils import deheap_sort

            init = self._initialise_search(
                self._rp_forest, train, test, int(k * queue_size), rng_state=rng_state
            )

            result = self._search(
                train, self._search_graph.indptr, self._search_graph.indices, init, test
            )
            indices, dists = deheap_sort(result)
            self._indices, self._distances = indices[:, :k], dists[:, :k]

    def _umap_transform(self):
        return self._umap.transform(self._obsm['rep'])

    def map_embedding(self, method):
        """\
        Map embeddings of `adata` to `adata_new`.

        This function infers embeddings, specified by `method`,
        for `adata_new` from existing embeddings in `adata`.
        `method` can be 'umap' or 'pca'.
        """
        if method == 'umap':
            self._obsm['X_umap'] = self._umap_transform()
        elif method == 'pca':
            self._obsm['X_pca'] = self._pca()
        else:
            raise NotImplementedError(
                'Ingest supports only umap and pca embeddings for now.'
            )

    def _knn_classify(self, labels):
        cat_array = self._adata_ref.obs[labels].astype(
            'category'
        )  # ensure it's categorical
        values = [cat_array[inds].mode()[0] for inds in self._indices]
        return pd.Categorical(values=values, categories=cat_array.cat.categories)

    def map_labels(self, labels, method):
        """\
        Map labels of `adata` to `adata_new`.

        This function infers `labels` for `adata_new.obs`
        from existing labels in `adata.obs`.
        `method` can be only 'knn'.
        """
        if method == 'knn':
            self._obs[labels] = self._knn_classify(labels)
        else:
            raise NotImplementedError('Ingest supports knn labeling for now.')

    def to_adata(self, inplace=False):
        """\
        Returns `adata_new` with mapped embeddings and labels.

        If `inplace=False` returns a copy of `adata_new`
        with mapped embeddings and labels in `obsm` and `obs` correspondingly.
        If `inplace=True` returns nothing and updates `adata_new.obsm`
        and `adata_new.obs` with mapped embeddings and labels.
        """
        adata = self._adata_new if inplace else self._adata_new.copy()

        adata.obsm.update(self._obsm)

        for key in self._obs:
            adata.obs[key] = self._obs[key]

        if not inplace:
            return adata

    def to_adata_joint(
        self, batch_key='batch', batch_categories=None, index_unique='-'
    ):
        """\
        Returns concatenated object.

        This function returns the new :class:`~anndata.AnnData` object
        with concatenated existing embeddings and labels of 'adata'
        and inferred embeddings and labels for `adata_new`.
        """
        adata = self._adata_ref.concatenate(
            self._adata_new,
            batch_key=batch_key,
            batch_categories=batch_categories,
            index_unique=index_unique,
        )

        obs_update = self._obs.copy()
        obs_update.index = adata[adata.obs[batch_key] == '1'].obs_names
        adata.obs.update(obs_update)

        for key in self._obsm:
            if key in self._adata_ref.obsm:
                adata.obsm[key] = np.vstack(
                    (self._adata_ref.obsm[key], self._obsm[key])
                )

        if self._use_rep not in ('X_pca', 'X'):
            adata.obsm[self._use_rep] = np.vstack(
                (self._adata_ref.obsm[self._use_rep], self._obsm['rep'])
            )

        if 'X_umap' in self._obsm:
            adata.uns['umap'] = self._adata_ref.uns['umap']
        if 'X_pca' in self._obsm:
            adata.uns['pca'] = self._adata_ref.uns['pca']
            adata.varm['PCs'] = self._adata_ref.varm['PCs']

        return adata
