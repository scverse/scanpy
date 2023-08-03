import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import pytest

@pytest.mark.parametrize(
    'use_layers',
    [
        False,
        True,
    ],
)
def test_groupby(use_layers):
    genes = ["A", "B"]
    cells = [
        "v0",
        "v1",
        "v2",
        "w0",
        "w1",
        "a1",
        "a2",
        "a3",
        "b1",
        "b2",
        "c1",
        "c2",
        "d0",
    ]

    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    obs["key"] = pd.Categorical([c[0] for c in cells])
    obs["key_superset"] = pd.Categorical([c[0] for c in cells]).map({'v': 'v', 'w': 'v', 'a': 'a', 'b': 'a', 'c': 'a', 'd': 'a'})
    obs["key_subset"] = pd.Categorical([c[1] for c in cells])
    obs["weight"] = 2.0

    var = pd.DataFrame(index=genes)

    X = np.array(
        [
            [0, -2],
            [1, 13],
            [2, 1],  # v
            [3, 12],
            [4, 2],  # w
            [5, 11],
            [6, 3],
            [7, 10],  # a
            [8, 4],
            [9, 9],  # b
            [10, 5],
            [11, 8],  # c
            [12, 6],  # d
        ],
        dtype=np.float32,
    )

    adata_sparse = ad.AnnData(obs=obs, var=var, X=csr_matrix(X), layers={ 'test': csr_matrix(X) })
    adata_dense = ad.AnnData(obs=obs, var=var, X=X, layers={ 'test': X.copy() }) # .copy needed?

    gb = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="key")
    stats_sparse = gb.count_mean_var()
    stats_dense = sc.get.GroupBy(adata_dense, data=(adata_dense.layers['test'] if use_layers else adata_dense.X), key="key").count_mean_var()

    # superset columns can be kept but not subsets
    assert 'key_superset' in stats_dense.obs
    assert 'key_subset' not in stats_dense.obs

    assert np.allclose(stats_sparse.obs['count'], stats_dense.obs['count'])
    assert np.allclose(stats_sparse.layers['mean'], stats_dense.layers['mean'])
    assert np.allclose(stats_sparse.layers['var'], stats_dense.layers['var'], equal_nan=True)

    gb_weight = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="key", weight="weight")
    stats_weight = gb_weight.count_mean_var()
    sum_ = gb.sum()
    sum_weight = gb_weight.sum()

    assert np.allclose(2 * sum_.X, sum_weight.X)
    assert np.allclose(stats_sparse.layers['mean'], stats_weight.layers['mean'])
    assert np.allclose(stats_sparse.layers['var'], stats_dense.layers['var'], equal_nan=True)

    key_set = ["v", "w"]
    mean_key_set = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="key", key_set=key_set).mean()
    assert np.allclose(stats_sparse[stats_sparse.obs.index.isin(key_set), :].layers['mean'], mean_key_set.X)
    