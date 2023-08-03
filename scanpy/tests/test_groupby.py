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
    pairs = [("v", "b"), ("v", "a"), ("w", "c"), ("w", "d"), ("w", "v")]

    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    obs["key"] = pd.Categorical([c[0] for c in cells])
    obs["tuple_key"] = pd.Categorical([(c[0],) for c in cells])
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
    stats_pd = sc.get.GroupBy(adata_dense, data=(adata_dense.layers['test'] if use_layers else adata_dense.X), key="key").pd_count_mean_var()

    assert stats_sparse.count.equals(stats_dense.count)
    assert np.allclose(stats_sparse.mean, stats_dense.mean)
    assert np.allclose(stats_sparse.var, stats_dense.var, equal_nan=True)

    assert stats_sparse.count.equals(stats_pd.count["A"])
    assert np.allclose(stats_sparse.mean, stats_pd.mean)
    assert np.allclose(stats_sparse.var, stats_pd.var, equal_nan=True)

    gb_weight = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="key", weight="weight")
    stats_weight = gb_weight.count_mean_var()
    sum_ = gb.sum()
    sum_weight = gb_weight.sum()

    assert np.allclose(2 * sum_, sum_weight)
    assert np.allclose(stats_sparse.mean, stats_weight.mean)
    assert np.allclose(stats_sparse.var, stats_dense.var, equal_nan=True)

    key_set = ["v", "w"]
    mean_key_set = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="key", key_set=key_set).mean()
    assert np.allclose(stats_sparse.mean.loc[key_set], mean_key_set)

    gb_explode = sc.get.GroupBy(adata_sparse, data=(adata_sparse.layers['test'] if use_layers else adata_sparse.X), key="tuple_key", explode=True)
    stats_explode = gb_explode.count_mean_var()

    assert stats_sparse.count.equals(stats_explode.count)
    assert np.allclose(stats_sparse.mean, stats_explode.mean)
    assert np.allclose(stats_sparse.var, stats_explode.var, equal_nan=True)
