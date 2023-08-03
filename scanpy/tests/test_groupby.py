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
@pytest.mark.parametrize(
    'df_key',
    [
        'obs',
        'var',
    ],
)
def test_groupby(use_layers, df_key):
    ax_base = ["A", "B"]
    ax_groupby = [
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

    df_groupby = pd.DataFrame(index=pd.Index(ax_groupby, name="cell"))
    df_groupby["key"] = pd.Categorical([c[0] for c in ax_groupby])
    df_groupby["key_superset"] = pd.Categorical([c[0] for c in ax_groupby]).map({'v': 'v', 'w': 'v', 'a': 'a', 'b': 'a', 'c': 'a', 'd': 'a'})
    df_groupby["key_subset"] = pd.Categorical([c[1] for c in ax_groupby])
    df_groupby["weight"] = 2.0

    df_base = pd.DataFrame(index=ax_base)

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
    if df_key == 'obs':
        adata_sparse = ad.AnnData(obs=df_groupby, var=df_base, X=csr_matrix(X), layers={ 'test': csr_matrix(X) })
        adata_dense = ad.AnnData(obs=df_groupby, var=df_base, X=X, layers={ 'test': X.copy() }) # .copy needed?
    else:
        adata_sparse = ad.AnnData(obs=df_base, var=df_groupby, X=csr_matrix(X.T), layers={ 'test': csr_matrix(X.T) })
        adata_dense = ad.AnnData(obs=df_base, var=df_groupby, X=X.T, layers={ 'test': X.T.copy() }) # .copy needed?

    data_sparse = adata_sparse.layers['test'] if use_layers else None
    if df_key == 'var' and use_layers:
        data_sparse = data_sparse.T
    gb = sc.get.GroupBy(adata_sparse, key="key", data=data_sparse, df_key=df_key)
    stats_sparse = gb.count_mean_var()
    data_dense = adata_dense.layers['test'] if use_layers else None
    if df_key == 'var' and use_layers:
        data_dense = data_dense.T
    stats_dense = sc.get.GroupBy(adata_dense, key="key", data=data_dense, df_key=df_key).count_mean_var()

    # superset columns can be kept but not subsets
    assert 'key_superset' in getattr(stats_sparse, df_key)
    assert 'key_subset' not in getattr(stats_sparse, df_key)

    assert np.allclose(getattr(stats_sparse, df_key)['count'], getattr(stats_sparse, df_key)['count'])
    assert np.allclose(stats_sparse.layers['mean'], stats_dense.layers['mean'])
    assert np.allclose(stats_sparse.layers['var'], stats_dense.layers['var'], equal_nan=True)
    gb_weight = sc.get.GroupBy(adata_sparse, key="key", data=data_sparse, weight="weight", df_key=df_key)
    stats_weight = gb_weight.count_mean_var()
    sum_ = gb.sum()
    sum_weight = gb_weight.sum()

    assert np.allclose(2 * sum_.X, sum_weight.X)
    assert np.allclose(stats_sparse.layers['mean'], stats_weight.layers['mean'])
    assert np.allclose(stats_sparse.layers['var'], stats_dense.layers['var'], equal_nan=True)

    key_set = ["v", "w"]
    mean_key_set = sc.get.GroupBy(adata_sparse, key="key", data=data_sparse, key_set=key_set,  df_key=df_key).mean()
    subset_idx = getattr(stats_sparse, df_key).index.isin(key_set)
    subset = stats_sparse[subset_idx, :] if df_key == 'obs' else stats_sparse[:, subset_idx]
    assert np.allclose(subset.layers['mean'], mean_key_set.X)

    df = pd.DataFrame(
        index=getattr(adata_dense, df_key)["key"],
        columns=getattr(adata_dense, f"{'var' if df_key == 'obs' else 'obs'}_names"),
        data=adata_dense.X if df_key == 'obs' else adata_dense.X.T,
    )
    grouped_agg_df = df.groupby('key').agg(["count", "mean", "var"]).swaplevel(axis=1).sort_index(axis=1)
    mean = stats_dense.layers['mean']
    if df_key == 'var':
        mean = mean.T
    assert np.allclose(mean, grouped_agg_df['mean'].values)
    var = stats_dense.layers['var']
    if df_key == 'var':
        var = var.T
    assert np.allclose(var, grouped_agg_df['var'].values, equal_nan=True)
    assert np.allclose(getattr(stats_dense, df_key)['count'], grouped_agg_df['count']['A'].values) # returns for both columns but counts only needs one because it is the same
    