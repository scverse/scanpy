from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import scanpy as sc
from scanpy.testing._helpers import assert_equal
from scanpy.testing._helpers.data import pbmc3k_processed


@pytest.fixture
def df_base():
    ax_base = ["A", "B"]
    return pd.DataFrame(index=ax_base)


@pytest.fixture
def df_groupby():
    ax_groupby = [
        *["v0", "v1", "v2"],
        *["w0", "w1"],
        *["a1", "a2", "a3"],
        *["b1", "b2"],
        *["c1", "c2"],
        "d0",
    ]

    df_groupby = pd.DataFrame(index=pd.Index(ax_groupby, name="cell"))
    df_groupby["key"] = pd.Categorical([c[0] for c in ax_groupby])
    df_groupby["key_superset"] = pd.Categorical([c[0] for c in ax_groupby]).map(
        {"v": "v", "w": "v", "a": "a", "b": "a", "c": "a", "d": "a"}
    )
    df_groupby["key_subset"] = pd.Categorical([c[1] for c in ax_groupby])
    df_groupby["weight"] = 2.0
    return df_groupby


@pytest.fixture
def X():
    data = [
        *[[0, -2], [1, 13], [2, 1]],  # v
        *[[3, 12], [4, 2]],  # w
        *[[5, 11], [6, 3], [7, 10]],  # a
        *[[8, 4], [9, 9]],  # b
        *[[10, 5], [11, 8]],  # c
        [12, 6],  # d
    ]
    return np.array(data, dtype=np.float32)


def gen_adata(data_key, dim, df_base, df_groupby, X):
    if (data_key == "varm" and dim == "obs") or (data_key == "obsm" and dim == "var"):
        pytest.skip("invalid parameter combination")

    obs_df, var_df = (df_groupby, df_base) if dim == "obs" else (df_base, df_groupby)
    data = X.T if dim == "var" and data_key != "varm" else X
    if data_key != "X":
        data_dict_sparse = {data_key: {"test": csr_matrix(data)}}
        data_dict_dense = {data_key: {"test": data}}
    else:
        data_dict_sparse = {data_key: csr_matrix(data)}
        data_dict_dense = {data_key: data}

    adata_sparse = ad.AnnData(obs=obs_df, var=var_df, **data_dict_sparse)
    adata_dense = ad.AnnData(obs=obs_df, var=var_df, **data_dict_dense)
    return adata_sparse, adata_dense


# TODO: There isn't an exact equivalent for our count operation in pandas I think (i.e. count non-zero values)
@pytest.mark.parametrize("metric", ["sum", "mean", "var"])
def test_aggregated_vs_pandas(metric):
    adata = pbmc3k_processed().raw.to_adata()
    adata.obs["percent_mito_binned"] = pd.cut(adata.obs["percent_mito"], bins=10)
    result = sc.get.aggregated(adata, ["louvain", "percent_mito_binned"], metric)

    expected = (
        adata.to_df()
        .astype(np.float64)
        .join(adata.obs[["louvain", "percent_mito_binned"]])
        .groupby(["louvain", "percent_mito_binned"], observed=True)
        .agg(metric)
    )
    # TODO: figure out the axis names
    expected.index = expected.index.to_frame().apply(
        lambda x: "_".join(map(str, x)), axis=1
    )
    expected.index.name = None
    expected.columns.name = None

    result_df = result.to_df(layer=metric)
    result_df.index.name = None
    result_df.columns.name = None
    # expected = adata.to_df().groupby(adata.obs[["louvain", "percent_mito_binned"]]).agg(metric)

    assert_equal(result_df, expected)


@pytest.mark.parametrize("data_key", ["layers", "obsm", "varm", "X"])
@pytest.mark.parametrize("dim", ["obs", "var"])
def test_groupby(data_key, dim, df_base, df_groupby, X):
    adata_sparse, adata_dense = gen_adata(data_key, dim, df_base, df_groupby, X)

    data_loc_dict = (
        {(data_key if data_key != "layers" else "layer"): "test"}
        if data_key != "X"
        else {}
    )
    # When `X` is not the `data_key`, the multi-aggregation data is colocated with the `data_key`.  Otherwise it is in `layers`.
    multi_agg_data_loc_key = data_key if data_key != "X" else "layers"

    stats_sparse, stats_dense = (
        sc.get.aggregated(
            adata,
            by="key",
            dim=dim,
            func=["count", "mean", "var"],
            **data_loc_dict,
        )
        for adata in [adata_sparse, adata_dense]
    )

    # superset columns can be kept but not subsets
    assert "key_superset" in getattr(stats_sparse, dim)
    assert "key_subset" not in getattr(stats_sparse, dim)

    assert np.allclose(
        getattr(stats_sparse, dim)["count"],
        getattr(stats_dense, dim)["count"],
    )
    assert np.allclose(
        getattr(stats_sparse, multi_agg_data_loc_key)["mean"],
        getattr(stats_dense, multi_agg_data_loc_key)["mean"],
    )
    assert np.allclose(
        getattr(stats_sparse, multi_agg_data_loc_key)["var"],
        getattr(stats_dense, multi_agg_data_loc_key)["var"],
        equal_nan=True,
    )

    stats_weight = sc.get.aggregated(
        adata_dense,
        by="key",
        dim=dim,
        func=["count", "mean", "var"],
        weight_key="weight",
        **data_loc_dict,
    )
    sum_ = sc.get.aggregated(
        adata_sparse, by="key", dim=dim, func="sum", **data_loc_dict
    )
    sum_weight = sc.get.aggregated(
        adata_dense,
        by="key",
        dim=dim,
        func="sum",
        weight_key="weight",
        **data_loc_dict,
    )

    def get_single_agg(adata, key, agg):
        # Get the data of the aggregation from the correct location when only one `func` is passed in to `aggregated`
        if (key != "obsm" and key != "varm") or data_key == "X":
            return adata.X
        return getattr(adata, key)[agg]

    assert np.allclose(
        2 * get_single_agg(sum_, data_key, "sum"),
        get_single_agg(sum_weight, data_key, "sum"),
    )
    assert np.allclose(
        getattr(stats_sparse, multi_agg_data_loc_key)["mean"],
        getattr(stats_weight, multi_agg_data_loc_key)["mean"],
    )
    assert np.allclose(
        getattr(stats_sparse, multi_agg_data_loc_key)["var"],
        getattr(stats_dense, multi_agg_data_loc_key)["var"],
        equal_nan=True,
    )

    key_set = ["v", "w"]
    mean_key_set_adata = sc.get.aggregated(
        adata_dense,
        by="key",
        dim=dim,
        func="mean",
        key_set=key_set,
        **data_loc_dict,
    )
    subset_idx = getattr(stats_sparse, dim).index.isin(key_set)
    subset_adata = (
        stats_sparse[subset_idx, :] if dim == "obs" else stats_sparse[:, subset_idx]
    )
    subset_mean = getattr(subset_adata, multi_agg_data_loc_key)["mean"]
    key_set_mean = get_single_agg(mean_key_set_adata, data_key, "mean")

    assert np.allclose(subset_mean, key_set_mean)

    df = pd.DataFrame(
        index=getattr(adata_dense, dim)["key"],
        columns=getattr(adata_dense, f"{'var' if dim == 'obs' else 'obs'}_names"),
        data=X,
    )
    grouped_agg_df = (
        df.groupby("key")
        .agg(["count", "mean", "var"])
        .swaplevel(axis=1)
        .sort_index(axis=1)
    )
    mean = getattr(stats_dense, multi_agg_data_loc_key)["mean"]
    if dim == "var" and data_key != "varm":
        mean = mean.T
    assert np.allclose(mean, grouped_agg_df["mean"].values)
    var = getattr(stats_dense, multi_agg_data_loc_key)["var"]
    if dim == "var" and multi_agg_data_loc_key != "varm":
        var = var.T
    assert np.allclose(var, grouped_agg_df["var"].values, equal_nan=True)
    assert np.allclose(
        getattr(stats_dense, dim)["count"],
        grouped_agg_df["count"]["A"].values,
    )  # returns for both columns but counts only needs one because it is the same


@pytest.mark.parametrize(
    "matrix,df,keys,metrics,expected",
    [
        (
            np.block(
                [
                    [np.ones((2, 2)), np.zeros((2, 2))],
                    [np.zeros((2, 2)), np.ones((2, 2))],
                ]
            ),
            pd.DataFrame(
                {
                    "a": ["a", "a", "b", "b"],
                    "b": ["c", "d", "d", "d"],
                }
            ),
            ["a", "b"],
            ["count"],  # , "sum", "mean"],
            ad.AnnData(
                obs=pd.DataFrame(index=["a_c", "a_d", "b_d"]),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "count": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]),
                    # "sum": np.array([[2, 0], [0, 2]]),
                    # "mean": np.array([[1, 0], [0, 1]]),
                },
            ),
        ),
        (
            np.block(
                [
                    [np.ones((2, 2)), np.zeros((2, 2))],
                    [np.zeros((2, 2)), np.ones((2, 2))],
                ]
            ),
            pd.DataFrame(
                {
                    "a": ["a", "a", "b", "b"],
                    "b": ["c", "d", "d", "d"],
                }
            ),
            ["a", "b"],
            ["sum", "mean", "count"],
            ad.AnnData(
                obs=pd.DataFrame(index=["a_c", "a_d", "b_d"]),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "sum": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]),
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                    "count": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]),
                },
            ),
        ),
        (
            np.block(
                [
                    [np.ones((2, 2)), np.zeros((2, 2))],
                    [np.zeros((2, 2)), np.ones((2, 2))],
                ]
            ),
            pd.DataFrame(
                {
                    "a": ["a", "a", "b", "b"],
                    "b": ["c", "d", "d", "d"],
                }
            ),
            ["a", "b"],
            ["mean"],  # , "sum", "mean"],
            ad.AnnData(
                obs=pd.DataFrame(index=["a_c", "a_d", "b_d"]),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                    # "sum": np.array([[2, 0], [0, 2]]),
                    # "mean": np.array([[1, 0], [0, 1]]),
                },
            ),
        ),
    ],
)
def test_aggregated_parameterized(matrix, df, keys, metrics, expected):
    adata = ad.AnnData(
        X=matrix,
        obs=df,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(matrix.shape[1])]),
    )
    result = sc.get.aggregated(adata, by=keys, func=metrics)

    print(result)
    print(expected)

    for k in metrics:
        assert_equal(result.layers[k], expected.layers[k])

    assert_equal(expected, result)


@pytest.mark.parametrize(
    "label_df,cols,expected",
    [
        (
            pd.DataFrame(
                {
                    "a": pd.Categorical(["a", "b", "c"]),
                    "b": pd.Categorical(["d", "d", "f"]),
                }
            ),
            ["a", "b"],
            pd.Categorical(["a_d", "b_d", "c_f"]),
        ),
        (
            pd.DataFrame(
                {
                    "a": pd.Categorical(["a", "b", "c"]),
                    "b": pd.Categorical(["d", "d", "f"]),
                    "c": pd.Categorical(["g", "h", "h"]),
                }
            ),
            ["a", "b", "c"],
            pd.Categorical(["a_d_g", "b_d_h", "c_f_h"]),
        ),
        (
            pd.DataFrame(
                {
                    "a": pd.Categorical(["a", "b", "c"]),
                    "b": pd.Categorical(["d", "d", "f"]),
                    "c": pd.Categorical(["g", "h", "h"]),
                }
            ),
            ["a", "c"],
            pd.Categorical(["a_g", "b_h", "c_h"]),
        ),
        (
            pd.DataFrame(
                {
                    "a": pd.Categorical(["a", "b", "c"]),
                    "b": pd.Categorical(["d", "d", "f"]),
                    "c": pd.Categorical(["g", "h", "h"]),
                }
            ),
            ["b", "c"],
            pd.Categorical(["d_g", "d_h", "f_h"]),
        ),
    ],
)
def test_combine_categories(label_df, cols, expected):
    from scanpy.get._aggregated import _combine_categories

    result = _combine_categories(label_df, cols)

    assert isinstance(result, pd.Categorical)

    # TODO: is there a better function here?
    pd.testing.assert_series_equal(pd.Series(result), pd.Series(expected))
