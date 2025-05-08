from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from packaging.version import Version
from scipy import sparse

import scanpy as sc
from scanpy._utils import _resolve_axis, get_literal_vals
from scanpy.get._aggregated import AggType
from testing.scanpy._helpers import assert_equal
from testing.scanpy._helpers.data import pbmc3k_processed
from testing.scanpy._pytest.params import ARRAY_TYPES_MEM


@pytest.fixture(params=get_literal_vals(AggType))
def metric(request: pytest.FixtureRequest) -> AggType:
    return request.param


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
        data_dict_sparse = {data_key: {"test": sparse.csr_matrix(data)}}  # noqa: TID251
        data_dict_dense = {data_key: {"test": data}}
    else:
        data_dict_sparse = {data_key: sparse.csr_matrix(data)}  # noqa: TID251
        data_dict_dense = {data_key: data}

    adata_sparse = ad.AnnData(obs=obs_df, var=var_df, **data_dict_sparse)
    adata_dense = ad.AnnData(obs=obs_df, var=var_df, **data_dict_dense)
    return adata_sparse, adata_dense


@pytest.mark.parametrize("axis", [0, 1])
def test_mask(axis):
    blobs = sc.datasets.blobs()
    mask = blobs.obs["blobs"] == 0
    blobs.obs["mask_col"] = mask
    if axis == 1:
        blobs = blobs.T
    by_name = sc.get.aggregate(blobs, "blobs", "sum", axis=axis, mask="mask_col")
    by_value = sc.get.aggregate(blobs, "blobs", "sum", axis=axis, mask=mask)

    assert_equal(by_name, by_value)

    assert np.all(by_name["0"].layers["sum"] == 0)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_aggregate_vs_pandas(metric, array_type):
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    adata.obs["percent_mito_binned"] = pd.cut(adata.obs["percent_mito"], bins=5)
    result = sc.get.aggregate(adata, ["louvain", "percent_mito_binned"], metric)

    if metric == "count_nonzero":
        expected = (
            (adata.to_df() != 0)
            .astype(np.float64)
            .join(adata.obs[["louvain", "percent_mito_binned"]])
            .groupby(["louvain", "percent_mito_binned"], observed=True)
            .agg("sum")
        )
    else:
        expected = (
            adata.to_df()
            .astype(np.float64)
            .join(adata.obs[["louvain", "percent_mito_binned"]])
            .groupby(["louvain", "percent_mito_binned"], observed=True)
            .agg(metric)
        )
    expected.index = expected.index.to_frame().apply(
        lambda x: "_".join(map(str, x)), axis=1
    )
    expected.index.name = None
    expected.columns.name = None

    result_df = result.to_df(layer=metric)
    result_df.index.name = None
    result_df.columns.name = None

    if Version(pd.__version__) < Version("2"):
        # Order of results returned by groupby changed in pandas 2
        assert expected.shape == result_df.shape
        assert expected.index.isin(result_df.index).all()

        expected = expected.loc[result_df.index]

    pd.testing.assert_frame_equal(result_df, expected, check_dtype=False, atol=1e-5)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_aggregate_axis(array_type, metric):
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    expected = sc.get.aggregate(adata, ["louvain"], metric)
    actual = sc.get.aggregate(adata.T, ["louvain"], metric, axis=1).T

    assert_equal(expected, actual)


def test_aggregate_entry():
    args = ("blobs", ["mean", "var", "count_nonzero"])

    adata = sc.datasets.blobs()
    X_result = sc.get.aggregate(adata, *args)
    # layer adata
    layer_adata = ad.AnnData(
        obs=adata.obs,
        var=adata.var,
        layers={"test": adata.X.copy()},
    )
    layer_result = sc.get.aggregate(layer_adata, *args, layer="test")
    obsm_adata = ad.AnnData(
        obs=adata.obs,
        var=adata.var,
        obsm={"test": adata.X.copy()},
    )
    obsm_result = sc.get.aggregate(obsm_adata, *args, obsm="test")
    varm_adata = ad.AnnData(
        obs=adata.var,
        var=adata.obs,
        varm={"test": adata.X.copy()},
    )
    varm_result = sc.get.aggregate(varm_adata, *args, varm="test")

    X_result_min = X_result.copy()
    del X_result_min.var
    X_result_min.var_names = [str(x) for x in np.arange(X_result_min.n_vars)]

    assert_equal(X_result, layer_result)
    assert_equal(X_result_min, obsm_result)
    assert_equal(X_result.layers, obsm_result.layers)
    assert_equal(X_result.layers, varm_result.T.layers)


def test_aggregate_incorrect_dim():
    adata = pbmc3k_processed().raw.to_adata()

    with pytest.raises(ValueError, match="was 'foo'"):
        sc.get.aggregate(adata, ["louvain"], "sum", axis="foo")


@pytest.mark.parametrize("axis_name", ["obs", "var"])
def test_aggregate_axis_specification(axis_name):
    axis, axis_name = _resolve_axis(axis_name)
    by = "blobs" if axis == 0 else "labels"

    adata = sc.datasets.blobs()
    adata.var["labels"] = np.tile(["a", "b"], adata.shape[1])[: adata.shape[1]]

    agg_index = sc.get.aggregate(adata, by=by, func="mean", axis=axis)
    agg_name = sc.get.aggregate(adata, by=by, func="mean", axis=axis_name)

    np.testing.assert_equal(agg_index.layers["mean"], agg_name.layers["mean"])

    if axis_name == "obs":
        agg_unspecified = sc.get.aggregate(adata, by=by, func="mean")
        np.testing.assert_equal(agg_name.layers["mean"], agg_unspecified.layers["mean"])


@pytest.mark.parametrize(
    ("matrix", "df", "keys", "metrics", "expected"),
    [
        pytest.param(
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
                },
                index=["a_c", "a_d", "b_d1", "b_d2"],
            ),
            ["a", "b"],
            ["count_nonzero"],  # , "sum", "mean"],
            ad.AnnData(
                obs=pd.DataFrame(
                    {"a": ["a", "a", "b"], "b": ["c", "d", "d"]},
                    index=["a_c", "a_d", "b_d"],
                ).astype("category"),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "count_nonzero": np.array(
                        [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]
                    ),
                    # "sum": np.array([[2, 0], [0, 2]]),
                    # "mean": np.array([[1, 0], [0, 1]]),
                },
            ),
            id="count_nonzero",
        ),
        pytest.param(
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
                },
                index=["a_c", "a_d", "b_d1", "b_d2"],
            ),
            ["a", "b"],
            ["sum", "mean", "count_nonzero"],
            ad.AnnData(
                obs=pd.DataFrame(
                    {"a": ["a", "a", "b"], "b": ["c", "d", "d"]},
                    index=["a_c", "a_d", "b_d"],
                ).astype("category"),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "sum": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]),
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                    "count_nonzero": np.array(
                        [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]
                    ),
                },
            ),
            id="sum-mean-count_nonzero",
        ),
        pytest.param(
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
                },
                index=["a_c", "a_d", "b_d1", "b_d2"],
            ),
            ["a", "b"],
            ["mean"],
            ad.AnnData(
                obs=pd.DataFrame(
                    {"a": ["a", "a", "b"], "b": ["c", "d", "d"]},
                    index=["a_c", "a_d", "b_d"],
                ).astype("category"),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                },
            ),
            id="mean",
        ),
    ],
)
def test_aggregate_examples(matrix, df, keys, metrics, expected):
    adata = ad.AnnData(
        X=matrix,
        obs=df,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(matrix.shape[1])]),
    )
    result = sc.get.aggregate(adata, by=keys, func=metrics)

    print(result)
    print(expected)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("label_cols", "cols", "expected"),
    [
        pytest.param(
            dict(
                a=pd.Categorical(["a", "b", "c"]),
                b=pd.Categorical(["d", "d", "f"]),
            ),
            ["a", "b"],
            pd.Categorical(["a_d", "b_d", "c_f"]),
            id="two_of_two",
        ),
        pytest.param(
            dict(
                a=pd.Categorical(["a", "b", "c"]),
                b=pd.Categorical(["d", "d", "f"]),
                c=pd.Categorical(["g", "h", "h"]),
            ),
            ["a", "b", "c"],
            pd.Categorical(["a_d_g", "b_d_h", "c_f_h"]),
            id="three_of_three",
        ),
        pytest.param(
            dict(
                a=pd.Categorical(["a", "b", "c"]),
                b=pd.Categorical(["d", "d", "f"]),
                c=pd.Categorical(["g", "h", "h"]),
            ),
            ["a", "c"],
            pd.Categorical(["a_g", "b_h", "c_h"]),
            id="two_of_three-1",
        ),
        pytest.param(
            dict(
                a=pd.Categorical(["a", "b", "c"]),
                b=pd.Categorical(["d", "d", "f"]),
                c=pd.Categorical(["g", "h", "h"]),
            ),
            ["b", "c"],
            pd.Categorical(["d_g", "d_h", "f_h"]),
            id="two_of_three-2",
        ),
    ],
)
def test_combine_categories(label_cols, cols, expected):
    from scanpy.get._aggregated import _combine_categories

    label_df = pd.DataFrame(label_cols)
    result, result_label_df = _combine_categories(label_df, cols)

    assert isinstance(result, pd.Categorical)

    pd.testing.assert_extension_array_equal(result, expected)

    pd.testing.assert_index_equal(
        pd.Index(result), result_label_df.index.astype("category")
    )

    reconstructed_df = pd.DataFrame(
        [x.split("_") for x in result], columns=cols, index=result.astype(str)
    ).astype("category")
    pd.testing.assert_frame_equal(reconstructed_df, result_label_df)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_aggregate_arraytype(array_type, metric):
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    aggregate = sc.get.aggregate(adata, ["louvain"], metric)
    assert isinstance(aggregate.layers[metric], np.ndarray)


def test_aggregate_obsm_varm():
    adata_obsm = sc.datasets.blobs()
    adata_obsm.obs["blobs"] = adata_obsm.obs["blobs"].astype(str)
    adata_obsm.obsm["test"] = adata_obsm.X[:, ::2].copy()
    adata_varm = adata_obsm.T.copy()

    result_obsm = sc.get.aggregate(adata_obsm, "blobs", ["sum", "mean"], obsm="test")
    result_varm = sc.get.aggregate(adata_varm, "blobs", ["sum", "mean"], varm="test")

    assert_equal(result_obsm, result_varm.T)

    expected_sum = (
        pd.DataFrame(adata_obsm.obsm["test"], index=adata_obsm.obs_names)
        .groupby(adata_obsm.obs["blobs"], observed=True)
        .sum()
    )
    expected_mean = (
        pd.DataFrame(adata_obsm.obsm["test"], index=adata_obsm.obs_names)
        .groupby(adata_obsm.obs["blobs"], observed=True)
        .mean()
    )

    assert_equal(expected_sum.values, result_obsm.layers["sum"])
    assert_equal(expected_mean.values, result_obsm.layers["mean"])


def test_aggregate_obsm_labels():
    from itertools import chain, repeat

    label_counts = [("a", 5), ("b", 3), ("c", 4)]
    blocks = [np.ones((n, 1)) for _, n in label_counts]
    obs_names = pd.Index(
        [f"cell_{i:02d}" for i in range(sum(b.shape[0] for b in blocks))]
    )
    entry = pd.DataFrame(
        sparse.block_diag(blocks).toarray(),
        columns=[f"dim_{i}" for i in range(len(label_counts))],
        index=obs_names,
    )

    adata = ad.AnnData(
        obs=pd.DataFrame(
            {
                "labels": list(
                    chain.from_iterable(repeat(l, n) for (l, n) in label_counts)
                )
            },
            index=obs_names,
        ),
        var=pd.DataFrame(index=["gene_0"]),
        obsm={"entry": entry},
    )

    expected = ad.AnnData(
        obs=pd.DataFrame({"labels": pd.Categorical(list("abc"))}, index=list("abc")),
        var=pd.DataFrame(index=[f"dim_{i}" for i in range(3)]),
        layers={
            "sum": np.diag([n for _, n in label_counts]),
        },
    )
    result = sc.get.aggregate(adata, by="labels", func="sum", obsm="entry")
    assert_equal(expected, result)


def test_dispatch_not_implemented():
    adata = sc.datasets.blobs()
    with pytest.raises(NotImplementedError):
        sc.get.aggregate(adata.X, adata.obs["blobs"], "sum")


def test_factors():
    from itertools import product

    obs = pd.DataFrame(product(range(5), repeat=4), columns=list("abcd"))
    obs.index = [f"cell_{i:04d}" for i in range(obs.shape[0])]
    adata = ad.AnnData(
        X=np.arange(obs.shape[0]).reshape(-1, 1),
        obs=obs,
    )

    res = sc.get.aggregate(adata, by=["a", "b", "c", "d"], func="sum")
    np.testing.assert_equal(res.layers["sum"], adata.X)
