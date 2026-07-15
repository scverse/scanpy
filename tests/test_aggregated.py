from __future__ import annotations

from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata.tests import helpers
from scipy import sparse

import scanpy as sc
from scanpy._compat import DaskArray
from scanpy._utils import _resolve_axis, get_literal_vals
from scanpy.get._aggregated import AggType
from testing.scanpy._helpers import assert_equal
from testing.scanpy._helpers.data import pbmc3k_processed
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES as ARRAY_TYPES_ALL
from testing.scanpy._pytest.params import ARRAY_TYPES_MEM, param_with

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from numpy.typing import NDArray

    from scanpy._compat import CSRBase

VALID_ARRAY_TYPES = [
    param_with(
        at,
        marks=[
            pytest.mark.xfail(reason="aggregate not implemented for array-api arrays")
        ],
    )
    if at.id == "jax_array"
    else at
    for at in ARRAY_TYPES_ALL
    if at.id
    not in {
        "dask_array_dense",
        "dask_array_sparse",
    }
]


@pytest.fixture(params=get_literal_vals(AggType))
def metric(request: pytest.FixtureRequest) -> AggType:
    return request.param


def xfail_dask_median(
    adata: ad.AnnData,
    metric: AggType,
    request: pytest.FixtureRequest,
) -> None:
    if isinstance(adata.X, DaskArray) and metric == "median":
        reason = "Median calculation not implemented for Dask"
        request.applymarker(pytest.mark.xfail(reason=reason))


@pytest.mark.parametrize("axis", [0, 1])
def test_mask(axis: Literal[0, 1]) -> None:
    blobs = sc.datasets.blobs()
    mask = blobs.obs["blobs"] == 0
    blobs.obs["mask_col"] = mask
    if axis == 1:
        blobs = blobs.T
    by_name = sc.get.aggregate(blobs, "blobs", "sum", axis=axis, mask="mask_col")
    by_value = sc.get.aggregate(blobs, "blobs", "sum", axis=axis, mask=mask)

    assert_equal(by_name, by_value)

    assert np.all(by_name["0"].layers["sum"] == 0)


@pytest.mark.parametrize("array_type", VALID_ARRAY_TYPES)
def test_aggregate_vs_pandas(
    metric: AggType, array_type, request: pytest.FixtureRequest
) -> None:
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    xfail_dask_median(adata, metric, request)
    adata.obs["percent_mito_binned"] = pd.cut(adata.obs["percent_mito"], bins=5)
    result = sc.get.aggregate(adata, ["louvain", "percent_mito_binned"], metric)
    if isinstance(adata.X, DaskArray):
        adata.X = adata.X.compute()
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
            adata
            .to_df()
            .astype(np.float64)
            .join(adata.obs[["louvain", "percent_mito_binned"]])
            .groupby(["louvain", "percent_mito_binned"], observed=True)
            .agg(metric)
        )
    expected.index = expected.index.to_frame().astype("string").agg("_".join, axis=1)
    expected.index.name = None
    expected.columns.name = None
    if isinstance(result.layers[metric], DaskArray):
        result.layers[metric] = result.layers[metric].compute()
    result_df = result.to_df(layer=metric)
    result_df.index.name = None
    result_df.columns.name = None

    pd.testing.assert_frame_equal(result_df, expected, check_dtype=False, atol=1e-5)


@pytest.mark.parametrize("array_type", VALID_ARRAY_TYPES)
def test_aggregate_axis(
    array_type, metric: AggType, request: pytest.FixtureRequest
) -> None:
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    xfail_dask_median(adata, metric, request)
    expected = sc.get.aggregate(adata, ["louvain"], metric)
    actual = sc.get.aggregate(adata.T, ["louvain"], metric, axis=1)
    actual = actual.T
    assert_equal(expected, actual)


def test_aggregate_entry() -> None:
    args = ("blobs", ["mean", "var", "count_nonzero"])

    adata = sc.datasets.blobs()
    x_result = sc.get.aggregate(adata, *args)
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

    x_result_min = x_result.copy()
    del x_result_min.var
    x_result_min.var_names = [str(x) for x in np.arange(x_result_min.n_vars)]

    assert_equal(x_result, layer_result)
    assert_equal(x_result_min, obsm_result)
    assert_equal(x_result.layers, obsm_result.layers)
    assert_equal(x_result.layers, varm_result.T.layers)


def test_aggregate_incorrect_dim() -> None:
    adata = pbmc3k_processed().raw.to_adata()

    with pytest.raises(ValueError, match="was 'foo'"):
        sc.get.aggregate(adata, ["louvain"], "sum", axis="foo")


def to_bad_chunking(x: CSRBase) -> DaskArray:
    import dask.array as da

    return da.from_array(
        x,
        chunks=(x.shape[0] // 2, x.shape[1] // 2),
        meta=sparse.csr_matrix(np.array([])),  # noqa: TID251
    )


def to_csc(x: CSRBase) -> DaskArray:
    import dask.array as da

    return da.from_array(
        x.tocsc(),
        chunks=(x.shape[0] // 2, x.shape[1]),
        meta=sparse.csc_matrix(np.array([])),  # noqa: TID251
    )


@needs.dask
@pytest.mark.parametrize(
    ("func", "error_msg"),
    [
        pytest.param(
            to_bad_chunking, r"Feature axis must be unchunked", id="bad_chunking"
        ),
    ],
)
def test_aggregate_bad_dask_array(
    func: Callable[[CSRBase], DaskArray], error_msg: str
) -> None:
    adata = pbmc3k_processed().raw.to_adata()
    adata.X = func(adata.X)
    with pytest.raises(ValueError, match=error_msg):
        sc.get.aggregate(adata, ["louvain"], "sum")


@pytest.mark.parametrize("axis_name", ["obs", "var"])
def test_aggregate_axis_specification(axis_name: Literal["obs", "var"]) -> None:
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
            np.block([
                [np.ones((2, 2)), np.zeros((2, 2))],
                [np.zeros((2, 2)), np.ones((2, 2))],
            ]),
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
                    {
                        "a": pd.Categorical(["a", "a", "b"]),
                        "b": pd.Categorical(["c", "d", "d"]),
                        "n_obs_aggregated": [1, 1, 2],
                    },
                    index=["a_c", "a_d", "b_d"],
                ),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "count_nonzero": np.array([
                        [1, 1, 0, 0],
                        [1, 1, 0, 0],
                        [0, 0, 2, 2],
                    ]),
                    # "sum": np.array([[2, 0], [0, 2]]),
                    # "mean": np.array([[1, 0], [0, 1]]),
                },
            ),
            id="count_nonzero",
        ),
        pytest.param(
            np.block([
                [np.ones((2, 2)), np.zeros((2, 2))],
                [np.zeros((2, 2)), np.ones((2, 2))],
            ]),
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
                    {
                        "a": pd.Categorical(["a", "a", "b"]),
                        "b": pd.Categorical(["c", "d", "d"]),
                        "n_obs_aggregated": [1, 1, 2],
                    },
                    index=["a_c", "a_d", "b_d"],
                ),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "sum": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2]]),
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                    "count_nonzero": np.array([
                        [1, 1, 0, 0],
                        [1, 1, 0, 0],
                        [0, 0, 2, 2],
                    ]),
                },
            ),
            id="sum-mean-count_nonzero",
        ),
        pytest.param(
            np.block([
                [np.ones((2, 2)), np.zeros((2, 2))],
                [np.zeros((2, 2)), np.ones((2, 2))],
            ]),
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
                    {
                        "a": pd.Categorical(["a", "a", "b"]),
                        "b": pd.Categorical(["c", "d", "d"]),
                        "n_obs_aggregated": [1, 1, 2],
                    },
                    index=["a_c", "a_d", "b_d"],
                ),
                var=pd.DataFrame(index=[f"gene_{i}" for i in range(4)]),
                layers={
                    "mean": np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]]),
                },
            ),
            id="mean",
        ),
    ],
)
def test_aggregate_examples(
    matrix: NDArray[np.int64],
    df: pd.DataFrame,
    keys: list[str],
    metrics: list[AggType],
    expected: ad.AnnData,
) -> None:
    adata = ad.AnnData(
        X=matrix,
        obs=df,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(matrix.shape[1])]),
    )
    result = sc.get.aggregate(adata, by=keys, func=metrics)

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
def test_combine_categories(
    label_cols: dict[str, pd.Categorical], cols: list[str], expected: pd.Categorical
) -> None:
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


@pytest.mark.parametrize("array_type", VALID_ARRAY_TYPES)
def test_aggregate_arraytype(
    array_type, metric: AggType, request: pytest.FixtureRequest
) -> None:
    adata = pbmc3k_processed().raw.to_adata()
    adata = adata[
        adata.obs["louvain"].isin(adata.obs["louvain"].cat.categories[:5]), :1_000
    ].copy()
    adata.X = array_type(adata.X)
    xfail_dask_median(adata, metric, request)
    aggregate = sc.get.aggregate(adata, ["louvain"], metric)
    assert isinstance(
        aggregate.layers[metric],
        DaskArray if isinstance(adata.X, DaskArray) else np.ndarray,
    )


def test_aggregate_obsm_varm() -> None:
    adata_obsm = sc.datasets.blobs()
    adata_obsm.obs["blobs"] = adata_obsm.obs["blobs"].astype(str)
    adata_obsm.obsm["test"] = adata_obsm.X[:, ::2].copy()
    adata_varm = adata_obsm.T.copy()

    result_obsm = sc.get.aggregate(adata_obsm, "blobs", ["sum", "mean"], obsm="test")
    result_varm = sc.get.aggregate(adata_varm, "blobs", ["sum", "mean"], varm="test")

    assert_equal(result_obsm, result_varm.T)

    expected_sum = (
        pd
        .DataFrame(adata_obsm.obsm["test"], index=adata_obsm.obs_names)
        .groupby(adata_obsm.obs["blobs"], observed=True)
        .sum()
    )
    expected_mean = (
        pd
        .DataFrame(adata_obsm.obsm["test"], index=adata_obsm.obs_names)
        .groupby(adata_obsm.obs["blobs"], observed=True)
        .mean()
    )

    assert_equal(expected_sum.values, result_obsm.layers["sum"])
    assert_equal(expected_mean.values, result_obsm.layers["mean"])


def test_aggregate_obsm_labels() -> None:
    from itertools import chain, repeat

    label_counts = [("a", 5), ("b", 3), ("c", 4)]
    blocks = [sparse.csr_array(np.ones((n, 1))) for _, n in label_counts]  # noqa: TID251
    obs_names = pd.Index([
        f"cell_{i:02d}" for i in range(sum(b.shape[0] for b in blocks))
    ])
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
        obs=pd.DataFrame(
            {
                "labels": pd.Categorical([lc[0] for lc in label_counts]),
                "n_obs_aggregated": [lc[1] for lc in label_counts],
            },
            index=[lc[0] for lc in label_counts],
        ),
        var=pd.DataFrame(index=[f"dim_{i}" for i in range(3)]),
        layers={
            "sum": np.diag([n for _, n in label_counts]),
        },
    )
    result = sc.get.aggregate(adata, by="labels", func="sum", obsm="entry")
    assert_equal(expected, result)


def test_dispatch_not_implemented() -> None:
    adata = sc.datasets.blobs()
    with pytest.raises(NotImplementedError):
        sc.get.aggregate(adata.X, adata.obs["blobs"], "sum")


def test_factors() -> None:
    from itertools import product

    obs = pd.DataFrame(product(range(5), repeat=4), columns=list("abcd"))
    obs.index = [f"cell_{i:04d}" for i in range(obs.shape[0])]
    adata = ad.AnnData(
        X=np.arange(obs.shape[0]).reshape(-1, 1),
        obs=obs,
    )

    res = sc.get.aggregate(adata, by=["a", "b", "c", "d"], func="sum")
    np.testing.assert_equal(res.layers["sum"], adata.X)


def test_nan() -> None:
    x = np.arange(6) + np.arange(6)[:, None]
    obs = pd.DataFrame(
        dict(
            cell_type=[np.nan, np.nan, "B", "C", "B", "B"],
            sample_id=["s1", "s1", "s1", "s2", "s2", "s2"],
            patient_type=[*(["responder"] * 3), *(["control"] * 3)],
        ),
        index=[f"cell{i}" for i in range(x.shape[0])],
    )
    adata = ad.AnnData(x, obs=obs)

    adata_agg = sc.get.aggregate(
        adata, by=["sample_id", "patient_type", "cell_type"], func="sum", layer=None
    )

    assert adata_agg.obs.index.tolist() == [
        "s1_responder_B",
        "s2_control_B",
        "s2_control_C",
    ]
    assert adata_agg.obs["n_obs_aggregated"].tolist() == [1, 2, 1]


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_var_no_catastrophic_cancellation(
    request: pytest.FixtureRequest, array_type
) -> None:
    # Values of the form `offset + tiny_noise` make the textbook two-pass
    # formula sum(x**2)/n - (sum(x)/n)**2 lose ~all precision: both terms are
    # ~n*offset**2 ≈ 1e19 in float64 (precision ~1e3) but their difference is
    # the variance ~1e-3, far below the rounding noise. Welford's online
    # algorithm avoids the subtraction entirely.
    if array_type is helpers.as_dense_jax_array:
        request.applymarker(
            pytest.mark.xfail(reason="aggregate not implemented for jax arrays")
        )
    n_per_group, n_features = 1000, 4
    offset, std = 1e8, 1e-3
    groups = ["a", "b"]
    x = np.vstack([
        offset
        + std * np.random.default_rng().standard_normal((n_per_group, n_features))
        for _ in groups
    ])
    obs = pd.DataFrame(
        {"group": pd.Categorical(np.repeat(groups, n_per_group))},
        index=[f"cell_{i}" for i in range(x.shape[0])],
    )
    adata = ad.AnnData(X=array_type(x), obs=obs)

    expected = np.vstack([
        np.var(x[i * n_per_group : (i + 1) * n_per_group], axis=0, ddof=0)
        for i in range(len(groups))
    ])
    # Sanity: textbook formula on this data is either catastrophically wrong by a large magnitude relative to the expected
    # or the sum-sq and sq-sum in naive are literally identical due to precision errors at the upper bound of the range.
    naive = np.vstack([
        (xg**2).mean(axis=0) - xg.mean(axis=0) ** 2
        for xg in (
            x[i * n_per_group : (i + 1) * n_per_group] for i in range(len(groups))
        )
    ])
    diff_magnitude = np.abs(naive - expected) / expected
    all_large = (diff_magnitude > 1e5).all()
    if not all_large:
        does_naive_fully_cancel = naive == 0
        assert does_naive_fully_cancel.any()
        assert (diff_magnitude[does_naive_fully_cancel] == 1).all()

    result = sc.get.aggregate(adata, by="group", func="var", dof=0).layers["var"]
    if isinstance(result, DaskArray):
        result = result.compute()
    np.testing.assert_allclose(result, expected, rtol=1e-4)
