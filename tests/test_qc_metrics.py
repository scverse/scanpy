from __future__ import annotations

from importlib.metadata import version

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from anndata.tests.helpers import assert_equal
from packaging.version import Version
from scipy import sparse

import scanpy as sc
from scanpy._compat import DaskArray
from scanpy._utils import axis_sum
from scanpy.preprocessing._qc import (
    describe_obs,
    describe_var,
    top_proportions,
    top_segment_proportions,
)
from testing.scanpy._helpers import as_sparse_dask_array, maybe_dask_process_context
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES, ARRAY_TYPES_MEM


@pytest.fixture
def adata() -> AnnData:
    a = np.random.binomial(100, 0.005, (1000, 1000))
    adata = AnnData(
        sparse.csr_matrix(a),  # noqa: TID251
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(a.shape[0])]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(a.shape[1])]),
    )
    return adata


def prepare_adata(adata: AnnData) -> AnnData:
    if isinstance(adata.X, DaskArray):
        adata.X = adata.X.rechunk((10, -1))
    adata.var["mito"] = np.concatenate(
        (np.ones(100, dtype=bool), np.zeros(900, dtype=bool))
    )
    adata.var["negative"] = False
    return adata


@pytest.fixture(params=ARRAY_TYPES)
def adata_prepared(request: pytest.FixtureRequest, adata: AnnData) -> AnnData:
    adata.X = request.param(adata.X)
    return prepare_adata(adata)


@pytest.mark.parametrize(
    "a",
    [np.ones((100, 100)), sparse.csr_matrix(np.ones((100, 100)))],  # noqa: TID251
    ids=["dense", "sparse"],
)
def test_proportions(a):
    prop = top_proportions(a, 100)
    assert (prop[:, -1] == 1).all()
    assert np.array_equal(np.sort(prop, axis=1), prop)
    assert np.apply_along_axis(lambda x: len(np.unique(x)) == 1, 0, prop).all()
    assert (prop[:, 49] == 0.5).all()


def test_segments_binary():
    a = np.concatenate([np.zeros((300, 50)), np.ones((300, 50))], 1)
    a = np.apply_along_axis(np.random.permutation, 1, a)
    seg = top_segment_proportions(a, [25, 50, 100])
    assert (seg[:, 0] == 0.5).all()
    assert (top_segment_proportions(a, [25]) == 0.5).all()
    assert (seg[:, 1] == 1.0).all()
    assert (seg[:, 2] == 1.0).all()
    segfull = top_segment_proportions(a, np.arange(100) + 1)
    propfull = top_proportions(a, 100)
    assert (segfull == propfull).all()


@pytest.mark.parametrize(
    "array_type", [*ARRAY_TYPES, pytest.param(sparse.coo_matrix, id="scipy_coo")]
)
def test_top_segments(request: pytest.FixtureRequest, array_type):
    if "dask" in array_type.__name__:
        reason = "DaskArray not yet supported"
        request.applymarker(pytest.mark.xfail(reason=reason))
    a = array_type(np.ones((300, 100)))
    seg = top_segment_proportions(a, [50, 100])
    assert (seg[:, 0] == 0.5).all()
    assert (seg[:, 1] == 1.0).all()
    segfull = top_segment_proportions(a, np.arange(100) + 1)
    propfull = top_proportions(a, 100)
    assert (segfull == propfull).all()


# While many of these are trivial,
# they’re also just making sure the metrics are there
def test_qc_metrics(adata_prepared: AnnData):
    with maybe_dask_process_context():
        sc.pp.calculate_qc_metrics(
            adata_prepared, qc_vars=["mito", "negative"], inplace=True
        )
    X = (
        adata_prepared.X.compute()
        if isinstance(adata_prepared.X, DaskArray)
        else adata_prepared.X
    )
    max_X = X.max(axis=0)
    if isinstance(max_X, sparse.coo_matrix):
        max_X = max_X.toarray()
    elif isinstance(max_X, DaskArray):
        max_X = max_X.compute()
    assert (adata_prepared.obs["n_genes_by_counts"] < adata_prepared.shape[1]).all()
    assert (
        adata_prepared.obs["n_genes_by_counts"]
        >= adata_prepared.obs["log1p_n_genes_by_counts"]
    ).all()
    assert (
        adata_prepared.obs["total_counts"]
        == np.ravel(axis_sum(adata_prepared.X, axis=1))
    ).all()
    assert (
        adata_prepared.obs["total_counts"] >= adata_prepared.obs["log1p_total_counts"]
    ).all()
    assert (
        adata_prepared.obs["total_counts_mito"]
        >= adata_prepared.obs["log1p_total_counts_mito"]
    ).all()
    assert (adata_prepared.obs["total_counts_negative"] == 0).all()
    assert (
        adata_prepared.obs["pct_counts_in_top_50_genes"]
        <= adata_prepared.obs["pct_counts_in_top_100_genes"]
    ).all()
    for col in filter(lambda x: "negative" not in x, adata_prepared.obs.columns):
        assert (adata_prepared.obs[col] >= 0).all()  # Values should be positive or zero
        assert (adata_prepared.obs[col] != 0).any().all()  # Nothing should be all zeros
        if col.startswith("pct_counts_in_top"):
            assert (adata_prepared.obs[col] <= 100).all()
            assert (adata_prepared.obs[col] >= 0).all()
    for col in adata_prepared.var.columns:
        assert (adata_prepared.var[col] >= 0).all()
    assert (adata_prepared.var["mean_counts"] < np.ravel(max_X)).all()
    assert (
        adata_prepared.var["mean_counts"] >= adata_prepared.var["log1p_mean_counts"]
    ).all()
    assert (
        adata_prepared.var["total_counts"] >= adata_prepared.var["log1p_total_counts"]
    ).all()


def test_qc_metrics_idempotent(adata_prepared: AnnData):
    with maybe_dask_process_context():
        sc.pp.calculate_qc_metrics(
            adata_prepared, qc_vars=["mito", "negative"], inplace=True
        )
        old_obs, old_var = adata_prepared.obs.copy(), adata_prepared.var.copy()
        sc.pp.calculate_qc_metrics(
            adata_prepared, qc_vars=["mito", "negative"], inplace=True
        )
    assert set(adata_prepared.obs.columns) == set(old_obs.columns)
    assert set(adata_prepared.var.columns) == set(old_var.columns)
    for col in adata_prepared.obs:
        assert np.allclose(adata_prepared.obs[col], old_obs[col])
    for col in adata_prepared.var:
        assert np.allclose(adata_prepared.var[col], old_var[col])


def test_qc_metrics_no_log1p(adata_prepared: AnnData):
    with maybe_dask_process_context():
        sc.pp.calculate_qc_metrics(
            adata_prepared, qc_vars=["mito", "negative"], log1p=False, inplace=True
        )
    assert not np.any(adata_prepared.obs.columns.str.startswith("log1p_"))
    assert not np.any(adata_prepared.var.columns.str.startswith("log1p_"))


@needs.dask
@pytest.mark.anndata_dask_support
@pytest.mark.parametrize("log1p", [True, False], ids=["log1p", "no_log1p"])
def test_dask_against_in_memory(adata, log1p):
    adata_as_dask = adata.copy()
    adata_as_dask.X = as_sparse_dask_array(adata.X)
    adata = prepare_adata(adata)
    adata_as_dask = prepare_adata(adata_as_dask)
    with maybe_dask_process_context():
        sc.pp.calculate_qc_metrics(
            adata_as_dask, qc_vars=["mito", "negative"], log1p=log1p, inplace=True
        )
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito", "negative"], log1p=log1p, inplace=True
    )
    assert_equal(adata, adata_as_dask)


def adata_mito():
    a = np.random.binomial(100, 0.005, (1000, 1000))
    init_var = pd.DataFrame(
        dict(mito=np.concatenate((np.ones(100, dtype=bool), np.zeros(900, dtype=bool))))
    )
    adata_dense = AnnData(X=a, var=init_var.copy())
    return adata_dense, init_var


skip_if_adata_0_12 = pytest.mark.skipif(
    Version(version("anndata")) >= Version("0.12.0.dev0"),
    reason="Newer AnnData removes implicit support for COO matrices",
)


@pytest.mark.parametrize(
    "cls",
    [
        *ARRAY_TYPES_MEM,
        pytest.param(sparse.coo_matrix, marks=[skip_if_adata_0_12], id="scipy_coo"),
    ],
)
def test_qc_metrics_format(cls):
    adata_dense, init_var = adata_mito()
    sc.pp.calculate_qc_metrics(adata_dense, qc_vars=["mito"], inplace=True)
    adata = AnnData(X=cls(adata_dense.X), var=init_var.copy())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    assert np.allclose(adata.obs, adata_dense.obs)
    for col in adata.var:  # np.allclose doesn't like mix of types
        assert np.allclose(adata.var[col], adata_dense.var[col])


def test_qc_metrics_format_str_qc_vars():
    adata_dense, init_var = adata_mito()
    sc.pp.calculate_qc_metrics(adata_dense, qc_vars="mito", inplace=True)
    adata = AnnData(X=adata_dense.X, var=init_var.copy())
    sc.pp.calculate_qc_metrics(adata, qc_vars="mito", inplace=True)
    assert np.allclose(adata.obs, adata_dense.obs)
    for col in adata.var:  # np.allclose doesn't like mix of types
        assert np.allclose(adata.var[col], adata_dense.var[col])


def test_qc_metrics_percentage():  # In response to #421
    adata_dense, init_var = adata_mito()
    sc.pp.calculate_qc_metrics(adata_dense, percent_top=[])
    sc.pp.calculate_qc_metrics(adata_dense, percent_top=())
    sc.pp.calculate_qc_metrics(adata_dense, percent_top=None)
    sc.pp.calculate_qc_metrics(adata_dense, percent_top=[1, 2, 3, 10])
    sc.pp.calculate_qc_metrics(adata_dense, percent_top=[1])
    with pytest.raises(IndexError):
        sc.pp.calculate_qc_metrics(adata_dense, percent_top=[1, 2, 3, -5])
    with pytest.raises(IndexError):
        sc.pp.calculate_qc_metrics(adata_dense, percent_top=[20, 30, 1001])


def test_layer_raw(adata: AnnData):
    adata = adata.copy()
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    obs_orig, var_orig = sc.pp.calculate_qc_metrics(adata)
    sc.pp.log1p(adata)  # To be sure they aren't reusing it
    obs_layer, var_layer = sc.pp.calculate_qc_metrics(adata, layer="counts")
    obs_raw, var_raw = sc.pp.calculate_qc_metrics(adata, use_raw=True)
    assert np.allclose(obs_orig, obs_layer)
    assert np.allclose(obs_orig, obs_raw)
    assert np.allclose(var_orig, var_layer)
    assert np.allclose(var_orig, var_raw)


def test_inner_methods(adata: AnnData):
    adata = adata.copy()
    full_inplace = adata.copy()
    partial_inplace = adata.copy()
    obs_orig, var_orig = sc.pp.calculate_qc_metrics(adata)
    assert np.all(obs_orig == describe_obs(adata))
    assert np.all(var_orig == describe_var(adata))
    sc.pp.calculate_qc_metrics(full_inplace, inplace=True)
    describe_obs(partial_inplace, inplace=True)
    describe_var(partial_inplace, inplace=True)
    assert np.all(full_inplace.obs == partial_inplace.obs)
    assert np.all(full_inplace.var == partial_inplace.var)
    assert np.all(partial_inplace.obs[obs_orig.columns] == obs_orig)
    assert np.all(partial_inplace.var[var_orig.columns] == var_orig)
