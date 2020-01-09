import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from scanpy.preprocessing._qc import (
    top_proportions,
    top_segment_proportions,
    describe_var,
    describe_obs,
)


@pytest.fixture
def anndata():
    a = np.random.binomial(100, 0.005, (1000, 1000))
    adata = AnnData(
        sparse.csr_matrix(a),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(a.shape[0])]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(a.shape[1])]),
    )
    return adata


@pytest.mark.parametrize(
    "a",
    [np.ones((100, 100)), sparse.csr_matrix(np.ones((100, 100)))],
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
    "cls", [np.asarray, sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix]
)
def test_top_segments(cls):
    a = cls(np.ones((300, 100)))
    seg = top_segment_proportions(a, [50, 100])
    assert (seg[:, 0] == 0.5).all()
    assert (seg[:, 1] == 1.0).all()
    segfull = top_segment_proportions(a, np.arange(100) + 1)
    propfull = top_proportions(a, 100)
    assert (segfull == propfull).all()


# While many of these are trivial,
# they’re also just making sure the metrics are there
def test_qc_metrics():
    adata = AnnData(X=sparse.csr_matrix(np.random.binomial(100, 0.005, (1000, 1000))))
    adata.var["mito"] = np.concatenate(
        (np.ones(100, dtype=bool), np.zeros(900, dtype=bool))
    )
    adata.var["negative"] = False
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "negative"], inplace=True)
    assert (adata.obs["n_genes_by_counts"] < adata.shape[1]).all()
    assert (
        adata.obs["n_genes_by_counts"] >= adata.obs["log1p_n_genes_by_counts"]
    ).all()
    assert (adata.obs["total_counts"] == np.ravel(adata.X.sum(axis=1))).all()
    assert (adata.obs["total_counts"] >= adata.obs["log1p_total_counts"]).all()
    assert (
        adata.obs["total_counts_mito"] >= adata.obs["log1p_total_counts_mito"]
    ).all()
    assert (adata.obs["total_counts_negative"] == 0).all()
    assert (
        adata.obs["pct_counts_in_top_50_genes"]
        <= adata.obs["pct_counts_in_top_100_genes"]
    ).all()
    for col in filter(lambda x: "negative" not in x, adata.obs.columns):
        assert (adata.obs[col] >= 0).all()  # Values should be positive or zero
        assert (adata.obs[col] != 0).any().all()  # Nothing should be all zeros
        if col.startswith("pct_counts_in_top"):
            assert (adata.obs[col] <= 100).all()
            assert (adata.obs[col] >= 0).all()
    for col in adata.var.columns:
        assert (adata.var[col] >= 0).all()
    assert (adata.var["mean_counts"] < np.ravel(adata.X.max(axis=0).todense())).all()
    assert (adata.var["mean_counts"] >= adata.var["log1p_mean_counts"]).all()
    assert (adata.var["total_counts"] >= adata.var["log1p_total_counts"]).all()
    # Should return the same thing if run again
    old_obs, old_var = adata.obs.copy(), adata.var.copy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "negative"], inplace=True)
    assert set(adata.obs.columns) == set(old_obs.columns)
    assert set(adata.var.columns) == set(old_var.columns)
    for col in adata.obs:
        assert np.allclose(adata.obs[col], old_obs[col])
    for col in adata.var:
        assert np.allclose(adata.var[col], old_var[col])


def adata_mito():
    a = np.random.binomial(100, 0.005, (1000, 1000))
    init_var = pd.DataFrame(
        dict(mito=np.concatenate((np.ones(100, dtype=bool), np.zeros(900, dtype=bool))))
    )
    adata_dense = AnnData(X=a, var=init_var.copy())
    return adata_dense, init_var


@pytest.mark.parametrize(
    "cls", [np.asarray, sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix]
)
def test_qc_metrics_format(cls):
    adata_dense, init_var = adata_mito()
    sc.pp.calculate_qc_metrics(adata_dense, qc_vars=["mito"], inplace=True)
    adata = AnnData(X=cls(adata_dense.X), var=init_var.copy())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
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


def test_layer_raw(anndata):
    adata = anndata.copy()
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


def test_inner_methods(anndata):
    adata = anndata.copy()
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
