import scanpy.api as sc
import scanpy
import numpy as np
from scipy import sparse

from scanpy.preprocessing.qc import top_proportions, top_segment_proportions

def test_proportions():
    a_dense = np.ones((100, 100))
    a_sparse = sparse.csr_matrix(a_dense)
    prop_dense = top_proportions(a_dense, 100)
    assert (prop_dense[:, -1] == 1).all()
    assert np.array_equal(np.sort(prop_dense, axis=1), prop_dense)
    assert np.apply_along_axis(
        lambda x: len(np.unique(x)) == 1, 0, prop_dense
    ).all()
    assert (prop_dense[:, 49] == .5).all()

    prop_sparse = top_proportions(a_sparse, 100)
    assert (prop_sparse[:, -1] == 1).all()
    assert np.array_equal(np.sort(prop_sparse, axis=1), prop_sparse)
    assert np.apply_along_axis(
        lambda x: len(np.unique(x)) == 1, 0, prop_sparse
    ).all()
    assert (prop_sparse[:, 49] == .5).all()


def test_top_segments_dense():
    a = np.ones((300, 100))
    seg = top_segment_proportions(a, [50, 100])
    assert (seg[:, 0] == .5).all()
    assert (seg[:, 1] == 1.).all()
    segfull = top_segment_proportions(a, np.arange(100)+1)
    propfull = top_proportions(a, 100)
    assert (segfull == propfull).all()


def test_segments_binary():
    a = np.concatenate([np.zeros((300, 50)), np.ones((300, 50))], 1)
    a = np.apply_along_axis(np.random.permutation, 1, a)
    seg = top_segment_proportions(a, [25, 50, 100])
    assert (seg[:, 0] == .5).all()
    assert (top_segment_proportions(a, [25]) == .5).all()
    assert (seg[:, 1] == 1.).all()
    assert (seg[:, 2] == 1.).all()
    segfull = top_segment_proportions(a, np.arange(100)+1)
    propfull = top_proportions(a, 100)
    assert (segfull == propfull).all()


def test_top_segments_sparse():
    a_dense = np.ones((300, 100))
    for fmt in [sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix]:
        print(fmt)
        a = fmt(a_dense)
        seg = top_segment_proportions(a, [50, 100])
        assert (seg[:, 0] == .5).all()
        assert (seg[:, 1] == 1.).all()
        segfull = top_segment_proportions(a, np.arange(100)+1)
        propfull = top_proportions(a, 100)
        assert (segfull == propfull).all()


def test_qc_metrics():
    adata = sc.AnnData(X=sparse.csr_matrix(
        np.random.binomial(100, .005, (1000, 1000))))
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    assert (adata.obs["total_features_by_counts"] > 0).all()
    assert (adata.obs["total_features_by_counts"] < adata.shape[1]).all()
    assert (adata.obs["total_counts"] == np.ravel(adata.X.sum(axis=1))).all()
    for col in adata.obs.columns:
        if col.startswith("pct_counts_in_top"):
            assert (adata.obs[col] <= 100).all()
            assert (adata.obs[col] >= 0).all()
    for col in adata.var.columns:
        assert (adata.var[col] >= 0).all()
    assert (adata.var["mean_counts"] < np.ravel(adata.X.max(axis=0).todense())).all()

def test_qc_metrics_format():
    a = np.random.binomial(100, .005, (1000, 1000))
    adata_dense = sc.AnnData(X=a)
    sc.pp.calculate_qc_metrics(adata_dense, inplace=True)
    for fmt in [sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix]:
        adata = sc.AnnData(X=fmt(a))
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        assert np.allclose(adata.obs, adata_dense.obs)
        assert np.allclose(adata.var, adata_dense.var)
