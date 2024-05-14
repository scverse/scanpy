from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from anndata.tests.helpers import asarray, assert_equal
from numpy.testing import assert_allclose
from scipy import sparse as sp
from scipy.sparse import issparse

import scanpy as sc
from testing.scanpy._helpers import (
    anndata_v0_8_constructor_compat,
    check_rep_mutation,
    check_rep_results,
)
from testing.scanpy._helpers.data import pbmc3k, pbmc68k_reduced
from testing.scanpy._pytest.params import ARRAY_TYPES


def test_log1p(tmp_path):
    A = np.random.rand(200, 10).astype(np.float32)
    A_l = np.log1p(A)
    ad = AnnData(A.copy())
    ad2 = AnnData(A.copy())
    ad3 = AnnData(A.copy())
    ad3.filename = tmp_path / "test.h5ad"
    sc.pp.log1p(ad)
    assert np.allclose(ad.X, A_l)
    sc.pp.log1p(ad2, chunked=True)
    assert np.allclose(ad2.X, ad.X)
    sc.pp.log1p(ad3, chunked=True)
    assert np.allclose(ad3.X, ad.X)

    # Test base
    ad4 = AnnData(A)
    sc.pp.log1p(ad4, base=2)
    assert np.allclose(ad4.X, A_l / np.log(2))


def test_log1p_deprecated_arg():
    A = np.random.rand(200, 10).astype(np.float32)
    with pytest.warns(FutureWarning, match=r".*`X` was renamed to `data`"):
        sc.pp.log1p(X=A)


@pytest.fixture(params=[None, 2])
def base(request):
    return request.param


def test_log1p_rep(count_matrix_format, base, dtype):
    X = count_matrix_format(
        np.abs(sp.random(100, 200, density=0.3, dtype=dtype)).toarray()
    )
    check_rep_mutation(sc.pp.log1p, X, base=base)
    check_rep_results(sc.pp.log1p, X, base=base)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_mean_var(array_type):
    pbmc = pbmc3k()
    pbmc.X = array_type(pbmc.X)

    true_mean = np.mean(asarray(pbmc.X), axis=0)
    true_var = np.var(asarray(pbmc.X), axis=0, dtype=np.float64, ddof=1)

    means, variances = sc.pp._utils._get_mean_var(pbmc.X)

    np.testing.assert_allclose(true_mean, means)
    np.testing.assert_allclose(true_var, variances)


def test_mean_var_sparse():
    from sklearn.utils.sparsefuncs import mean_variance_axis

    csr64 = sp.random(10000, 1000, format="csr", dtype=np.float64)
    csc64 = csr64.tocsc()

    # Test that we're equivalent for 64 bit
    for mtx, ax in product((csr64, csc64), (0, 1)):
        scm, scv = sc.pp._utils._get_mean_var(mtx, axis=ax)
        skm, skv = mean_variance_axis(mtx, ax)
        skv *= mtx.shape[ax] / (mtx.shape[ax] - 1)

        assert np.allclose(scm, skm)
        assert np.allclose(scv, skv)

    csr32 = csr64.astype(np.float32)
    csc32 = csc64.astype(np.float32)

    # Test whether ours is more accurate for 32 bit
    for mtx32, mtx64 in [(csc32, csc64), (csr32, csr64)]:
        scm32, scv32 = sc.pp._utils._get_mean_var(mtx32)
        scm64, scv64 = sc.pp._utils._get_mean_var(mtx64)
        skm32, skv32 = mean_variance_axis(mtx32, 0)
        skm64, skv64 = mean_variance_axis(mtx64, 0)
        skv32 *= mtx.shape[0] / (mtx.shape[0] - 1)
        skv64 *= mtx.shape[0] / (mtx.shape[0] - 1)

        m_resid_sc = np.mean(np.abs(scm64 - scm32))
        m_resid_sk = np.mean(np.abs(skm64 - skm32))
        v_resid_sc = np.mean(np.abs(scv64 - scv32))
        v_resid_sk = np.mean(np.abs(skv64 - skv32))

        assert m_resid_sc < m_resid_sk
        assert v_resid_sc < v_resid_sk


def test_normalize_per_cell():
    A = np.array([[1, 0], [3, 0], [5, 6]], dtype=np.float32)
    adata = AnnData(A.copy())
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1, key_n_counts="n_counts2")
    assert adata.X.sum(axis=1).tolist() == [1.0, 1.0, 1.0]
    # now with copy option
    adata = AnnData(A.copy())
    # note that sc.pp.normalize_per_cell is also used in
    # pl.highest_expr_genes with parameter counts_per_cell_after=100
    adata_copy = sc.pp.normalize_per_cell(adata, counts_per_cell_after=1, copy=True)
    assert adata_copy.X.sum(axis=1).tolist() == [1.0, 1.0, 1.0]
    # now sparse
    adata = AnnData(A.copy())
    adata_sparse = AnnData(sp.csr_matrix(A.copy()))
    sc.pp.normalize_per_cell(adata)
    sc.pp.normalize_per_cell(adata_sparse)
    assert adata.X.sum(axis=1).tolist() == adata_sparse.X.sum(axis=1).A1.tolist()


def test_subsample():
    adata = AnnData(np.ones((200, 10)))
    sc.pp.subsample(adata, n_obs=40)
    assert adata.n_obs == 40
    sc.pp.subsample(adata, fraction=0.1)
    assert adata.n_obs == 4


def test_subsample_copy():
    adata = AnnData(np.ones((200, 10)))
    assert sc.pp.subsample(adata, n_obs=40, copy=True).shape == (40, 10)
    assert sc.pp.subsample(adata, fraction=0.1, copy=True).shape == (20, 10)


def test_subsample_copy_backed(tmp_path):
    A = np.random.rand(200, 10).astype(np.float32)
    adata_m = AnnData(A.copy())
    adata_d = AnnData(A.copy())
    filename = tmp_path / "test.h5ad"
    adata_d.filename = filename
    # This should not throw an error
    assert sc.pp.subsample(adata_d, n_obs=40, copy=True).shape == (40, 10)
    np.testing.assert_array_equal(
        sc.pp.subsample(adata_m, n_obs=40, copy=True).X,
        sc.pp.subsample(adata_d, n_obs=40, copy=True).X,
    )
    with pytest.raises(NotImplementedError):
        sc.pp.subsample(adata_d, n_obs=40, copy=False)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("zero_center", [True, False])
@pytest.mark.parametrize("max_value", [None, 1.0])
def test_scale_matrix_types(array_type, zero_center, max_value):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    sc.pp.scale(adata, zero_center=zero_center, max_value=max_value)
    sc.pp.scale(adata_casted, zero_center=zero_center, max_value=max_value)
    X = adata_casted.X
    if "dask" in array_type.__name__:
        X = X.compute()
    if issparse(X):
        X = X.todense()
    if issparse(adata.X):
        adata.X = adata.X.todense()
    assert_allclose(X, adata.X, rtol=1e-5, atol=1e-5)


ARRAY_TYPES_DASK_SPARSE = [
    a for a in ARRAY_TYPES if "sparse" in a.id and "dask" in a.id
]


@pytest.mark.parametrize("array_type", ARRAY_TYPES_DASK_SPARSE)
def test_scale_zero_center_warns_dask_sparse(array_type):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    with pytest.warns(UserWarning, match="zero-center being used with `DaskArray`*"):
        sc.pp.scale(adata_casted)
    sc.pp.scale(adata)
    assert_allclose(adata_casted.X, adata.X, rtol=1e-5, atol=1e-5)


def test_scale():
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    v = adata[:, 0 : adata.shape[1] // 2]
    # Should turn view to copy https://github.com/scverse/anndata/issues/171#issuecomment-508689965
    assert v.is_view
    with pytest.warns(Warning, match="view"):
        sc.pp.scale(v)
    assert not v.is_view
    assert_allclose(v.X.var(axis=0), np.ones(v.shape[1]), atol=0.01)
    assert_allclose(v.X.mean(axis=0), np.zeros(v.shape[1]), atol=0.00001)


@pytest.fixture(params=[True, False])
def zero_center(request):
    return request.param


def test_scale_rep(count_matrix_format, zero_center):
    """
    Test that it doesn't matter where the array being scaled is in the anndata object.
    """
    X = count_matrix_format(sp.random(100, 200, density=0.3).toarray())
    check_rep_mutation(sc.pp.scale, X, zero_center=zero_center)
    check_rep_results(sc.pp.scale, X, zero_center=zero_center)


def test_scale_array(count_matrix_format, zero_center):
    """
    Test that running sc.pp.scale on an anndata object and an array returns the same results.
    """
    X = count_matrix_format(sp.random(100, 200, density=0.3).toarray())
    adata = anndata_v0_8_constructor_compat(X=X.copy())

    sc.pp.scale(adata, zero_center=zero_center)
    scaled_X = sc.pp.scale(X, zero_center=zero_center, copy=True)
    np.testing.assert_equal(asarray(scaled_X), asarray(adata.X))


def test_recipe_plotting():
    sc.settings.autoshow = False
    adata = AnnData(np.random.randint(0, 1000, (1000, 1000)))
    # These shouldn't throw an error
    sc.pp.recipe_seurat(adata.copy(), plot=True)
    sc.pp.recipe_zheng17(adata.copy(), plot=True)


def test_regress_out_ordinal():
    from scipy.sparse import random

    adata = AnnData(random(1000, 100, density=0.6, format="csr"))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    # results using only one processor
    single = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )
    assert adata.X.shape == single.X.shape

    # results using 8 processors
    multi = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=8, copy=True
    )

    np.testing.assert_array_equal(single.X, multi.X)


def test_regress_out_layer():
    from scipy.sparse import random

    adata = AnnData(random(1000, 100, density=0.6, format="csr"))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.layers["counts"] = adata.X.copy()

    single = sc.pp.regress_out(
        adata, keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )
    assert adata.X.shape == single.X.shape

    layer = sc.pp.regress_out(
        adata, layer="counts", keys=["n_counts", "percent_mito"], n_jobs=1, copy=True
    )

    np.testing.assert_array_equal(single.X, layer.layers["counts"])


def test_regress_out_view():
    from scipy.sparse import random

    adata = AnnData(random(500, 1100, density=0.2, format="csr"))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    subset_adata = adata[:, :1050]
    subset_adata_copy = subset_adata.copy()

    sc.pp.regress_out(subset_adata, keys=["n_counts", "percent_mito"])
    sc.pp.regress_out(subset_adata_copy, keys=["n_counts", "percent_mito"])
    assert_equal(subset_adata, subset_adata_copy)
    assert not subset_adata.is_view


def test_regress_out_categorical():
    import pandas as pd
    from scipy.sparse import random

    adata = AnnData(random(1000, 100, density=0.6, format="csr"))
    # create a categorical column
    adata.obs["batch"] = pd.Categorical(np.random.randint(1, 4, size=adata.X.shape[0]))

    multi = sc.pp.regress_out(adata, keys="batch", n_jobs=8, copy=True)
    assert adata.X.shape == multi.X.shape


def test_regress_out_constants():
    adata = AnnData(np.hstack((np.full((10, 1), 0.0), np.full((10, 1), 1.0))))
    adata.obs["percent_mito"] = np.random.rand(adata.X.shape[0])
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata_copy = adata.copy()

    sc.pp.regress_out(adata, keys=["n_counts", "percent_mito"])
    assert_equal(adata, adata_copy)


def test_regress_out_constants_equivalent():
    # Tests that constant values don't change results
    # (since support for constant values is implemented by us)
    from sklearn.datasets import make_blobs

    X, cat = make_blobs(100, 20)
    a = sc.AnnData(np.hstack([X, np.zeros((100, 5))]), obs={"cat": pd.Categorical(cat)})
    b = sc.AnnData(X, obs={"cat": pd.Categorical(cat)})

    sc.pp.regress_out(a, "cat")
    sc.pp.regress_out(b, "cat")

    np.testing.assert_equal(a[:, b.var_names].X, b.X)


@pytest.fixture(params=[lambda x: x.copy(), sp.csr_matrix, sp.csc_matrix])
def count_matrix_format(request):
    return request.param


@pytest.fixture(params=[True, False])
def replace(request):
    return request.param


@pytest.fixture(params=[np.int64, np.float32, np.float64])
def dtype(request):
    return request.param


def test_downsample_counts_per_cell(count_matrix_format, replace, dtype):
    TARGET = 1000
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata = anndata_v0_8_constructor_compat(X=count_matrix_format(X).astype(dtype))
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(
            adata, counts_per_cell=TARGET, total_counts=TARGET, replace=replace
        )
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(adata, replace=replace)
    initial_totals = np.ravel(adata.X.sum(axis=1))
    adata = sc.pp.downsample_counts(
        adata, counts_per_cell=TARGET, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGET)
    assert all(initial_totals >= new_totals)
    assert all(
        initial_totals[initial_totals <= TARGET] == new_totals[initial_totals <= TARGET]
    )
    if not replace:
        assert np.all(X >= adata.X)
    assert X.dtype == adata.X.dtype


def test_downsample_counts_per_cell_multiple_targets(
    count_matrix_format, replace, dtype
):
    TARGETS = np.random.randint(500, 1500, 1000)
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata = anndata_v0_8_constructor_compat(X=count_matrix_format(X).astype(dtype))
    initial_totals = np.ravel(adata.X.sum(axis=1))
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(adata, counts_per_cell=[40, 10], replace=replace)
    adata = sc.pp.downsample_counts(
        adata, counts_per_cell=TARGETS, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGETS)
    assert all(initial_totals >= new_totals)
    assert all(
        initial_totals[initial_totals <= TARGETS]
        == new_totals[initial_totals <= TARGETS]
    )
    if not replace:
        assert np.all(X >= adata.X)
    assert X.dtype == adata.X.dtype


def test_downsample_total_counts(count_matrix_format, replace, dtype):
    X = np.random.randint(0, 100, (1000, 100)) * np.random.binomial(1, 0.3, (1000, 100))
    X = X.astype(dtype)
    adata_orig = anndata_v0_8_constructor_compat(X=count_matrix_format(X))
    total = X.sum()
    target = np.floor_divide(total, 10)
    initial_totals = np.ravel(adata_orig.X.sum(axis=1))
    adata = sc.pp.downsample_counts(
        adata_orig, total_counts=target, replace=replace, copy=True
    )
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert adata.X.sum() == target
    assert all(initial_totals >= new_totals)
    if not replace:
        assert np.all(X >= adata.X)
        adata = sc.pp.downsample_counts(
            adata_orig, total_counts=total + 10, replace=False, copy=True
        )
        assert (adata.X == X).all()
    assert X.dtype == adata.X.dtype


def test_recipe_weinreb():
    # Just tests for failure for now
    adata = pbmc68k_reduced().raw.to_adata()
    adata.X = adata.X.toarray()

    orig = adata.copy()
    sc.pp.recipe_weinreb17(adata, log=False, copy=True)
    assert_equal(orig, adata)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize(
    "max_cells,max_counts,min_cells,min_counts",
    [
        [100, None, None, None],
        [None, 100, None, None],
        [None, None, 20, None],
        [None, None, None, 20],
    ],
)
def test_filter_genes(array_type, max_cells, max_counts, min_cells, min_counts):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    sc.pp.filter_genes(
        adata,
        max_cells=max_cells,
        max_counts=max_counts,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    sc.pp.filter_genes(
        adata_casted,
        max_cells=max_cells,
        max_counts=max_counts,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    X = adata_casted.X
    if "dask" in array_type.__name__:
        X = X.compute()
    if issparse(X):
        X = X.todense()
    if issparse(adata.X):
        adata.X = adata.X.todense()
    assert_allclose(X, adata.X, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize(
    "max_genes,max_counts,min_genes,min_counts",
    [
        [100, None, None, None],
        [None, 100, None, None],
        [None, None, 20, None],
        [None, None, None, 20],
    ],
)
def test_filter_cells(array_type, max_genes, max_counts, min_genes, min_counts):
    adata = pbmc68k_reduced()
    adata.X = adata.raw.X
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X)
    sc.pp.filter_cells(
        adata,
        max_genes=max_genes,
        max_counts=max_counts,
        min_genes=min_genes,
        min_counts=min_counts,
    )
    sc.pp.filter_cells(
        adata_casted,
        max_genes=max_genes,
        max_counts=max_counts,
        min_genes=min_genes,
        min_counts=min_counts,
    )
    X = adata_casted.X
    if "dask" in array_type.__name__:
        X = X.compute()
    if issparse(X):
        X = X.todense()
    if issparse(adata.X):
        adata.X = adata.X.todense()
    assert_allclose(X, adata.X, rtol=1e-5, atol=1e-5)
