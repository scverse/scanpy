from itertools import product
import numpy as np
from scipy import sparse as sp
import scanpy as sc
import pytest
from anndata import AnnData


def test_log1p_chunked():
    A = np.random.rand(200, 10)
    ad = AnnData(A)
    ad2 = AnnData(A)
    ad3 = AnnData(A)
    ad3.filename = 'test.h5ad'
    sc.pp.log1p(ad)
    sc.pp.log1p(ad2, chunked=True)
    assert np.allclose(ad2.X, ad.X)
    sc.pp.log1p(ad3, chunked=True)
    assert np.allclose(ad3.X, ad.X)


def test_normalize_per_cell():
    adata = AnnData(
        np.array([[1, 0], [3, 0], [5, 6]]))
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1,
                             key_n_counts='n_counts2')
    assert adata.X.sum(axis=1).tolist() == [1., 1., 1.]
    # now with copy option
    adata = AnnData(
        np.array([[1, 0], [3, 0], [5, 6]]))
    # note that sc.pp.normalize_per_cell is also used in
    # pl.highest_expr_genes with parameter counts_per_cell_after=100
    adata_copy = sc.pp.normalize_per_cell(
        adata, counts_per_cell_after=1, copy=True)
    assert adata_copy.X.sum(axis=1).tolist() == [1., 1., 1.]
    # now sparse
    adata = AnnData(
        np.array([[1, 0], [3, 0], [5, 6]]))
    adata_sparse = AnnData(
        sp.csr_matrix([[1, 0], [3, 0], [5, 6]]))
    sc.pp.normalize_per_cell(adata)
    sc.pp.normalize_per_cell(adata_sparse)
    assert adata.X.sum(axis=1).tolist() == adata_sparse.X.sum(
        axis=1).A1.tolist()


def test_subsample():
    adata = AnnData(np.ones((200, 10)))
    sc.pp.subsample(adata, n_obs=40)
    assert adata.n_obs == 40
    sc.pp.subsample(adata, fraction=0.1)
    assert adata.n_obs == 4


def test_recipe_plotting():
    sc.settings.autoshow = False
    adata = AnnData(np.random.randint(0, 1000, (1000, 1000)))
    # These shouldn't throw an error
    sc.pp.recipe_seurat(adata.copy(), plot=True)
    sc.pp.recipe_zheng17(adata.copy(), plot=True)


def test_regress_out_ordinal():
    from scipy.sparse import random
    adata = AnnData(random(1000, 100, density=0.6, format='csr'))
    adata.obs['percent_mito'] = np.random.rand(adata.X.shape[0])
    adata.obs['n_counts'] = adata.X.sum(axis=1)

    # results using only one processor
    single = sc.pp.regress_out(
        adata, keys=['n_counts', 'percent_mito'], n_jobs=1, copy=True)
    assert adata.X.shape == single.X.shape

    # results using 8 processors
    multi = sc.pp.regress_out(
        adata, keys=['n_counts', 'percent_mito'], n_jobs=8, copy=True)

    np.testing.assert_array_equal(single.X, multi.X)


def test_regress_out_categorical():
    from scipy.sparse import random
    import pandas as pd
    adata = AnnData(random(1000, 100, density=0.6, format='csr'))
    # create a categorical column
    adata.obs['batch'] = pd.Categorical(
        np.random.randint(1, 4, size=adata.X.shape[0]))

    multi = sc.pp.regress_out(adata, keys='batch', n_jobs=8, copy=True)
    assert adata.X.shape == multi.X.shape


@pytest.fixture(params=[lambda x: x.copy(), sp.csr_matrix, sp.csc_matrix])
def count_matrix_format(request):
    return request.param


@pytest.fixture(params=[True, False])
def replace(request):
    return request.param


def test_downsample_counts_per_cell(count_matrix_format, replace):
    TARGET = 1000
    X = np.random.randint(0, 100, (1000, 100)) * \
        np.random.binomial(1, .3, (1000, 100))
    adata = AnnData(X=count_matrix_format(X))
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(adata, counts_per_cell=TARGET, total_counts=TARGET, replace=replace)
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(adata, replace=replace)
    initial_totals = np.ravel(adata.X.sum(axis=1))
    adata = sc.pp.downsample_counts(adata, counts_per_cell=TARGET,
                                    replace=replace, copy=True)
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGET)
    assert all(initial_totals >= new_totals)
    assert all(initial_totals[initial_totals <= TARGET]
                == new_totals[initial_totals <= TARGET])
    if not replace:
        assert np.all(X >= adata.X)


def test_downsample_counts_per_cell_multiple_targets(count_matrix_format, replace):
    TARGETS = np.random.randint(500, 1500, 1000)
    X = np.random.randint(0, 100, (1000, 100)) * \
        np.random.binomial(1, .3, (1000, 100))
    adata = AnnData(X=count_matrix_format(X))
    initial_totals = np.ravel(adata.X.sum(axis=1))
    with pytest.raises(ValueError):
        sc.pp.downsample_counts(adata, counts_per_cell=[40, 10], replace=replace)
    adata = sc.pp.downsample_counts(adata, counts_per_cell=TARGETS, replace=replace, copy=True)
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert all(new_totals <= TARGETS)
    assert all(initial_totals >= new_totals)
    assert all(initial_totals[initial_totals <= TARGETS]
                == new_totals[initial_totals <= TARGETS])
    if not replace:
        assert np.all(X >= adata.X)


def test_downsample_total_counts(count_matrix_format, replace):
    X = np.random.randint(0, 100, (1000, 100)) * \
        np.random.binomial(1, .3, (1000, 100))
    adata_orig = AnnData(X=count_matrix_format(X))
    total = X.sum()
    target = np.floor_divide(total, 10)
    initial_totals = np.ravel(adata_orig.X.sum(axis=1))
    adata = sc.pp.downsample_counts(adata_orig, total_counts=target,
                                    replace=replace, copy=True)
    new_totals = np.ravel(adata.X.sum(axis=1))
    if sp.issparse(adata.X):
        assert all(adata.X.toarray()[X == 0] == 0)
    else:
        assert all(adata.X[X == 0] == 0)
    assert adata.X.sum() == target
    assert all(initial_totals >= new_totals)
    if not replace:
        assert np.all(X >= adata.X)
        adata = sc.pp.downsample_counts(adata_orig, total_counts=total + 10,
                                        replace=False, copy=True)
        assert (adata.X == X).all()
