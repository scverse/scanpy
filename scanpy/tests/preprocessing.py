import numpy as np
from scipy import sparse as sp
import scanpy.api as sc
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
