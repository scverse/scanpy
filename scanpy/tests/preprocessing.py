import numpy as np
from numpy import ma
from scipy import sparse as sp
# we don’t need this in requirements.txt, as it’s only needed for testing
from pytest import mark
import scanpy.api as sc

from anndata import AnnData


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


def test_recipe_plotting():
    sc.settings.autoshow = False
    adata = AnnData(np.random.randint(0, 1000, (1000, 1000)))
    # These shouldn't throw an error
    sc.pp.recipe_seurat(adata.copy(), plot=True)
    sc.pp.recipe_zheng17(adata.copy(), plot=True)


def test_regress_out_ordinal():
    from scipy.sparse import random
    adata = AnnData(random(5000, 2000, density=0.6, format='csr'))
    adata.obs['percent_mito'] = np.random.rand(adata.X.shape[0])
    adata.obs['n_counts'] = adata.X.sum(axis=1)

    # results using only one processor
    single = sc.pp.regress_out(adata, keys=['n_counts', 'percent_mito'], n_jobs=1, copy=True)
    assert adata.X.shape == single.X.shape

    # results using 8 processors
    multi = sc.pp.regress_out(adata, keys=['n_counts', 'percent_mito'], n_jobs=8, copy=True)

    np.testing.assert_array_equal(single.X, multi.X)


def test_regress_out_categorical():
    from scipy.sparse import random
    import pandas as pd
    adata = AnnData(random(2500, 1000, density=0.6, format='csr'))
    # create a categorical column
    adata.obs['batch'] = pd.Categorical(np.random.randint(1, 10, size=adata.X.shape[0]))

    multi = sc.pp.regress_out(adata, keys='batch', n_jobs=8, copy=True)
    assert adata.X.shape == multi.X.shape
