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
