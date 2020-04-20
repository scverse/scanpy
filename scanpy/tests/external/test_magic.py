import pytest
import numpy as np
from anndata import AnnData

import scanpy as sc

pytest.importorskip("magic", minversion=sc.external.pp._magic.MIN_VERSION)

A_list = [
    [0, 0, 7, 0, 0],
    [8, 5, 0, 2, 0],
    [6, 0, 0, 2, 5],
    [0, 0, 0, 1, 0],
    [8, 8, 2, 1, 0],
    [0, 0, 0, 4, 5],
]


def test_magic_default():
    A = np.array(A_list, dtype='float32')
    adata = AnnData(A)
    sc.external.pp.magic(adata, knn=1)
    # check raw unchanged
    np.testing.assert_array_equal(adata.raw.X, A)
    # check .X changed
    assert not np.all(adata.X == A)
    # check .X shape unchanged
    assert adata.X.shape == A.shape


def test_magic_pca_only():
    A = np.array(A_list, dtype='float32')
    # pca only
    adata = AnnData(A)
    n_pca = 3
    sc.external.pp.magic(adata, knn=1, name_list='pca_only', n_pca=n_pca)
    # check raw unchanged
    np.testing.assert_array_equal(adata.X, A)
    # check .X shape consistent with n_pca
    assert adata.obsm['X_magic'].shape == (A.shape[0], n_pca)


def test_magic_copy():
    A = np.array(A_list, dtype='float32')
    adata = AnnData(A)
    adata_copy = sc.external.pp.magic(adata, knn=1, copy=True)
    # check adata unchanged
    np.testing.assert_array_equal(adata.X, A)
    # check copy raw unchanged
    np.testing.assert_array_equal(adata_copy.raw.X, A)
    # check .X changed
    assert not np.all(adata_copy.X == A)
    # check .X shape unchanged
    assert adata_copy.X.shape == A.shape
