import pytest
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix

import scanpy as sc

# test "data" for 3 cells * 4 genes
X = [
    [-1, 2, 0, 0],
    [1, 2, 4, 0],
    [0, 2, 2, 0],
]  # with gene std 1,0,2,0 and center 0,2,2,0
X_scaled = [
    [-1, 2, 0, 0],
    [1, 2, 2, 0],
    [0, 2, 1, 0],
]  # with gene std 1,0,1,0 and center 0,2,1,0
X_centered = [
    [-1, 0, -1, 0],
    [1, 0, 1, 0],
    [0, 0, 0, 0],
]  # with gene std 1,0,1,0 and center 0,0,0,0


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_scale(typ, dtype):
    ## test AnnData arguments
    # test scaling with default zero_center == True
    adata0 = AnnData(typ(X), dtype=dtype)
    sc.pp.scale(adata0)
    assert np.allclose(csr_matrix(adata0.X).toarray(), X_centered)
    # test scaling with explicit zero_center == True
    adata1 = AnnData(typ(X), dtype=dtype)
    sc.pp.scale(adata1, zero_center=True)
    assert np.allclose(csr_matrix(adata1.X).toarray(), X_centered)
    # test scaling with explicit zero_center == False
    adata2 = AnnData(typ(X), dtype=dtype)
    sc.pp.scale(adata2, zero_center=False)
    assert np.allclose(csr_matrix(adata2.X).toarray(), X_scaled)
    ## test bare count arguments, for simplicity only with explicit copy=True
    # test scaling with default zero_center == True
    data0 = typ(X, dtype=dtype)
    cdata0 = sc.pp.scale(data0, copy=True)
    assert np.allclose(csr_matrix(cdata0).toarray(), X_centered)
    # test scaling with explicit zero_center == True
    data1 = typ(X, dtype=dtype)
    cdata1 = sc.pp.scale(data1, zero_center=True, copy=True)
    assert np.allclose(csr_matrix(cdata1).toarray(), X_centered)
    # test scaling with explicit zero_center == False
    data2 = typ(X, dtype=dtype)
    cdata2 = sc.pp.scale(data2, zero_center=False, copy=True)
    assert np.allclose(csr_matrix(cdata2).toarray(), X_scaled)
