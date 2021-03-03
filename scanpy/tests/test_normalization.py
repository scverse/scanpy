import pytest
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy import sparse

import scanpy as sc
from scanpy.tests.helpers import check_rep_mutation, check_rep_results
from anndata.tests.helpers import assert_equal, asarray

X_total = [[1, 0], [3, 0], [5, 6]]
X_frac = [[1, 0, 1], [3, 0, 1], [5, 6, 1]]


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    sc.pp.normalize_total(adata, key_added='n_counts')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [3.0, 3.0, 3.0])
    sc.pp.normalize_total(adata, target_sum=1, key_added='n_counts2')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [1.0, 1.0, 1.0])

    adata = AnnData(typ(X_frac, dtype=dtype))
    sc.pp.normalize_total(adata, exclude_highly_expressed=True, max_fraction=0.7)
    assert np.allclose(np.ravel(adata.X[:, 1:3].sum(axis=1)), [1.0, 1.0, 1.0])


@pytest.mark.parametrize('typ', [asarray, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_rep(typ, dtype):
    # Test that layer kwarg works
    X = typ(sparse.random(100, 50, format="csr", density=0.2, dtype=dtype))
    check_rep_mutation(sc.pp.normalize_total, X, fields=["layer"])
    check_rep_results(sc.pp.normalize_total, X, fields=["layer"])


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_layers(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    adata.layers["layer"] = adata.X.copy()
    with pytest.warns(FutureWarning, match=r".*layers.*deprecated"):
        sc.pp.normalize_total(adata, layers=["layer"])
    assert np.allclose(adata.layers["layer"].sum(axis=1), [3.0, 3.0, 3.0])


@pytest.mark.parametrize('typ', [np.array, csr_matrix], ids=lambda x: x.__name__)
@pytest.mark.parametrize('dtype', ['float32', 'int64'])
def test_normalize_total_view(typ, dtype):
    adata = AnnData(typ(X_total), dtype=dtype)
    v = adata[:, :]

    sc.pp.normalize_total(v)
    sc.pp.normalize_total(adata)

    assert not v.is_view
    assert_equal(adata, v)
