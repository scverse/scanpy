import pytest
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix

import scanpy as sc

X_total = [[1, 0], [3, 0], [5, 6]]
X_frac = [[1, 0, 1], [3, 0, 1], [5, 6, 1]]


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_normalize_total(typ):
    adata = AnnData(typ(X_total, dtype='float32'))
    sc.pp.normalize_total(adata, key_added='n_counts')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [3.0, 3.0, 3.0])
    sc.pp.normalize_total(adata, target_sum=1, key_added='n_counts2')
    assert np.allclose(np.ravel(adata.X.sum(axis=1)), [1.0, 1.0, 1.0])

    adata = AnnData(typ(X_frac, dtype='float32'))
    sc.pp.normalize_total(
        adata, exclude_highly_expressed=True, max_fraction=0.7
    )
    assert np.allclose(np.ravel(adata.X[:, 1:3].sum(axis=1)), [1.0, 1.0, 1.0])
