from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy.sparse import csc_matrix, csr_matrix

import scanpy as sc

# test "data" for 3 cells * 4 genes
X_original = [
    [-1, 2, 0, 0],
    [1, 2, 4, 0],
    [0, 2, 2, 0],
]  # with gene std 1,0,2,0 and center 0,2,2,0
X_scaled_original = [
    [-1, 2, 0, 0],
    [1, 2, 2, 0],
    [0, 2, 1, 0],
]  # with gene std 1,0,1,0 and center 0,2,1,0
X_centered_original = [
    [-1, 0, -1, 0],
    [1, 0, 1, 0],
    [0, 0, 0, 0],
]  # with gene std 1,0,1,0 and center 0,0,0,0
X_scaled_original_clipped = [
    [-1, 1, 0, 0],
    [1, 1, 1, 0],
    [0, 1, 1, 0],
]  # with gene std 1,0,1,0 and center 0,2,1,0


X_for_mask = [
    [27, 27, 27, 27],
    [27, 27, 27, 27],
    [-1, 2, 0, 0],
    [1, 2, 4, 0],
    [0, 2, 2, 0],
    [27, 27, 27, 27],
    [27, 27, 27, 27],
]
X_scaled_for_mask = [
    [27, 27, 27, 27],
    [27, 27, 27, 27],
    [-1, 2, 0, 0],
    [1, 2, 2, 0],
    [0, 2, 1, 0],
    [27, 27, 27, 27],
    [27, 27, 27, 27],
]
X_centered_for_mask = [
    [27, 27, 27, 27],
    [27, 27, 27, 27],
    [-1, 0, -1, 0],
    [1, 0, 1, 0],
    [0, 0, 0, 0],
    [27, 27, 27, 27],
    [27, 27, 27, 27],
]
X_scaled_for_mask_clipped = [
    [27, 27, 27, 27],
    [27, 27, 27, 27],
    [-1, 1, 0, 0],
    [1, 1, 1, 0],
    [0, 1, 1, 0],
    [27, 27, 27, 27],
    [27, 27, 27, 27],
]


@pytest.mark.parametrize(
    "typ", [np.array, csr_matrix, csc_matrix], ids=lambda x: x.__name__
)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
@pytest.mark.parametrize(
    ("mask_obs", "X", "X_centered", "X_scaled"),
    [
        (None, X_original, X_centered_original, X_scaled_original),
        (
            np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool),
            X_for_mask,
            X_centered_for_mask,
            X_scaled_for_mask,
        ),
    ],
)
def test_scale(*, typ, dtype, mask_obs, X, X_centered, X_scaled):
    # test AnnData arguments
    # test scaling with default zero_center == True
    adata0 = AnnData(typ(X).astype(dtype))
    sc.pp.scale(adata0, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(adata0.X).toarray(), X_centered)
    # test scaling with explicit zero_center == True
    adata1 = AnnData(typ(X).astype(dtype))
    sc.pp.scale(adata1, zero_center=True, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(adata1.X).toarray(), X_centered)
    # test scaling with explicit zero_center == False
    adata2 = AnnData(typ(X).astype(dtype))
    sc.pp.scale(adata2, zero_center=False, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(adata2.X).toarray(), X_scaled)
    # test bare count arguments, for simplicity only with explicit copy=True
    # test scaling with default zero_center == True
    data0 = typ(X, dtype=dtype)
    cdata0 = sc.pp.scale(data0, copy=True, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(cdata0).toarray(), X_centered)
    # test scaling with explicit zero_center == True
    data1 = typ(X, dtype=dtype)
    cdata1 = sc.pp.scale(data1, zero_center=True, copy=True, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(cdata1).toarray(), X_centered)
    # test scaling with explicit zero_center == False
    data2 = typ(X, dtype=dtype)
    cdata2 = sc.pp.scale(data2, zero_center=False, copy=True, mask_obs=mask_obs)
    assert np.allclose(csr_matrix(cdata2).toarray(), X_scaled)


def test_mask_string():
    with pytest.raises(ValueError):
        sc.pp.scale(np.array(X_original), mask_obs="mask")
    adata = AnnData(np.array(X_for_mask, dtype="float32"))
    adata.obs["some cells"] = np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool)
    sc.pp.scale(adata, mask_obs="some cells")
    assert np.array_equal(adata.X, X_centered_for_mask)
    assert "mean of some cells" in adata.var.keys()


@pytest.mark.parametrize("zero_center", [True, False])
def test_clip(zero_center):
    adata = sc.datasets.pbmc3k()
    sc.pp.scale(adata, max_value=1, zero_center=zero_center)
    if zero_center:
        assert adata.X.min() >= -1
    assert adata.X.max() <= 1


@pytest.mark.parametrize(
    ("mask_obs", "X", "X_scaled", "X_clipped"),
    [
        (None, X_original, X_scaled_original, X_scaled_original_clipped),
        (
            np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool),
            X_for_mask,
            X_scaled_for_mask,
            X_scaled_for_mask_clipped,
        ),
    ],
)
def test_scale_sparse(*, mask_obs, X, X_scaled, X_clipped):
    adata0 = AnnData(csr_matrix(X).astype(np.float32))
    sc.pp.scale(adata0, mask_obs=mask_obs, zero_center=False)
    assert np.allclose(csr_matrix(adata0.X).toarray(), X_scaled)
    # test scaling with explicit zero_center == True
    adata1 = AnnData(csr_matrix(X).astype(np.float32))
    sc.pp.scale(adata1, zero_center=False, mask_obs=mask_obs, max_value=1)
    assert np.allclose(csr_matrix(adata1.X).toarray(), X_clipped)
