from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

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
    "typ",
    [np.array, sparse.csr_matrix, sparse.csc_matrix],  # noqa: TID251
    ids=lambda x: x.__name__,
)
@pytest.mark.parametrize("container", ["anndata", "array"])
@pytest.mark.parametrize("dtype", [np.float32, np.int64])
@pytest.mark.parametrize("zero_center", [True, False], ids=["center", "no_center"])
@pytest.mark.parametrize(
    ("mask_obs", "x", "x_centered", "x_scaled"),
    [
        pytest.param(
            None, X_original, X_centered_original, X_scaled_original, id="no_mask"
        ),
        pytest.param(
            np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool),
            X_for_mask,
            X_centered_for_mask,
            X_scaled_for_mask,
            id="mask",
        ),
    ],
)
def test_scale(
    *, typ, container, zero_center, dtype, mask_obs, x, x_centered, x_scaled
):
    x = AnnData(typ(x, dtype=dtype)) if container == "anndata" else typ(x, dtype=dtype)
    scaled = sc.pp.scale(
        x, zero_center=zero_center, copy=container == "array", mask_obs=mask_obs
    )
    received = sparse.csr_matrix(  # noqa: TID251
        x.X if scaled is None else scaled
    ).toarray()
    expected = x_centered if zero_center else x_scaled
    assert np.allclose(received, expected)


def test_mask_string():
    with pytest.raises(ValueError, match=r"Cannot refer to mask.* without.*anndata"):
        sc.pp.scale(np.array(X_original), mask_obs="mask")
    adata = AnnData(np.array(X_for_mask, dtype="float32"))
    adata.obs["some cells"] = np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool)
    sc.pp.scale(adata, mask_obs="some cells")
    assert np.array_equal(adata.X, X_centered_for_mask)
    assert "mean of some cells" in adata.var.columns


@pytest.mark.parametrize("zero_center", [True, False])
def test_clip(zero_center):
    adata = sc.datasets.pbmc3k()
    sc.pp.scale(adata, max_value=1, zero_center=zero_center)
    if zero_center:
        assert adata.X.min() >= -1
    assert adata.X.max() <= 1


@pytest.mark.parametrize(
    ("mask_obs", "x", "x_scaled", "x_clipped"),
    [
        pytest.param(
            None, X_original, X_scaled_original, X_scaled_original_clipped, id="no_mask"
        ),
        pytest.param(
            np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool),
            X_for_mask,
            X_scaled_for_mask,
            X_scaled_for_mask_clipped,
            id="mask",
        ),
    ],
)
@pytest.mark.parametrize("clip", [False, True], ids=["no_clip", "clip"])
def test_scale_sparse(*, mask_obs, x, x_scaled, x_clipped, clip):
    max_value, expected = (1, x_clipped) if clip else (None, x_scaled)
    adata = AnnData(sparse.csr_matrix(x).astype(np.float32))  # noqa: TID251
    sc.pp.scale(adata, mask_obs=mask_obs, zero_center=False, max_value=max_value)
    assert np.allclose(sparse.csr_matrix(adata.X).toarray(), expected)  # noqa: TID251
