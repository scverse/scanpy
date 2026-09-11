from __future__ import annotations

import warnings
from contextlib import nullcontext

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from testing.scanpy._pytest.marks import needs

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
    ("mask", "x", "x_centered", "x_scaled"),
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
def test_scale(*, typ, container, zero_center, dtype, mask, x, x_centered, x_scaled):
    x = AnnData(typ(x, dtype=dtype)) if container == "anndata" else typ(x, dtype=dtype)
    with warnings.catch_warnings():
        # TODO: fix setting slices of sparse matrices in scale()
        warnings.filterwarnings("always", category=sparse.SparseEfficiencyWarning)

        with (
            pytest.warns(UserWarning, match=r"zero-center.*densifies")
            if zero_center and any(f in typ.__name__ for f in ("csr", "csc"))
            else nullcontext()
        ):
            scaled = sc.pp.scale(
                x, zero_center=zero_center, copy=container == "array", mask=mask
            )
    received = sparse.csr_matrix(  # noqa: TID251
        x.X if scaled is None else scaled
    ).toarray()
    expected = x_centered if zero_center else x_scaled
    assert np.allclose(received, expected)


@pytest.fixture
def adata_masked() -> AnnData:
    adata = AnnData(np.array(X_for_mask, dtype="float32"))
    adata.obs["some cells"] = np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool)
    return adata


@needs.anndata_acc
@pytest.mark.parametrize("as_ref", [True, False], ids=["ref", "str"])
def test_mask_ref(adata_masked: AnnData, *, as_ref: bool) -> None:
    from anndata.acc import A

    mask = A.obs["some cells"] if as_ref else "obs.some cells"
    with pytest.raises(ValueError, match=r"Cannot.*refer.*mask.*without.*anndata"):
        sc.pp.scale(np.array(X_original), mask=mask)
    sc.pp.scale(adata_masked, mask=mask)
    assert np.array_equal(adata_masked.X, X_centered_for_mask)
    assert "mean of A.obs['some cells']" in adata_masked.var.columns


@needs.anndata_acc
def test_mask_wrong_dim(adata_masked: AnnData) -> None:
    with pytest.raises(ValueError, match=r"Dimension of .* \(var\) does not match"):
        sc.pp.scale(adata_masked, mask="var.some genes")


@pytest.mark.parametrize(
    ("preset", "col"),
    [
        pytest.param(sc.Preset.ScanpyV1, "mean of some cells", id="v1"),
        pytest.param(
            sc.Preset.ScanpyV2Preview,
            "mean of A.obs['some cells']",
            id="v2",
            marks=needs.scanpy2,
        ),
    ],
)
def test_mask_obs_deprecated(
    adata_masked: AnnData, *, preset: sc.Preset, col: str
) -> None:
    r"""The deprecated `mask_obs` interprets :class:`str`\ s as `.obs` column names."""
    with (
        sc.settings.override(preset=preset),
        pytest.warns(FutureWarning, match=r"argument mask_obs is deprecated"),
    ):
        sc.pp.scale(adata_masked, zero_center=True, mask_obs="some cells")
    assert np.array_equal(adata_masked.X, X_centered_for_mask)
    assert col in adata_masked.var.columns


def test_mask_obs_deprecated_fallback() -> None:
    # extra test for the singledispatch’s fallback branch calling `scale_array`
    # (registered types like `np.ndarray` never run it)
    mask = np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool)
    scale_fallback = sc.pp.scale.dispatch(object)
    with pytest.warns(FutureWarning, match=r"argument mask_obs is deprecated"):
        scaled = scale_fallback(
            np.array(X_for_mask, dtype="float32"), copy=True, mask_obs=mask
        )
    assert np.array_equal(scaled, X_centered_for_mask)


def test_mask_both() -> None:
    mask = np.array((0, 0, 1, 1, 1, 0, 0), dtype=bool)
    with (
        pytest.warns(FutureWarning, match=r"argument mask_obs is deprecated"),
        pytest.raises(TypeError, match=r"Pass either `mask` or `mask_obs`, not both"),
    ):
        sc.pp.scale(np.array(X_for_mask, dtype="float32"), mask=mask, mask_obs=mask)


@pytest.mark.parametrize("zero_center", [True, False], ids=["center", "no_center"])
def test_clip(*, zero_center: bool) -> None:
    adata = sc.datasets.pbmc3k()
    with (
        (pytest.warns(UserWarning, match=r"zero-center.*densifies"))
        if zero_center
        else nullcontext()
    ):
        sc.pp.scale(adata, max_value=1, zero_center=zero_center)
    if zero_center:
        assert adata.X.min() >= -1
    assert adata.X.max() <= 1


@pytest.mark.parametrize(
    ("mask", "x", "x_scaled", "x_clipped"),
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
def test_scale_sparse(*, mask, x, x_scaled, x_clipped, clip):
    max_value, expected = (1, x_clipped) if clip else (None, x_scaled)
    adata = AnnData(sparse.csr_matrix(x).astype(np.float32))  # noqa: TID251
    sc.pp.scale(adata, mask=mask, zero_center=False, max_value=max_value)
    assert np.allclose(sparse.csr_matrix(adata.X).toarray(), expected)  # noqa: TID251
