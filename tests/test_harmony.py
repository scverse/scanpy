from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pooch
import pytest
from anndata import AnnData
from scipy.stats import pearsonr

from scanpy.preprocessing import harmony_integrate
from scanpy.preprocessing._harmony.core import _SUPPRESS_PENALTY, _compute_lambda_kb
from testing.scanpy._helpers.data import pbmc68k_reduced

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike


_HARMONY_DATA_BASE = "https://exampledata.scverse.org/rapids-singlecell/harmony_data"

DATA = dict(
    pca=("pbmc_3500_pcs.tsv.gz", "md5:27e319b3ddcc0c00d98e70aa8e677b10"),
    pca_harmonized=(
        "pbmc_3500_pcs_harmonized.tsv.gz",
        "md5:a7c4ce4b98c390997c66d63d48e09221",
    ),
    meta=("pbmc_3500_meta.tsv.gz", "md5:8c7ca20e926513da7cf0def1211baecb"),
)


def _get_measure(
    x: np.ndarray, base: np.ndarray, norm: Literal["r", "L2"]
) -> np.ndarray:
    """Compute correlation or L2 distance between arrays."""
    if norm == "r":  # Compute per-column correlation
        if x.ndim == 1:
            corr, _ = pearsonr(x, base)
            return corr
        return np.array([pearsonr(x[:, i], base[:, i])[0] for i in range(x.shape[1])])
    if norm == "L2":
        # L2 distance normalized by base norm
        if x.ndim == 1:
            return np.linalg.norm(x - base) / np.linalg.norm(base)
        return np.array([
            np.linalg.norm(x[:, i] - base[:, i]) / np.linalg.norm(base[:, i])
            for i in range(x.shape[1])
        ])
    pytest.fail(f"Unknown {norm=!r}")


@pytest.fixture(scope="module")
def _adata_reference() -> AnnData:
    """Load reference data once per module (avoids re-reading CSV)."""
    paths = {
        f: pooch.retrieve(
            f"{_HARMONY_DATA_BASE}/{name}",
            known_hash=hash_,
        )
        for f, (name, hash_) in DATA.items()
    }
    dfs = {f: pd.read_csv(path, delimiter="\t") for f, path in paths.items()}
    # Create unique index using row number + cell name
    dfs["meta"].index = [f"{i}_{cell}" for i, cell in enumerate(dfs["meta"]["cell"])]
    return AnnData(
        X=None,
        obs=dfs["meta"],
        obsm={"X_pca": dfs["pca"].values, "harmony_org": dfs["pca_harmonized"].values},
    )


@pytest.fixture
def adata_reference(_adata_reference: AnnData) -> AnnData:
    """Return a fresh copy per test so tests don't mutate shared state."""
    return _adata_reference.copy()


@pytest.mark.parametrize("correction_method", ["fast", "original"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("flavor", ["harmony1", "harmony2"])
def test_harmony_integrate(
    correction_method: Literal["fast", "original"],
    dtype: DTypeLike,
    flavor: Literal["harmony1", "harmony2"],
) -> None:
    """Test that Harmony integrate works."""
    adata = pbmc68k_reduced()
    harmony_integrate(
        adata,
        "bulk_labels",
        correction_method=correction_method,
        dtype=dtype,
        flavor=flavor,
    )
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate_algos(subtests: pytest.Subtests, dtype: DTypeLike) -> None:
    """Test that both correction methods produce similar results."""
    adata = pbmc68k_reduced()

    harmony_integrate(
        adata, "bulk_labels", correction_method="fast", dtype=dtype, flavor="harmony1"
    )
    fast = adata.obsm["X_pca_harmony"].copy()
    harmony_integrate(
        adata,
        "bulk_labels",
        correction_method="original",
        dtype=dtype,
        flavor="harmony1",
    )
    slow = adata.obsm["X_pca_harmony"].copy()

    with subtests.test("r"):
        assert _get_measure(fast, slow, "r").min() > 0.99
    with subtests.test("L2"):
        assert _get_measure(fast, slow, "L2").max() < 0.1


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("correction_method", ["fast", "original"])
def test_harmony_integrate_reference(
    *,
    subtests: pytest.Subtests,
    adata_reference: AnnData,
    dtype: DTypeLike,
    correction_method: Literal["fast", "original"],
) -> None:
    """Test that Harmony1 produces results similar to the reference implementation."""
    harmony_integrate(
        adata_reference,
        "donor",
        correction_method=correction_method,
        dtype=dtype,
        max_iter_harmony=20,
        flavor="harmony1",
    )
    x, base = adata_reference.obsm["harmony_org"], adata_reference.obsm["X_pca_harmony"]

    with subtests.test("r"):
        assert _get_measure(x, base, "r").min() > 0.95
    with subtests.test("L2"):
        assert _get_measure(x, base, "L2").max() < 0.1


@pytest.mark.parametrize("correction_method", ["fast", "original"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony2_correction_methods_agree(
    subtests: pytest.Subtests,
    adata_reference: AnnData,
    correction_method: Literal["fast", "original"],
    dtype: DTypeLike,
) -> None:
    """Harmony2 default path: correction methods produce consistent results."""
    harmony_integrate(
        adata_reference,
        "donor",
        correction_method=correction_method,
        dtype=dtype,
        max_iter_harmony=20,
    )
    h2 = adata_reference.obsm["X_pca_harmony"]

    # Run the other method for comparison
    other = "original" if correction_method == "fast" else "fast"
    adata_ref2 = adata_reference.copy()
    harmony_integrate(
        adata_ref2,
        "donor",
        correction_method=other,
        dtype=dtype,
        max_iter_harmony=20,
    )
    h2_ref = adata_ref2.obsm["X_pca_harmony"]

    with subtests.test("r"):
        assert _get_measure(h2, h2_ref, "r").min() > 0.99
    with subtests.test("L2"):
        assert _get_measure(h2, h2_ref, "L2").max() < 0.05


def test_harmony_multiple_keys() -> None:
    """Test Harmony with multiple batch keys."""
    adata = pbmc68k_reduced()
    # Create a second batch key
    adata.obs["batch2"] = np.random.choice(["A", "B", "C"], size=adata.n_obs)

    harmony_integrate(adata, ["bulk_labels", "batch2"], correction_method="original")

    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_custom_parameters() -> None:
    """Test Harmony with custom parameters."""
    adata = pbmc68k_reduced()
    harmony_integrate(
        adata,
        "bulk_labels",
        theta=1.5,
        sigma=0.15,
        n_clusters=50,
        max_iter_harmony=5,
        ridge_lambda=0.5,
        flavor="harmony1",
    )
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_no_nan_output() -> None:
    """Test that Harmony output contains no NaN values."""
    adata = pbmc68k_reduced()
    harmony_integrate(adata, "bulk_labels")
    assert not np.isnan(adata.obsm["X_pca_harmony"]).any()


def test_harmony_input_validation(subtests) -> None:
    """Test that Harmony raises errors for invalid inputs."""
    adata = pbmc68k_reduced()

    with subtests.test("no basis"), pytest.raises(ValueError, match="not available"):
        harmony_integrate(adata, "bulk_labels", basis="nonexistent")
    with subtests.test("no key"), pytest.raises(KeyError):
        harmony_integrate(adata, "nonexistent_key")


def test_harmony_invalid_flavor() -> None:
    """Test that invalid flavor raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="flavor must be"):
        harmony_integrate(adata, "bulk_labels", flavor="harmony3")


@pytest.mark.parametrize("bad_alpha", [-0.1, 0.0, float("inf"), float("nan")])
def test_harmony_integrate_bad_alpha(bad_alpha: float) -> None:
    """Non-positive or non-finite alpha with flavor='harmony2' raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="alpha must be a finite positive"):
        harmony_integrate(adata, "bulk_labels", alpha=bad_alpha)


@pytest.mark.parametrize("bad_threshold", [-0.1, 1.5, 2.0])
def test_harmony_integrate_bad_prune_threshold(bad_threshold: float) -> None:
    """batch_prune_threshold outside [0, 1] raises ValueError."""
    adata = pbmc68k_reduced()
    with pytest.raises(ValueError, match="batch_prune_threshold must be in"):
        harmony_integrate(adata, "bulk_labels", batch_prune_threshold=bad_threshold)


def test_harmony_flavor_warnings() -> None:
    """Test that flavor-incompatible parameter warnings are raised."""
    adata = pbmc68k_reduced()

    # harmony2 with ridge_lambda should warn
    with pytest.warns(UserWarning, match="ridge_lambda is ignored"):
        harmony_integrate(
            adata,
            "bulk_labels",
            flavor="harmony2",
            ridge_lambda=2.0,
            max_iter_harmony=1,
        )

    # harmony1 with alpha should warn
    with pytest.warns(UserWarning, match="alpha is ignored"):
        harmony_integrate(
            adata, "bulk_labels", flavor="harmony1", alpha=0.5, max_iter_harmony=1
        )

    # harmony1 with batch_prune_threshold should warn
    with pytest.warns(UserWarning, match="batch_prune_threshold is ignored"):
        harmony_integrate(
            adata,
            "bulk_labels",
            flavor="harmony1",
            batch_prune_threshold=0.01,
            max_iter_harmony=1,
        )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_pruning(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb suppresses correction for N_b==0 and below-threshold pairs."""
    n_batches, n_clusters = 4, 3
    alpha = 0.2
    threshold = 1e-5
    sentinel = dtype(_SUPPRESS_PENALTY)

    n_b = np.array([0, 100, 1, 50], dtype=dtype)
    o = np.array(
        [[0, 0, 0], [30, 40, 30], [0, 0, 1], [20, 15, 15]],
        dtype=dtype,
    )
    e = np.ones((n_batches, n_clusters), dtype=dtype) * 10

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=alpha,
        threshold=threshold,
        ridge_lambda=1.0,
        dynamic_lambda=True,
    )

    # batch 0 (N_b==0): all clusters must be sentinel
    assert np.all(result[0] == sentinel)
    # batch 1 (well-represented): should be alpha * E = 2.0
    np.testing.assert_allclose(result[1], np.full(n_clusters, alpha * 10, dtype=dtype))
    # batch 2, clusters 0,1 (O/N_b = 0/1 < threshold): sentinel
    assert result[2, 0] == sentinel
    assert result[2, 1] == sentinel
    # batch 2, cluster 2 (O/N_b = 1/1 = 1.0 >= threshold): alpha * E
    np.testing.assert_allclose(result[2, 2], dtype(alpha * 10))
    # batch 3: all alpha * E
    np.testing.assert_allclose(result[3], np.full(n_clusters, alpha * 10, dtype=dtype))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_dynamic_false(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb returns uniform ridge_lambda when dynamic_lambda=False."""
    n_batches, n_clusters = 3, 5
    e = np.ones((n_batches, n_clusters), dtype=dtype)
    o = np.ones((n_batches, n_clusters), dtype=dtype)
    n_b = np.ones(n_batches, dtype=dtype)

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=0.5,
        threshold=1e-5,
        ridge_lambda=1.0,
        dynamic_lambda=False,
    )
    np.testing.assert_array_equal(result, np.full_like(e, 1.0))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_compute_lambda_kb_zero_denom(dtype: type[np.floating]) -> None:
    """_compute_lambda_kb guards against O==0 and E==0 (zero-denominator)."""
    sentinel = dtype(_SUPPRESS_PENALTY)
    e = np.array([[0.0, 5.0]], dtype=dtype)
    o = np.array([[0.0, 10.0]], dtype=dtype)
    n_b = np.array([100.0], dtype=dtype)

    result = _compute_lambda_kb(
        e,
        o=o,
        n_b=n_b,
        alpha=0.2,
        threshold=None,
        ridge_lambda=1.0,
        dynamic_lambda=True,
    )
    # (0,0): O+lambda_kb = 0+0 = 0 -> sentinel
    assert result[0, 0] == sentinel
    # (0,1): normal -> alpha * E = 1.0
    np.testing.assert_allclose(result[0, 1], dtype(1.0))
