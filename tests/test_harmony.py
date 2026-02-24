from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pooch
import pytest
from anndata import AnnData
from scipy.stats import pearsonr

from scanpy.preprocessing import harmony_integrate
from testing.scanpy._helpers.data import pbmc68k_reduced

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike


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


@pytest.fixture
def adata_reference() -> AnnData:
    """Load reference data from harmonypy repository."""
    paths = {
        f: pooch.retrieve(
            f"https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/{name}",
            known_hash=hash_,
        )
        for f, (name, hash_) in DATA.items()
    }
    dfs = {f: pd.read_csv(path, delimiter="\t") for f, path in paths.items()}
    # Create unique index using row number + cell name
    dfs["meta"].index = [f"{i}_{cell}" for i, cell in enumerate(dfs["meta"]["cell"])]
    adata = AnnData(
        X=None,
        obs=dfs["meta"],
        obsm={"X_pca": dfs["pca"].values, "harmony_org": dfs["pca_harmonized"].values},
    )
    return adata


@pytest.mark.parametrize("correction_method", ["fast", "original"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate(
    correction_method: Literal["fast", "original"], dtype: DTypeLike
) -> None:
    """Test that Harmony integrate works."""
    adata = pbmc68k_reduced()
    harmony_integrate(
        adata, "bulk_labels", correction_method=correction_method, dtype=dtype
    )
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate_algos(subtests: pytest.Subtests, dtype: DTypeLike) -> None:
    """Test that both correction methods produce similar results."""
    adata = pbmc68k_reduced()

    harmony_integrate(adata, "bulk_labels", correction_method="fast", dtype=dtype)
    fast = adata.obsm["X_pca_harmony"].copy()
    harmony_integrate(adata, "bulk_labels", correction_method="original", dtype=dtype)
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
    """Test that Harmony produces results similar to the reference implementation."""
    harmony_integrate(
        adata_reference,
        "donor",
        correction_method=correction_method,
        dtype=dtype,
        max_iter_harmony=20,
    )
    x, base = adata_reference.obsm["harmony_org"], adata_reference.obsm["X_pca_harmony"]

    with subtests.test("r"):
        assert _get_measure(x, base, "r").min() > 0.95
    with subtests.test("L2"):
        assert _get_measure(x, base, "L2").max() < 0.05


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
