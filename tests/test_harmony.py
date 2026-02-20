from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pooch
import pytest
from anndata import AnnData
from scipy.stats import pearsonr

import scanpy as sc
from scanpy.preprocessing import harmony_integrate

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import DTypeLike


def _get_measure(
    x: np.ndarray, base: np.ndarray, norm: Literal["r", "L2"]
) -> np.ndarray:
    """Compute correlation or L2 distance between arrays."""
    if norm == "r":
        # Compute per-column correlation
        if x.ndim == 1:
            corr, _ = pearsonr(x, base)
            return corr
        corrs = []
        for i in range(x.shape[1]):
            corr, _ = pearsonr(x[:, i], base[:, i])
            corrs.append(corr)
        return np.array(corrs)

    assert norm == "L2"
    # L2 distance normalized by base norm
    if x.ndim == 1:
        return np.linalg.norm(x - base) / np.linalg.norm(base)
    dists = []
    for i in range(x.shape[1]):
        dist = np.linalg.norm(x[:, i] - base[:, i]) / np.linalg.norm(base[:, i])
        dists.append(dist)
    return np.array(dists)


@pytest.fixture
def adata_reference() -> AnnData:
    """Load reference data from harmonypy repository."""
    x_pca_file = pooch.retrieve(
        "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_pcs.tsv.gz",
        known_hash="md5:27e319b3ddcc0c00d98e70aa8e677b10",
    )
    x_pca = pd.read_csv(x_pca_file, delimiter="\t")
    x_pca_harmony_file = pooch.retrieve(
        "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_pcs_harmonized.tsv.gz",
        known_hash="md5:a7c4ce4b98c390997c66d63d48e09221",
    )
    x_pca_harmony = pd.read_csv(x_pca_harmony_file, delimiter="\t")
    meta_file = pooch.retrieve(
        "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_meta.tsv.gz",
        known_hash="md5:8c7ca20e926513da7cf0def1211baecb",
    )
    meta = pd.read_csv(meta_file, delimiter="\t")
    # Create unique index using row number + cell name
    meta.index = [f"{i}_{cell}" for i, cell in enumerate(meta["cell"])]
    adata = AnnData(
        X=None,
        obs=meta,
        obsm={"X_pca": x_pca.values, "harmony_org": x_pca_harmony.values},
    )
    return adata


@pytest.mark.parametrize("correction_method", ["fast", "original"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate(
    correction_method: Literal["fast", "original"], dtype: DTypeLike
) -> None:
    """Test that Harmony integrate works."""
    adata = sc.datasets.pbmc68k_reduced()

    harmony_integrate(
        adata, "bulk_labels", correction_method=correction_method, dtype=dtype
    )

    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_harmony_integrate_algos(subtests: pytest.Subtests, dtype: DTypeLike) -> None:
    """Test that both correction methods produce similar results."""
    adata = sc.datasets.pbmc68k_reduced()

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
    adata = sc.datasets.pbmc68k_reduced()
    # Create a second batch key
    adata.obs["batch2"] = np.random.choice(["A", "B", "C"], size=adata.n_obs)

    harmony_integrate(adata, ["bulk_labels", "batch2"], correction_method="original")

    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_harmony_custom_parameters() -> None:
    """Test Harmony with custom parameters."""
    adata = sc.datasets.pbmc68k_reduced()
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
    adata = sc.datasets.pbmc68k_reduced()
    harmony_integrate(adata, "bulk_labels")
    assert not np.isnan(adata.obsm["X_pca_harmony"]).any()


def test_harmony_input_validation(subtests) -> None:
    """Test that Harmony raises errors for invalid inputs."""
    adata = sc.datasets.pbmc68k_reduced()

    with subtests.test("no basis"), pytest.raises(ValueError, match="not available"):
        harmony_integrate(adata, "bulk_labels", basis="nonexistent")
    with subtests.test("no key"), pytest.raises(KeyError):
        harmony_integrate(adata, "nonexistent_key")
