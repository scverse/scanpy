from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced


def test_embedding_density():
    # Test that density values are scaled
    # Test that the highest value is in the middle for a grid layout
    test_data = AnnData(X=np.ones((9, 10)))
    test_data.obsm["X_test"] = np.array([[x, y] for x in range(3) for y in range(3)])
    sc.tl.embedding_density(test_data, "test")

    max_dens = np.max(test_data.obs["test_density"])
    min_dens = np.min(test_data.obs["test_density"])
    max_idx = test_data.obs["test_density"].idxmax()

    assert max_idx == "4"
    assert max_dens == 1
    assert min_dens == 0


@pytest.mark.parametrize(
    ("basis", "components"), [("umap", [0, 1]), ("diffmap", [1, 2])]
)
def test_embedding_density_groupby(basis: str, components: list[int]) -> None:
    """Density is estimated per group, on the default components of `basis`."""
    rng = np.random.default_rng(0)
    adata = AnnData(np.ones((60, 1)))
    adata.obs["group"] = pd.Categorical(np.repeat(["a", "b"], 30))
    adata.obsm[f"X_{basis}"] = rng.normal(size=(60, 3))

    sc.tl.embedding_density(adata, basis, groupby="group")

    for idx in adata.obs.groupby("group", observed=True).indices.values():
        subset = AnnData(np.ones((len(idx), 1)))
        subset.obsm["X_test"] = adata.obsm[f"X_{basis}"][idx][:, components]
        sc.tl.embedding_density(subset, "test")
        np.testing.assert_allclose(
            adata.obs[f"{basis}_density_group"].iloc[idx], subset.obs["test_density"]
        )


def test_embedding_density_plot():
    # Test that sc.pl.embedding_density() runs without error
    adata = pbmc68k_reduced()
    sc.tl.embedding_density(adata, "umap")
    sc.pl.embedding_density(adata, "umap", key="umap_density", show=False)
