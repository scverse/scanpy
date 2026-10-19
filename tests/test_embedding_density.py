from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


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


@pytest.fixture
def adata_rand() -> AnnData:
    adata = AnnData(X=np.ones((40, 1)))
    adata.obs["group"] = np.repeat(["a", "b"], 20)  # str, not yet categorical
    adata.obsm["X_test"] = np.random.default_rng(0).normal(size=(40, 2))
    return adata


@needs.anndata_acc
@pytest.mark.parametrize("as_str", [True, False], ids=["str", "acc"])
def test_embedding_density_out_acc(adata_rand: AnnData, *, as_str: bool) -> None:
    from anndata.acc import A

    sc.tl.embedding_density(adata_rand, "test", out="obs.d" if as_str else A.obs["d"])
    assert "d" in adata_rand.obs
    assert "test_density" not in adata_rand.obs
    assert "d_params" in adata_rand.uns


def test_embedding_density_prefix(adata_rand: AnnData) -> None:
    sc.tl.embedding_density(adata_rand, "test", groupby="group", out=sc.Prefix("p"))
    assert "p_test_density_group" in adata_rand.obs
    assert "p_test_density_group_params" in adata_rand.uns


def test_embedding_density_out_true(adata_rand: AnnData) -> None:
    """`out=True` returns the densities and writes nothing (not even categoricals)."""
    before = adata_rand.copy()
    dens = sc.tl.embedding_density(adata_rand, "test", groupby="group", out=True)
    pd.testing.assert_frame_equal(adata_rand.obs, before.obs)
    assert not adata_rand.uns
    sc.tl.embedding_density(adata_rand, "test", groupby="group")
    np.testing.assert_array_equal(dens, adata_rand.obs["test_density_group"])


@needs.anndata_acc
def test_embedding_density_key_added_deprecated(adata_rand: AnnData) -> None:
    """`key_added='k'` works like `out=A.obs['k']`, but warns."""
    from anndata.acc import A

    expected = adata_rand.copy()
    sc.tl.embedding_density(expected, "test", out=A.obs["k"])
    with pytest.warns(FutureWarning, match=r"argument key_added is deprecated"):
        sc.tl.embedding_density(adata_rand, "test", key_added="k")
    pd.testing.assert_frame_equal(adata_rand.obs, expected.obs)
    assert adata_rand.uns.keys() == expected.uns.keys() == {"k_params"}
