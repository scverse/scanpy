from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_raises

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


@pytest.mark.parametrize(
    ("key_added", "key_obsm", "key_uns"),
    [
        pytest.param(None, "X_tsne", "tsne", id="None"),
        pytest.param("custom_key", "custom_key", "custom_key", id="custom_key"),
    ],
)
def test_tsne(key_added: str | None, key_obsm: str, key_uns: str):
    pbmc = pbmc68k_reduced()[:200].copy()

    euclidean1 = sc.tl.tsne(pbmc, metric="euclidean", copy=True)
    with pytest.warns(UserWarning, match="In previous versions of scanpy"):
        euclidean2 = sc.tl.tsne(
            pbmc, metric="euclidean", n_jobs=2, key_added=key_added, copy=True
        )
    cosine = sc.tl.tsne(pbmc, metric="cosine", copy=True)

    # Reproducibility
    np.testing.assert_equal(euclidean1.obsm["X_tsne"], euclidean2.obsm[key_obsm])
    # Metric has some effect
    assert not np.array_equal(euclidean1.obsm["X_tsne"], cosine.obsm["X_tsne"])

    # Params are recorded
    assert euclidean1.uns["tsne"]["params"]["n_jobs"] == 1
    assert euclidean2.uns[key_uns]["params"]["n_jobs"] == 2
    assert cosine.uns["tsne"]["params"]["n_jobs"] == 1
    assert euclidean1.uns["tsne"]["params"]["metric"] == "euclidean"
    assert euclidean2.uns[key_uns]["params"]["metric"] == "euclidean"
    assert cosine.uns["tsne"]["params"]["metric"] == "cosine"


@pytest.mark.parametrize(
    ("key_added", "key_obsm", "key_uns"),
    [
        pytest.param(None, "X_umap", "umap", id="None"),
        pytest.param("custom_key", "custom_key", "custom_key", id="custom_key"),
    ],
)
def test_umap_init_dtype(key_added: str | None, key_obsm: str, key_uns: str):
    pbmc1 = pbmc68k_reduced()[:100, :].copy()
    pbmc2 = pbmc1.copy()
    for pbmc, dtype, k in [(pbmc1, np.float32, None), (pbmc2, np.float64, key_added)]:
        sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(dtype), key_added=k)

    # check that embeddings are close for different dtypes
    assert_array_almost_equal(pbmc1.obsm["X_umap"], pbmc2.obsm[key_obsm])

    # check that params are recorded
    assert pbmc1.uns["umap"]["params"]["a"] == pbmc2.uns[key_uns]["params"]["a"]
    assert pbmc1.uns["umap"]["params"]["b"] == pbmc2.uns[key_uns]["params"]["b"]


@pytest.mark.parametrize(
    "layout",
    [
        pytest.param("fa", marks=needs.fa2),
        pytest.param("fr", marks=needs.igraph),
    ],
)
def test_umap_init_paga(layout):
    pbmc = pbmc68k_reduced()[:100, :].copy()
    sc.tl.paga(pbmc)
    sc.pl.paga(pbmc, layout=layout, show=False)
    sc.tl.umap(pbmc, init_pos="paga")


def test_diffmap():
    pbmc = pbmc68k_reduced()

    sc.tl.diffmap(pbmc)
    d1 = pbmc.obsm["X_diffmap"].copy()
    sc.tl.diffmap(pbmc)
    d2 = pbmc.obsm["X_diffmap"].copy()
    assert_array_equal(d1, d2)

    # Checking if specifying random_state  works, arrays shouldn't be equal
    sc.tl.diffmap(pbmc, random_state=1234)
    d3 = pbmc.obsm["X_diffmap"].copy()
    assert_raises(AssertionError, assert_array_equal, d1, d3)
