from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_raises

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs


def test_tsne():
    pbmc = pbmc68k_reduced()

    euclidean1 = sc.tl.tsne(pbmc, metric="euclidean", copy=True)
    with pytest.warns(UserWarning, match="In previous versions of scanpy"):
        euclidean2 = sc.tl.tsne(pbmc, metric="euclidean", n_jobs=2, copy=True)
    cosine = sc.tl.tsne(pbmc, metric="cosine", copy=True)

    # Reproducibility
    np.testing.assert_equal(euclidean1.obsm["X_tsne"], euclidean2.obsm["X_tsne"])
    # Metric has some effect
    assert not np.array_equal(euclidean1.obsm["X_tsne"], cosine.obsm["X_tsne"])

    # Params are recorded
    assert euclidean1.uns["tsne"]["params"]["n_jobs"] == 1
    assert euclidean2.uns["tsne"]["params"]["n_jobs"] == 2
    assert cosine.uns["tsne"]["params"]["n_jobs"] == 1
    assert euclidean1.uns["tsne"]["params"]["metric"] == "euclidean"
    assert euclidean2.uns["tsne"]["params"]["metric"] == "euclidean"
    assert cosine.uns["tsne"]["params"]["metric"] == "cosine"


def test_umap_init_dtype():
    pbmc = pbmc68k_reduced()[:100, :].copy()
    sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(np.float32))
    embed1 = pbmc.obsm["X_umap"].copy()
    sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(np.float64))
    embed2 = pbmc.obsm["X_umap"].copy()
    assert_array_almost_equal(embed1, embed2)
    assert_array_almost_equal(embed1, embed2)


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
