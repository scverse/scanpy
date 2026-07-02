from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_raises

import scanpy as sc
from scanpy._settings import Default
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from typing import Literal


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
    euclidean2 = sc.tl.tsne(
        pbmc, metric="euclidean", n_jobs=2, key_added=key_added, copy=True
    )
    cosine = sc.tl.tsne(pbmc, metric="cosine", copy=True)

    # Reproducibility
    np.testing.assert_equal(euclidean1.obsm["X_tsne"], euclidean2.obsm[key_obsm])
    # Metric has some effect
    assert not np.array_equal(euclidean1.obsm["X_tsne"], cosine.obsm["X_tsne"])

    # Params are recorded
    assert euclidean1.uns["tsne"]["params"]["n_jobs"] == sc.settings.n_jobs
    assert euclidean2.uns[key_uns]["params"]["n_jobs"] == 2
    assert cosine.uns["tsne"]["params"]["n_jobs"] == sc.settings.n_jobs
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


def test_umap_preserves_connectivities():
    # https://github.com/scverse/scanpy/issues/4028
    pbmc = pbmc68k_reduced()[:100, :].copy()
    conn = pbmc.obsp["connectivities"]
    data_before = conn.data.copy()
    nnz_before = conn.nnz

    sc.tl.umap(pbmc)

    assert_array_equal(data_before, conn.data)
    assert conn.nnz == nnz_before
    assert (conn.data == 0).sum() == 0, "CSR should have no explicit zeros"
    assert "X_umap" in pbmc.obsm


@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_diffmap(
    subtests: pytest.Subtests, rng_arg: Literal["rng", "random_state"]
) -> None:
    pbmc = pbmc68k_reduced()

    d1, d2, d3 = (
        sc.tl.diffmap(pbmc, copy=True, **{rng_arg: seed}).obsm["X_diffmap"].copy()
        for seed in (0, 0, 1234)
    )

    with subtests.test("reproducible"):
        assert_array_equal(d1, d2)
    with subtests.test("different embedding"):
        assert_raises(AssertionError, assert_array_equal, d1, d3)


@pytest.mark.parametrize(
    ("key_added", "key_obsm", "key_uns"),
    [
        pytest.param(None, "X_diffmap", "diffmap_evals", id="None"),
        pytest.param("custom_key", "custom_key", "custom_key_evals", id="custom_key"),
        pytest.param(sc.Preset.ScanpyV1, "X_diffmap", "diffmap_evals", id="v1"),
        pytest.param(
            *(sc.Preset.ScanpyV2Preview, "diffmap", "diffmap_evals"),
            marks=[needs.igraph, needs.skmisc],
            id="v2",
        ),
    ],
)
def test_diffmap_key_added(
    key_added: str | None | Default | sc.Preset, key_obsm: str, key_uns: str
) -> None:
    pbmc = pbmc68k_reduced()[:300, :100].copy()
    if isinstance(key_added, sc.Preset):
        sc.settings.preset = key_added
        key_added = Default()
    adata = sc.tl.diffmap(pbmc, key_added=key_added, copy=True)
    assert key_obsm in adata.obsm
    assert key_uns in adata.uns


@needs.igraph
@pytest.mark.parametrize(
    ("key_added", "key_obsm", "key_uns"),
    [
        pytest.param(None, "X_draw_graph_fr", "draw_graph", id="None"),
        pytest.param("custom_{layout}", "custom_fr", "custom_fr", id="custom_template"),
        pytest.param(sc.Preset.ScanpyV1, "X_draw_graph_fr", "draw_graph", id="v1"),
        pytest.param(
            *(sc.Preset.ScanpyV2Preview, "graph_fr", "graph_fr"),
            marks=needs.skmisc,
            id="v2",
        ),
    ],
)
def test_draw_graph_key_added(
    key_added: str | None | Default | sc.Preset, key_obsm: str, key_uns: str
) -> None:
    pbmc = pbmc68k_reduced()[:100, :100].copy()
    if isinstance(key_added, sc.Preset):
        sc.settings.preset = key_added
        key_added = Default()
    adata = sc.tl.draw_graph(pbmc, layout="fr", key_added=key_added, copy=True)
    assert key_obsm in adata.obsm
    assert key_uns in adata.uns
