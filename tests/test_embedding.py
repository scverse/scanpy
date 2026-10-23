from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_raises

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from anndata import AnnData


@pytest.mark.parametrize(
    ("out", "key_obsm", "key_uns"),
    [
        pytest.param(False, "X_tsne", "tsne", id="default"),
        pytest.param(sc.Prefix("custom"), "custom_tsne", "custom_tsne", id="prefix"),
    ],
)
def test_tsne(*, out: sc.Prefix | bool, key_obsm: str, key_uns: str):
    pbmc = pbmc68k_reduced()[:200].copy()

    euclidean1, euclidean2, cosine = (pbmc.copy() for _ in range(3))
    sc.tl.tsne(euclidean1, metric="euclidean")
    sc.tl.tsne(euclidean2, metric="euclidean", n_jobs=2, out=out)
    sc.tl.tsne(cosine, metric="cosine")

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
    ("out", "key_obsm", "key_uns"),
    [
        pytest.param(False, "X_umap", "umap", id="default"),
        pytest.param(sc.Prefix("custom"), "custom_umap", "custom_umap", id="prefix"),
    ],
)
def test_umap_init_dtype(*, out: sc.Prefix | bool, key_obsm: str, key_uns: str):
    pbmc1 = pbmc68k_reduced()[:100, :].copy()
    pbmc2 = pbmc1.copy()
    for pbmc, dtype, o in [(pbmc1, np.float32, False), (pbmc2, np.float64, out)]:
        sc.tl.umap(pbmc, init_pos=pbmc.obsm["X_pca"][:, :2].astype(dtype), out=o)

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
        sc.tl.diffmap(pbmc, out=True, **{rng_arg: seed}) for seed in (0, 0, 1234)
    )

    with subtests.test("reproducible"):
        assert_array_equal(d1, d2)
    with subtests.test("different embedding"):
        assert_raises(AssertionError, assert_array_equal, d1, d3)


@pytest.mark.parametrize(
    ("out", "key_obsm", "key_uns", "is_dict"),
    [
        pytest.param(False, "X_diffmap", "diffmap_evals", False, id="default"),
        pytest.param(
            *(sc.Prefix("custom"), "custom_diffmap", "custom_diffmap", True),
            id="prefix",
        ),
        pytest.param(sc.Preset.ScanpyV1, "X_diffmap", "diffmap_evals", False, id="v1"),
        pytest.param(
            *(sc.Preset.ScanpyV2Preview, "diffmap", "diffmap", True),
            marks=needs.scanpy2,
            id="v2",
        ),
    ],
)
def test_diffmap_out(
    *,
    out: sc.Prefix | sc.Preset | bool,
    key_obsm: str,
    key_uns: str,
    is_dict: bool,
) -> None:
    adata = pbmc68k_reduced()[:300, :100].copy()
    if isinstance(out, sc.Preset):
        sc.settings.preset = out
        out = False
    sc.tl.diffmap(adata, out=out)
    assert key_obsm in adata.obsm
    assert key_uns in adata.uns
    assert isinstance(adata.uns[key_uns], dict if is_dict else np.ndarray)


@needs.igraph
@pytest.mark.parametrize(
    ("out", "key_obsm", "key_uns"),
    [
        pytest.param(False, "X_draw_graph_fr", "draw_graph", id="default"),
        pytest.param(
            sc.Prefix("custom"), "custom_graph_fr", "custom_graph_fr", id="prefix"
        ),
        pytest.param(sc.Preset.ScanpyV1, "X_draw_graph_fr", "draw_graph", id="v1"),
        pytest.param(
            *(sc.Preset.ScanpyV2Preview, "graph_fr", "graph_fr"),
            marks=needs.scanpy2,
            id="v2",
        ),
    ],
)
def test_draw_graph_out(
    *, out: sc.Prefix | sc.Preset | bool, key_obsm: str, key_uns: str
) -> None:
    adata = pbmc68k_reduced()[:100, :100].copy()
    if isinstance(out, sc.Preset):
        sc.settings.preset = out
        out = False
    sc.tl.draw_graph(adata, layout="fr", out=out)
    assert key_obsm in adata.obsm
    assert key_uns in adata.uns


EMBEDDINGS = [
    pytest.param(sc.tl.tsne, id="tsne"),
    pytest.param(sc.tl.diffmap, id="diffmap"),
    pytest.param(
        lambda a, **kw: sc.tl.draw_graph(a, layout="fr", **kw),
        marks=needs.igraph,
        id="draw_graph",
    ),
]


@pytest.fixture
def pbmc_small() -> AnnData:
    adata = pbmc68k_reduced()[:100, :50].copy()
    sc.pp.neighbors(adata)  # the subset’s graph has isolated cells
    return adata


@needs.anndata_acc
@pytest.mark.parametrize("func", EMBEDDINGS)
@pytest.mark.parametrize("as_str", [True, False], ids=["str", "acc"])
def test_embedding_out_acc(
    pbmc_small: AnnData, func: Callable[..., object], *, as_str: bool
) -> None:
    """An accessor `out` is used literally, and `.uns` follows it."""
    from anndata.acc import A

    obsm_before = set(pbmc_small.obsm)
    func(pbmc_small, out="obsm.m_{layout}" if as_str else A.obsm["m_{layout}"])
    assert set(pbmc_small.obsm) - obsm_before == {"m_{layout}"}
    assert "m_{layout}" in pbmc_small.uns


@pytest.mark.parametrize("func", EMBEDDINGS)
def test_embedding_out_true(pbmc_small: AnnData, func: Callable[..., object]) -> None:
    """`out=True` returns the embedding and writes nothing."""
    before = pbmc_small.copy()
    emb = func(pbmc_small, out=True)
    assert isinstance(emb, np.ndarray)
    assert emb.shape[0] == pbmc_small.n_obs
    assert pbmc_small.obsm.keys() == before.obsm.keys()
    assert pbmc_small.uns.keys() == before.uns.keys()
    with (
        pytest.warns(FutureWarning, match=r"argument copy is deprecated"),
        pytest.raises(ValueError, match=r"`copy=True` cannot be used"),
    ):
        func(pbmc_small, out=True, copy=True)


@needs.anndata_acc
def test_tsne_key_added_deprecated(pbmc_small: AnnData) -> None:
    """`key_added='k'` works like `out=A.obsm['k']`, but warns."""
    from anndata.acc import A

    expected = pbmc_small.copy()
    sc.tl.tsne(expected, out=A.obsm["k"])
    with pytest.warns(FutureWarning, match=r"argument key_added is deprecated"):
        sc.tl.tsne(pbmc_small, key_added="k")
    assert pbmc_small.obsm.keys() == expected.obsm.keys()
    assert pbmc_small.uns.keys() == expected.uns.keys()
