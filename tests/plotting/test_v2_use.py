"""`use` of the holoviews plots finds what their tools wrote with the same `out`."""

from __future__ import annotations

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc3k, pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.scanpy2]


@pytest.fixture
def pl():
    from scanpy.plotting import _v2

    return _v2


@pytest.fixture
def acc(pl):
    """`hv_anndata`’s accessor, which holoviews understands."""
    return pl.hv_init("bokeh")


def test_highest_expr_genes(pl, acc) -> None:
    adata = pbmc68k_reduced()
    adata.var["symbol"] = adata.var_names.str.lower()
    plot = pl.highest_expr_genes(
        adata, n_top=5, use=acc.layers["counts"], gene_symbols=acc.var["symbol"]
    )
    assert set(plot.data["gene"]) <= set(adata.var["symbol"])


@pytest.mark.parametrize("kind", ["default", "prefix", "varm"])
def test_highly_variable_genes(pl, acc, kind: str) -> None:
    adata = pbmc3k()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.log1p(adata)
    out = {"default": False, "prefix": sc.Prefix("p"), "varm": acc.varm["hvg"]}[kind]
    sc.pp.highly_variable_genes(adata, n_top_genes=100, out=out)
    layout = pl.highly_variable_genes(adata, use=None if out is False else out)
    [scatter, _] = layout
    assert np.sum(scatter.data["highly_variable"]) == 100


@pytest.mark.parametrize("kind", ["default", "prefix", "mapping"])
def test_scrublet_score_distribution(pl, acc, kind: str) -> None:
    adata = pbmc3k()[200:400].copy()
    out = {
        "default": False,
        "prefix": sc.Prefix("p"),
        "mapping": {"doublet_score": acc.obs["score"]},
    }[kind]
    sc.pp.scrublet(adata, use_approx_neighbors=False, out=out)
    pl.scrublet_score_distribution(adata, use=None if out is False else out)


def test_embedding_density(pl, acc) -> None:
    adata = pbmc68k_reduced()
    sc.tl.embedding_density(adata, basis="umap", out=sc.Prefix("p"))
    scatter = pl.embedding_density(adata, acc.obsm["X_umap"], use=sc.Prefix("p"))
    assert acc.obs["p_umap_density"] in scatter.vdims


def test_draw_graph(pl, acc) -> None:
    adata = pbmc68k_reduced()
    sc.pp.neighbors(adata, out=sc.Prefix("p"))
    graph = pl.draw_graph(adata, acc.obsm["X_umap"], acc.obsp["p_connectivities"])
    assert len(graph.data) == adata.obsp["p_connectivities"].nnz
