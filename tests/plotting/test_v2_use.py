"""`use` of the holoviews plots finds what their tools wrote with the same `out`."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc3k, pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

with contextlib.suppress(ImportError):
    from scanpy.plotting import _v2 as pl

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata.acc import AdAcc, LayerAcc, MultiAcc


pytestmark = [needs.scanpy2]


@pytest.fixture
def acc():
    """`hv_anndata`’s accessor, which holoviews understands."""
    return pl.hv_init("bokeh")


def test_highest_expr_genes(acc: AdAcc) -> None:
    adata = pbmc68k_reduced()
    adata.var["symbol"] = adata.var_names.str.lower()
    plot = pl.highest_expr_genes(
        adata, n_top=5, use=acc.layers["counts"], gene_symbols=acc.var["symbol"]
    )
    assert set(plot.data["gene"]) <= set(adata.var["symbol"])


@pytest.mark.usefixtures("acc")  # loads the holoviews plotting extension
def test_highly_variable_genes() -> None:
    adata = pbmc3k()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=100)
    layout = pl.highly_variable_genes(adata)
    [scatter, _] = layout
    assert np.sum(scatter.data["highly_variable"]) == 100


@pytest.mark.parametrize("kind", ["default", "mapping"])
def test_scrublet_score_distribution(acc: AdAcc, kind: str) -> None:
    adata = pbmc3k()[200:400].copy()
    out = {
        "default": False,
        "mapping": {"doublet_score": acc.obs["score"]},
    }[kind]
    sc.pp.scrublet(adata, use_approx_neighbors=False, out=out)
    pl.scrublet_score_distribution(adata, use=None if out is False else out)


def test_embedding_density(acc: AdAcc) -> None:
    adata = pbmc68k_reduced()
    sc.tl.embedding_density(adata, "umap", out="obs.p")
    scatter = pl.embedding_density(adata, acc.obsm["X_umap"], use=acc.obs["p"])
    assert acc.obs["p"] in scatter.vdims


def test_embedding_density_groupby(acc: AdAcc) -> None:
    adata = pbmc68k_reduced()
    sc.tl.embedding_density(adata, "umap", groupby=acc.obs["phase"])
    layout = pl.embedding_density(adata, acc.obsm["X_umap"], groupby=acc.obs["phase"])
    assert len(layout) == adata.obs["phase"].nunique()
    # strings are always resolved, independent of the preset
    assert len(
        pl.embedding_density(adata, acc.obsm["X_umap"], groupby="obs.phase")
    ) == len(layout)


def test_draw_graph(acc: AdAcc) -> None:
    adata = pbmc68k_reduced()
    sc.pp.neighbors(adata, out="obsp.p")
    graph = pl.draw_graph(adata, acc.obsm["X_umap"], acc.obsp["p"])
    assert len(graph.data) == adata.obsp["p"].nnz
    graph = pl.draw_graph(adata, acc.obsm["X_umap"], "obsp.p")
    assert len(graph.data) == adata.obsp["p"].nnz


@pytest.mark.parametrize(
    "use",
    [
        pytest.param(lambda _: None, id="default"),
        pytest.param(lambda acc: acc.layers["counts"], id="accessor"),
        pytest.param(lambda _: "layers.counts", id="str"),
    ],
)
def test_dotplot(
    acc: AdAcc, use: Callable[[AdAcc], LayerAcc | MultiAcc | str | None]
) -> None:
    adata = pbmc68k_reduced()[:, :5]
    points = pl.dotplot(adata, acc.obs["bulk_labels"], use=use(acc))
    expected = sc.get.aggregate(adata, acc.obs["bulk_labels"], "mean", use=use(acc))
    means = expected.to_df("mean")
    np.testing.assert_allclose(
        points.data["mean"],
        [
            means.loc[group, marker]
            for group, marker in zip(
                points.data["group"], points.data["marker"], strict=True
            )
        ],
    )
