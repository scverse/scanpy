from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata import AnnData
    from anndata.acc import AdAcc


pytestmark = [needs.scanpy2]


def test_supported_opts_keeps_dim_expressions() -> None:
    import holoviews as hv

    from scanpy.plotting._v2._core import _supported_opts

    hv.extension("bokeh")
    size = hv.dim("x") * 2
    opts = _supported_opts(hv.Points, size=size, s=size)  # "s" is matplotlib-only

    assert opts.options == {"size": size}


@pytest.fixture(scope="module")
def v2_adata() -> AnnData:
    """Tiny AnnData, built once per module (nothing below mutates it)."""
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(
        np.arange(12, dtype=np.float32).reshape(4, 3),
        obs=pd.DataFrame(
            dict(
                grp=pd.Categorical(list("aabb")),
                n=[1.0, 2.0, 3.0, 4.0],
                pca_density=[0.1, 0.2, 0.3, 0.4],
            ),
            index=[f"c{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
        obsm=dict(X_pca=np.arange(8, dtype=np.float32).reshape(4, 2)),
    )
    adata.layers["counts"] = adata.X
    sc.pp.neighbors(adata, n_neighbors=3, use_rep="X_pca")
    return adata


@pytest.fixture
def v2_acc() -> AdAcc:
    """Switch to the v2 preset (reset per test by `original_settings`) and get `A`."""
    sc.settings.preset = sc.Preset.ScanpyV2Preview
    return sc.pl.hv_init("bokeh")


@pytest.mark.parametrize(
    ("plot", "expected"),
    [
        pytest.param(
            lambda adata, acc: sc.pl.violin(adata, "obs.n", color="obs.grp"),
            "Violin",
            id="violin",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.violin(adata, ["obs.n", "X[:,g0]"]),
            "Layout",
            id="violin-multi",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.heatmap(adata, "layers.counts", ["obs.n"]),
            "HeatMap",
            id="heatmap",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.tracksplot(adata, ["X[:,g0]"], color="obs.grp"),
            "NdLayout",
            id="tracksplot",
        ),
        # `acc.var.index` has no string form: `"var.index"` means the column named “index”
        pytest.param(
            lambda adata, acc: sc.pl.stacked_violin(adata, acc.var.index, "obs.grp"),
            "GridSpace",
            id="stacked_violin",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.dotplot(adata, "obs.grp"), "Points", id="dotplot"
        ),
        pytest.param(
            lambda adata, acc: sc.pl.matrixplot(adata, "obs.grp", data="layers.counts"),
            "HeatMap",
            id="matrixplot",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.ranking(adata, "obs.n", 2), "Overlay", id="ranking"
        ),
        pytest.param(
            lambda adata, acc: sc.pl.embedding_density(adata, "obsm.X_pca"),
            "Scatter",
            id="embedding_density",
        ),
        pytest.param(
            lambda adata, acc: sc.pl.draw_graph(
                adata, "obsm.X_pca", node_vdims="obs.n"
            ),
            "Graph",
            id="draw_graph",
        ),
    ],
)
def test_plots_accept_str_specs(
    v2_adata: AnnData,
    v2_acc: AdAcc,
    plot: Callable[[AnnData, AdAcc], object],
    expected: str,
) -> None:
    """Every v2 plot taking dims/accessors should also take string specs."""
    assert type(plot(v2_adata, v2_acc)).__name__ == expected


def test_str_specs_match_refs(v2_adata: AnnData, v2_acc: AdAcc) -> None:
    adata, acc = v2_adata, v2_acc
    from_str = sc.pl.scatter(
        adata, ["obsm.X_pca.0", acc.obsm["X_pca"][1]], color="obs.grp"
    )
    from_ref = sc.pl.scatter(adata, acc.obsm["X_pca"][:, [0, 1]], color=acc.obs["grp"])
    assert from_str.dimensions() == from_ref.dimensions()
