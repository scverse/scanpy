from __future__ import annotations

from functools import partial
from importlib.util import find_spec

import numpy as np
import pandas as pd
import pytest
from matplotlib import colormaps
from matplotlib import pyplot as plt
from packaging.version import Version
from scipy import sparse

import scanpy as sc
from scanpy._compat import pkg_version
from testing.scanpy._helpers.data import pbmc3k_processed, pbmc68k_reduced
from testing.scanpy._pytest import needs

pytestmark = [needs.igraph]


SKIP_IF_OLD_IGRAPH = pytest.mark.skipif(
    not find_spec("igraph") or pkg_version("igraph") < Version("1"),
    reason="igraph 0.x has different RNG behavior",
)


@pytest.fixture(scope="module", params=["dense", "sparse"])
def pbmc_session(request: pytest.FixtureRequest) -> sc.AnnData:
    pbmc = pbmc68k_reduced()
    pbmc.obs["cool_feature"] = pbmc[:, "CST3"].X.squeeze().copy()
    if request.param == "sparse":
        pbmc.X = pbmc.raw.X.tocsr()
    sc.tl.paga(pbmc, groups="bulk_labels")
    assert not pbmc.obs["cool_feature"].isna().all()
    return pbmc


@pytest.fixture
def pbmc(pbmc_session):
    return pbmc_session.copy()


@SKIP_IF_OLD_IGRAPH
@pytest.mark.parametrize(
    ("test_id", "func"),
    [
        ("", sc.pl.paga),
        ("continuous", partial(sc.pl.paga, color="CST3")),
        ("continuous_obs", partial(sc.pl.paga, color="cool_feature")),
        ("continuous_multiple", partial(sc.pl.paga, color=["CST3", "GATA2"])),
        ("compare", partial(sc.pl.paga_compare, legend_fontoutline=2)),
        pytest.param(
            "compare_continuous",
            partial(sc.pl.paga_compare, color="CST3", legend_fontsize=5),
            marks=pytest.mark.xfail(reason="expects .uns['paga']['pos']"),
        ),
        (
            "compare_pca",
            partial(sc.pl.paga_compare, basis="X_pca", legend_fontweight="normal"),
        ),
    ],
)
def test_paga_plots(plot_cmp, pbmc, test_id, func):
    common = dict(threshold=0.5, max_edge_width=1.0, random_state=0, show=False)

    func(pbmc, **common)
    plot_cmp(f"paga_{test_id}" if test_id else "paga")


@SKIP_IF_OLD_IGRAPH
def test_paga_pie(plot_cmp, pbmc) -> None:
    colors = {
        c: {colormaps["Set1"](i): 0.33 for i in range(3)}
        for c in pbmc.obs["bulk_labels"].cat.categories
    }
    colors["Dendritic"] = {colormaps["Set2"](i): 0.25 for i in range(4)}

    sc.pl.paga(pbmc, color=colors, colorbar=False, show=False)
    plot_cmp("paga_pie")


def test_paga_path(plot_cmp, pbmc) -> None:
    pbmc.uns["iroot"] = 0
    sc.tl.dpt(pbmc)
    sc.pl.paga_path(
        pbmc,
        nodes=["Dendritic"],
        keys=["HES4", "SRM", "CSTB"],
        show=False,
    )
    plot_cmp("paga_path")


def test_paga_compare(plot_cmp):
    # Tests that https://github.com/scverse/scanpy/issues/1887 is fixed
    pbmc = pbmc3k_processed()
    sc.tl.paga(pbmc, groups="louvain")

    sc.pl.paga_compare(pbmc, basis="umap", show=False)

    plot_cmp("paga_compare_pbmc3k")


def test_paga_ncols() -> None:
    # Tests that https://github.com/scverse/scanpy/issues/1203 is fixed
    rng = np.random.default_rng(0)
    adata = sc.AnnData(rng.random((80, 20)))
    adata.obs["group"] = pd.Categorical(rng.choice(["a", "b", "c", "d", "e"], 80))
    for i in range(4):
        adata.obs[f"c{i}"] = pd.Categorical(rng.choice(["x", "y"], 80))

    k = 5
    rows = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 0])
    cols = np.array([1, 0, 2, 1, 3, 2, 4, 3, 0, 4])
    connectivities = sparse.csr_matrix(  # noqa: TID251
        (np.ones(len(rows)), (rows, cols)), shape=(k, k)
    )
    adata.uns["paga"] = {
        "groups": "group",
        "connectivities": connectivities,
        "connectivities_tree": connectivities.copy(),
    }
    pos = rng.random((k, 2))
    colors = ["c0", "c1", "c2", "c3"]

    # `ncols` wraps the panels into a grid
    axs = sc.pl.paga(adata, color=colors, ncols=2, pos=pos, show=False)
    assert len(axs) == 4
    gridspec = axs[0].get_subplotspec().get_gridspec()
    assert (gridspec.nrows, gridspec.ncols) == (2, 2)

    # the default layout keeps all panels in a single row
    axs = sc.pl.paga(adata, color=colors, pos=pos, show=False)
    assert len(axs) == 4

    # continuous colors (with colorbars) also wrap into a grid
    gene_colors = adata.var_names[:3].tolist()
    axs = sc.pl.paga(adata, color=gene_colors, ncols=2, pos=pos, show=False)
    assert len(axs) == 3
    gridspec = axs[0].get_subplotspec().get_gridspec()
    assert gridspec.ncols == 2

    # `ncols` cannot be combined with a pre-supplied `ax`
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="`ncols` cannot be combined"):
        sc.pl.paga(adata, color=colors, ncols=2, ax=ax, pos=pos, show=False)
