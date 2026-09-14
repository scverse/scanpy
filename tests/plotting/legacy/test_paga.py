from __future__ import annotations

from functools import partial
from importlib.util import find_spec

import pytest
from matplotlib import colormaps
from packaging.version import Version

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
