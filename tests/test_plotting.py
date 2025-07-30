from __future__ import annotations

from functools import partial
from itertools import chain, combinations, repeat
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns
from anndata import AnnData
from matplotlib.testing.compare import compare_images
from packaging.version import Version

import scanpy as sc
from scanpy._compat import pkg_version
from testing.scanpy._helpers.data import (
    krumsiek11,
    pbmc3k,
    pbmc3k_processed,
    pbmc68k_reduced,
)
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from matplotlib.axes import Axes


HERE: Path = Path(__file__).parent
ROOT = HERE / "_images"


# Test images are saved in the directory ./_images/<test-name>/
# If test images need to be updated, simply copy actual.png to expected.png.


@pytest.mark.parametrize("col", [None, "symb"])
@pytest.mark.parametrize("layer", [None, "layer_name"])
def test_highest_expr_genes(image_comparer, col, layer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=5)

    adata = pbmc3k()
    if layer is not None:
        adata.layers[layer] = adata.X
        del adata.X
    # check that only existing categories are shown
    adata.var["symb"] = adata.var_names.astype("category")

    sc.pl.highest_expr_genes(adata, 20, gene_symbols=col, layer=layer, show=False)

    save_and_compare_images("highest_expr_genes")


@needs.leidenalg
@pytest.mark.parametrize(
    ("params", "key"),
    [
        pytest.param({}, "heatmap", id="default"),
        pytest.param(
            dict(swap_axes=True, figsize=(10, 3), cmap="YlGnBu"),
            "heatmap_swap_axes",
            id="swap",
        ),
        pytest.param(
            dict(
                groupby="numeric_value",
                num_categories=4,
                figsize=(4.5, 5),
                dendrogram=False,
            ),
            "heatmap2",
            id="numeric",
        ),
        pytest.param(
            dict(standard_scale="var", layer="test"),
            "heatmap_std_scale_var",
            id="std_scale=var",
        ),
        pytest.param(
            dict(standard_scale="obs"),
            "heatmap_std_scale_obs",
            id="std_scale=obs",
        ),
    ],
)
def test_heatmap(image_comparer, params: dict[str, Any], key: str) -> None:
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = krumsiek11()
    adata.obs["numeric_value"] = adata.X[:, 0]
    adata.layers["test"] = -1 * adata.X.copy()

    params = dict(groupby="cell_type", dendrogram=True) | params
    sc.pl.heatmap(adata, adata.var_names, **params, use_raw=False, show=False)
    save_and_compare_images(key)


@needs.leidenalg
def test_heatmap_var_as_dict(image_comparer) -> None:
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.tl.leiden(
        pbmc,
        key_added="clusters",
        resolution=0.5,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    # call umap to trigger colors for the clusters
    sc.pl.umap(pbmc, color="clusters")
    marker_genes_dict = {
        "3": ["GNLY", "NKG7"],
        "1": ["FCER1A"],
        "2": ["CD3D"],
        "0": ["FCGR3A"],
        "4": ["CD79A", "MS4A1"],
    }
    sc.pl.heatmap(
        adata=pbmc,
        var_names=marker_genes_dict,
        groupby="clusters",
        vmin=-2,
        vmax=2,
        cmap="RdBu_r",
        dendrogram=True,
        swap_axes=True,
    )
    save_and_compare_images("heatmap_var_as_dict")


@needs.leidenalg
@pytest.mark.parametrize("swap_axes", [True, False])
def test_heatmap_alignment(*, image_comparer, swap_axes: bool) -> None:
    """Test that plot elements are well aligned."""
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    a = AnnData(
        np.array([[0, 0.3, 0.5], [1, 1.3, 1.5], [2, 2.3, 2.5]]),
        obs={"foo": ["a", "b", "c"]},
        var=pd.DataFrame({"genes": ["g1", "g2", "g3"]}).set_index("genes"),
    )
    a.obs["foo"] = a.obs["foo"].astype("category")
    sc.pl.heatmap(
        a, var_names=a.var_names, groupby="foo", swap_axes=swap_axes, figsize=(4, 4)
    )
    save_and_compare_images(f"heatmap_small{'_swap' if swap_axes else ''}_alignment")


@pytest.mark.skipif(
    pkg_version("matplotlib") < Version("3.1"),
    reason="https://github.com/mwaskom/seaborn/issues/1953",
)
@pytest.mark.parametrize(
    ("obs_keys", "name"),
    [(None, "clustermap"), ("cell_type", "clustermap_withcolor")],
)
def test_clustermap(image_comparer, obs_keys, name):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = krumsiek11()
    sc.pl.clustermap(adata, obs_keys)
    save_and_compare_images(name)


params_dotplot_matrixplot_stacked_violin = [
    pytest.param(id, fn, id=id)
    for id, fn in [
        (
            "dotplot",
            partial(
                sc.pl.dotplot, groupby="cell_type", title="dotplot", dendrogram=True
            ),
        ),
        (
            "dotplot2",
            partial(
                sc.pl.dotplot,
                groupby="numeric_column",
                use_raw=False,
                num_categories=7,
                title="non categorical obs",
                figsize=(7, 2.5),
            ),
        ),
        (
            "dotplot3",
            partial(
                sc.pl.dotplot,
                groupby="cell_type",
                dot_max=0.7,
                dot_min=0.1,
                cmap="hot_r",
                title="dot_max=0.7 dot_min=0.1, var_groups",
                var_group_positions=[(0, 1), (9, 10)],
                var_group_labels=["A", "B"],
                dendrogram=True,
            ),
        ),
        (
            "dotplot_std_scale_group",
            partial(
                sc.pl.dotplot,
                groupby="cell_type",
                use_raw=False,
                dendrogram=True,
                layer="test",
                swap_axes=True,
                title="swap_axes, layer=-1*X, scale=group\nsmallest_dot=10",
                standard_scale="group",
                smallest_dot=10,
            ),
        ),
        (
            "dotplot_dict",
            partial(
                sc.pl.dotplot,
                groupby="cell_type",
                dot_max=0.7,
                dot_min=0.1,
                color_map="winter",
                title="var as dict",
                dendrogram=True,
            ),
        ),
        (
            "matrixplot",
            partial(
                sc.pl.matrixplot,
                groupby="cell_type",
                use_raw=False,
                title="matrixplot",
                dendrogram=True,
            ),
        ),
        (
            "matrixplot_std_scale_var_dict",
            partial(
                sc.pl.matrixplot,
                groupby="cell_type",
                dendrogram=True,
                standard_scale="var",
                layer="test",
                cmap="Blues_r",
                title='scale var, custom colorbar_title, layer="test"',
                colorbar_title="Scaled expression",
            ),
        ),
        (
            "matrixplot_std_scale_group",
            partial(
                sc.pl.matrixplot,
                groupby="cell_type",
                use_raw=False,
                standard_scale="group",
                title="scale_group, swap_axes",
                swap_axes=True,
            ),
        ),
        (
            "matrixplot2",
            partial(
                sc.pl.matrixplot,
                groupby="numeric_column",
                use_raw=False,
                num_categories=4,
                title="non-categorical obs, custom figsize",
                figsize=(8, 2.5),
                cmap="RdBu_r",
            ),
        ),
        (
            "stacked_violin",
            partial(
                sc.pl.stacked_violin,
                groupby="cell_type",
                use_raw=False,
                title="stacked_violin",
                dendrogram=True,
            ),
        ),
        (
            "stacked_violin_std_scale_var_dict",
            partial(
                sc.pl.stacked_violin,
                groupby="cell_type",
                dendrogram=True,
                standard_scale="var",
                layer="test",
                title='scale var, layer="test"',
            ),
        ),
        (
            "stacked_violin_std_scale_group",
            partial(
                sc.pl.stacked_violin,
                groupby="cell_type",
                use_raw=False,
                standard_scale="group",
                title="scale_group\nswap_axes",
                swap_axes=True,
                cmap="Blues",
            ),
        ),
        (
            "stacked_violin_no_cat_obs",
            partial(
                sc.pl.stacked_violin,
                groupby="numeric_column",
                use_raw=False,
                num_categories=4,
                title="non-categorical obs, custom figsize",
                figsize=(8, 2.5),
            ),
        ),
    ]
]


@pytest.mark.parametrize(("id", "fn"), params_dotplot_matrixplot_stacked_violin)
def test_dotplot_matrixplot_stacked_violin(image_comparer, id, fn):
    save_and_compare_images = partial(image_comparer, ROOT, tol=5)

    adata = krumsiek11()
    adata.obs["numeric_column"] = adata.X[:, 0]
    adata.layers["test"] = -1 * adata.X.copy()
    genes_dict = {
        "group a": ["Gata2", "Gata1"],
        "group b": ["Fog1", "EKLF", "Fli1", "SCL"],
        "group c": ["Cebpa", "Pu.1", "cJun", "EgrNab", "Gfi1"],
    }

    if id.endswith("dict"):
        fn(adata, genes_dict, show=False)
    else:
        fn(adata, adata.var_names, show=False)
    save_and_compare_images(id)


def test_dotplot_obj(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    # test dotplot dot_min, dot_max, color_map, and var_groups
    pbmc = pbmc68k_reduced()
    genes = [
        *["CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3"],
        *["S100A8", "GNLY", "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3"],
    ]
    # test layer, var standardization, smallest_dot,
    # color title, size_title return_fig and dot_edge
    pbmc.layers["test"] = pbmc.X * -1
    plot = sc.pl.dotplot(
        pbmc,
        genes,
        "bulk_labels",
        layer="test",
        dendrogram=True,
        return_fig=True,
        standard_scale="var",
        smallest_dot=40,
        colorbar_title="scaled column max",
        size_title="Fraction of cells",
    )
    plot.style(dot_edge_color="black", dot_edge_lw=0.1, cmap="Reds").show()

    save_and_compare_images("dotplot_std_scale_var")


def test_dotplot_style_no_reset():
    pbmc = pbmc68k_reduced()
    plot = sc.pl.dotplot(pbmc, "CD79A", "bulk_labels", return_fig=True)
    assert isinstance(plot, sc.pl.DotPlot)
    assert plot.cmap == sc.pl.DotPlot.DEFAULT_COLORMAP
    plot.style(cmap="winter")
    assert plot.cmap == "winter"
    plot.style(color_on="square")
    assert plot.cmap == "winter", "style() should not reset unspecified parameters"


def test_dotplot_add_totals(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=5)

    pbmc = pbmc68k_reduced()
    markers = {"T-cell": "CD3D", "B-cell": "CD79A", "myeloid": "CST3"}
    sc.pl.dotplot(pbmc, markers, "bulk_labels", return_fig=True).add_totals().show()
    save_and_compare_images("dotplot_totals")


def test_matrixplot_obj(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = pbmc68k_reduced()
    marker_genes_dict = {
        "3": ["GNLY", "NKG7"],
        "1": ["FCER1A"],
        "2": ["CD3D"],
        "0": ["FCGR3A"],
        "4": ["CD79A", "MS4A1"],
    }

    plot = sc.pl.matrixplot(
        adata,
        marker_genes_dict,
        "bulk_labels",
        use_raw=False,
        title="added totals",
        return_fig=True,
    )
    plot.add_totals(sort="descending").style(edge_color="white", edge_lw=0.5).show()
    save_and_compare_images("matrixplot_with_totals")

    axes = plot.get_axes()
    assert "mainplot_ax" in axes, "mainplot_ax not found in returned axes dict"


def test_stacked_violin_obj(image_comparer, plt):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    markers = {
        "T-cell": ["CD3D", "CD3E", "IL32"],
        "B-cell": ["CD79A", "CD79B", "MS4A1"],
        "myeloid": ["CST3", "LYZ"],
    }
    plot = sc.pl.stacked_violin(
        pbmc,
        markers,
        "bulk_labels",
        use_raw=False,
        title="return_fig. add_totals",
        return_fig=True,
    )
    plot.add_totals().style(row_palette="tab20").show()
    save_and_compare_images("stacked_violin_return_fig")


# checking for https://github.com/scverse/scanpy/issues/3152
def test_stacked_violin_swap_axes_match(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=10)
    pbmc = pbmc68k_reduced()
    sc.tl.rank_genes_groups(
        pbmc,
        "bulk_labels",
        method="wilcoxon",
        tie_correct=True,
        pts=True,
        key_added="wilcoxon",
    )
    swapped_ax = sc.pl.rank_genes_groups_stacked_violin(
        pbmc,
        n_genes=2,
        key="wilcoxon",
        groupby="bulk_labels",
        swap_axes=True,
        return_fig=True,
    )
    swapped_ax.show()
    save_and_compare_images("stacked_violin_swap_axes_pbmc68k_reduced")


def test_tracksplot(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = krumsiek11()
    sc.pl.tracksplot(
        adata, adata.var_names, "cell_type", dendrogram=True, use_raw=False
    )
    save_and_compare_images("tracksplot")


def test_multiple_plots(image_comparer):
    # only testing stacked_violin, matrixplot and dotplot
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = pbmc68k_reduced()
    markers = {
        "T-cell": ["CD3D", "CD3E", "IL32"],
        "B-cell": ["CD79A", "CD79B", "MS4A1"],
        "myeloid": ["CST3", "LYZ"],
    }
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(20, 5), gridspec_kw={"wspace": 0.7}
    )
    _ = sc.pl.stacked_violin(
        adata,
        markers,
        groupby="bulk_labels",
        ax=ax1,
        title="stacked_violin",
        dendrogram=True,
        show=False,
    )
    _ = sc.pl.dotplot(
        adata,
        markers,
        groupby="bulk_labels",
        ax=ax2,
        title="dotplot",
        dendrogram=True,
        show=False,
    )
    _ = sc.pl.matrixplot(
        adata,
        markers,
        groupby="bulk_labels",
        ax=ax3,
        title="matrixplot",
        dendrogram=True,
        show=False,
    )
    save_and_compare_images("multiple_plots")


def test_violin(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=40)

    with plt.rc_context():
        sc.pl.set_rcParams_defaults()
        sc.set_figure_params(dpi=50, color_map="viridis")

        pbmc = pbmc68k_reduced()
        sc.pl.violin(
            pbmc,
            ["n_genes", "percent_mito", "n_counts"],
            stripplot=True,
            multi_panel=True,
            jitter=True,
            show=False,
        )
        save_and_compare_images("violin_multi_panel")

        sc.pl.violin(
            pbmc,
            ["n_genes", "percent_mito", "n_counts"],
            ylabel=["foo", "bar", "baz"],
            groupby="bulk_labels",
            stripplot=True,
            multi_panel=True,
            jitter=True,
            show=False,
            rotation=90,
        )
        save_and_compare_images("violin_multi_panel_with_groupby")

        # test use of layer
        pbmc.layers["negative"] = pbmc.X * -1
        sc.pl.violin(
            pbmc,
            "CST3",
            groupby="bulk_labels",
            stripplot=True,
            multi_panel=True,
            jitter=True,
            show=False,
            layer="negative",
            use_raw=False,
            rotation=90,
        )
        save_and_compare_images("violin_multi_panel_with_layer")


# TODO: Generalize test to more plotting types
def test_violin_without_raw(tmp_path):
    # https://github.com/scverse/scanpy/issues/1546
    has_raw_pth = tmp_path / "has_raw.png"
    no_raw_pth = tmp_path / "no_raw.png"

    pbmc = pbmc68k_reduced()
    pbmc_no_raw = pbmc.raw.to_adata().copy()

    sc.pl.violin(pbmc, "CST3", groupby="bulk_labels", show=False, jitter=False)
    plt.savefig(has_raw_pth)
    plt.close()

    sc.pl.violin(pbmc_no_raw, "CST3", groupby="bulk_labels", show=False, jitter=False)
    plt.savefig(no_raw_pth)
    plt.close()

    assert compare_images(has_raw_pth, no_raw_pth, tol=5) is None


def test_dendrogram(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=10)

    pbmc = pbmc68k_reduced()
    sc.pl.dendrogram(pbmc, "bulk_labels")
    save_and_compare_images("dendrogram")


def test_correlation(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.pl.correlation_matrix(pbmc, "bulk_labels")
    save_and_compare_images("correlation")


_RANK_GENES_GROUPS_PARAMS = [
    (
        "sharey",
        partial(sc.pl.rank_genes_groups, n_genes=12, n_panels_per_row=3, show=False),
    ),
    (
        "basic",
        partial(
            sc.pl.rank_genes_groups,
            n_genes=12,
            n_panels_per_row=3,
            sharey=False,
            show=False,
        ),
    ),
    (
        "heatmap",
        partial(sc.pl.rank_genes_groups_heatmap, n_genes=4, cmap="YlGnBu", show=False),
    ),
    (
        "heatmap_swap_axes",
        partial(
            sc.pl.rank_genes_groups_heatmap,
            n_genes=20,
            swap_axes=True,
            use_raw=False,
            show_gene_labels=False,
            show=False,
            vmin=-3,
            vmax=3,
            cmap="bwr",
        ),
    ),
    (
        "heatmap_swap_axes_vcenter",
        partial(
            sc.pl.rank_genes_groups_heatmap,
            n_genes=20,
            swap_axes=True,
            use_raw=False,
            show_gene_labels=False,
            show=False,
            vmin=-3,
            vcenter=1,
            vmax=3,
            cmap="RdBu_r",
        ),
    ),
    (
        "stacked_violin",
        partial(
            sc.pl.rank_genes_groups_stacked_violin,
            n_genes=3,
            show=False,
            groups=["3", "0", "5"],
        ),
    ),
    (
        "dotplot",
        partial(sc.pl.rank_genes_groups_dotplot, n_genes=4, show=False),
    ),
    (
        "dotplot_gene_names",
        partial(
            sc.pl.rank_genes_groups_dotplot,
            var_names={
                "T-cell": ["CD3D", "CD3E", "IL32"],
                "B-cell": ["CD79A", "CD79B", "MS4A1"],
                "myeloid": ["CST3", "LYZ"],
            },
            values_to_plot="logfoldchanges",
            cmap="bwr",
            vmin=-3,
            vmax=3,
            show=False,
        ),
    ),
    (
        "dotplot_logfoldchange",
        partial(
            sc.pl.rank_genes_groups_dotplot,
            n_genes=4,
            values_to_plot="logfoldchanges",
            vmin=-5,
            vmax=5,
            min_logfoldchange=3,
            cmap="RdBu_r",
            swap_axes=True,
            title="log fold changes swap_axes",
            show=False,
        ),
    ),
    (
        "dotplot_logfoldchange_vcenter",
        partial(
            sc.pl.rank_genes_groups_dotplot,
            n_genes=4,
            values_to_plot="logfoldchanges",
            vmin=-5,
            vcenter=1,
            vmax=5,
            min_logfoldchange=3,
            cmap="RdBu_r",
            swap_axes=True,
            title="log fold changes swap_axes",
            show=False,
        ),
    ),
    (
        "matrixplot",
        partial(
            sc.pl.rank_genes_groups_matrixplot,
            n_genes=5,
            show=False,
            title="matrixplot",
            gene_symbols="symbol",
            use_raw=False,
        ),
    ),
    (
        "matrixplot_gene_names_symbol",
        partial(
            sc.pl.rank_genes_groups_matrixplot,
            var_names={
                "T-cell": ["CD3D__", "CD3E__", "IL32__"],
                "B-cell": ["CD79A__", "CD79B__", "MS4A1__"],
                "myeloid": ["CST3__", "LYZ__"],
            },
            values_to_plot="logfoldchanges",
            cmap="bwr",
            vmin=-3,
            vmax=3,
            gene_symbols="symbol",
            use_raw=False,
            show=False,
        ),
    ),
    (
        "matrixplot_n_genes_negative",
        partial(
            sc.pl.rank_genes_groups_matrixplot,
            n_genes=-5,
            show=False,
            title="matrixplot n_genes=-5",
        ),
    ),
    (
        "matrixplot_swap_axes",
        partial(
            sc.pl.rank_genes_groups_matrixplot,
            n_genes=5,
            show=False,
            swap_axes=True,
            values_to_plot="logfoldchanges",
            vmin=-6,
            vmax=6,
            cmap="bwr",
            title="log fold changes swap_axes",
        ),
    ),
    (
        "matrixplot_swap_axes_vcenter",
        partial(
            sc.pl.rank_genes_groups_matrixplot,
            n_genes=5,
            show=False,
            swap_axes=True,
            values_to_plot="logfoldchanges",
            vmin=-6,
            vcenter=1,
            vmax=6,
            cmap="bwr",
            title="log fold changes swap_axes",
        ),
    ),
    (
        "tracksplot",
        partial(
            sc.pl.rank_genes_groups_tracksplot,
            n_genes=3,
            show=False,
            groups=["3", "2", "1"],
        ),
    ),
    (
        "violin",
        partial(
            sc.pl.rank_genes_groups_violin,
            groups="0",
            n_genes=5,
            use_raw=True,
            jitter=False,
            strip=False,
            show=False,
        ),
    ),
    (
        "violin_not_raw",
        partial(
            sc.pl.rank_genes_groups_violin,
            groups="0",
            n_genes=5,
            use_raw=False,
            jitter=False,
            strip=False,
            show=False,
        ),
    ),
]


@pytest.mark.parametrize(
    ("name", "fn"),
    [pytest.param(name, fn, id=name) for name, fn in _RANK_GENES_GROUPS_PARAMS],
)
def test_rank_genes_groups(image_comparer, name, fn):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.tl.rank_genes_groups(pbmc, "louvain", n_genes=pbmc.raw.shape[1])

    # add gene symbol
    pbmc.var["symbol"] = pbmc.var.index + "__"

    with plt.rc_context({"axes.grid": True, "figure.figsize": (4, 4)}):
        fn(pbmc)
    key = "ranked_genes" if name == "basic" else f"ranked_genes_{name}"
    save_and_compare_images(key)
    plt.close()


def test_rank_genes_group_axes(image_comparer):
    fn = next(fn for name, fn in _RANK_GENES_GROUPS_PARAMS if name == "basic")

    save_and_compare_images = partial(image_comparer, ROOT, tol=23)

    pbmc = pbmc68k_reduced()
    sc.tl.rank_genes_groups(pbmc, "louvain", n_genes=pbmc.raw.shape[1])

    pbmc.var["symbol"] = pbmc.var.index + "__"

    fig, ax = plt.subplots(figsize=(12, 16))
    ax.set_axis_off()
    with plt.rc_context({"axes.grid": True}):
        axes: list[Axes] = fn(pbmc, ax=ax, show=False)

    assert len(axes) == 11
    fig.show()
    save_and_compare_images("ranked_genes")
    plt.close()


@pytest.fixture(scope="session")
def gene_symbols_adatas_session() -> tuple[AnnData, AnnData]:
    """Create two anndata objects which are equivalent except for var_names.

    Both have ensembl ids and hgnc symbols as columns in var. The first has ensembl
    ids as var_names, the second has symbols.
    """
    pbmc = pbmc3k_processed().raw.to_adata()
    pbmc_counts = pbmc3k()

    pbmc.layers["counts"] = pbmc_counts[pbmc.obs_names, pbmc.var_names].X.copy()
    pbmc.var["gene_symbol"] = pbmc.var_names
    pbmc.var["ensembl_id"] = pbmc_counts.var["gene_ids"].loc[pbmc.var_names]

    pbmc.var = pbmc.var.set_index("ensembl_id", drop=False)

    # Cutting down on size for plotting, tracksplot and stacked_violin are slow
    pbmc = pbmc[pbmc.obs["louvain"].isin(pbmc.obs["louvain"].cat.categories[:4])]
    pbmc = pbmc[::3].copy()

    # Creating variations
    a = pbmc.copy()
    b = pbmc.copy()
    a.var = a.var.set_index("ensembl_id")
    b.var = b.var.set_index("gene_symbol")

    # Computing DE
    sc.tl.rank_genes_groups(a, groupby="louvain")
    sc.tl.rank_genes_groups(b, groupby="louvain")

    return a, b


@pytest.fixture
def gene_symbols_adatas(gene_symbols_adatas_session) -> tuple[AnnData, AnnData]:
    a, b = gene_symbols_adatas_session
    return a.copy(), b.copy()


@pytest.mark.parametrize(
    "func",
    [
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        sc.pl.rank_genes_groups_tracksplot,
        # TODO: add other rank_genes_groups plots here once they work
    ],
)
def test_plot_rank_genes_groups_gene_symbols(
    gene_symbols_adatas, func, tmp_path, check_same_image
):
    a, b = gene_symbols_adatas

    pth_1_a = tmp_path / f"{func.__name__}_equivalent_gene_symbols_1_a.png"
    pth_1_b = tmp_path / f"{func.__name__}_equivalent_gene_symbols_1_b.png"

    func(a, gene_symbols="gene_symbol")
    plt.savefig(pth_1_a)
    plt.close()

    func(b)
    plt.savefig(pth_1_b)

    check_same_image(pth_1_a, pth_1_b, tol=1, root=tmp_path)

    pth_2_a = tmp_path / f"{func.__name__}_equivalent_gene_symbols_2_a.png"
    pth_2_b = tmp_path / f"{func.__name__}_equivalent_gene_symbols_2_b.png"

    func(a)
    plt.savefig(pth_2_a)
    plt.close()

    func(b, gene_symbols="ensembl_id")
    plt.savefig(pth_2_b)
    plt.close()

    check_same_image(pth_2_a, pth_2_b, tol=1, root=tmp_path)


@pytest.mark.parametrize(
    "func",
    [
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        sc.pl.rank_genes_groups_tracksplot,
        # TODO: add other rank_genes_groups plots here once they work
    ],
)
def test_rank_genes_groups_plots_n_genes_vs_var_names(tmp_path, func, check_same_image):
    """Checks that once can pass a negative value for n_genes and var_names as a dict."""
    N = 3
    pbmc = pbmc68k_reduced().raw.to_adata()
    groups = pbmc.obs["louvain"].cat.categories[:3]
    pbmc = pbmc[pbmc.obs["louvain"].isin(groups)][::3].copy()

    sc.tl.rank_genes_groups(pbmc, groupby="louvain")

    top_genes = {}
    bottom_genes = {}
    for g, subdf in sc.get.rank_genes_groups_df(pbmc, group=groups).groupby(
        "group", observed=True
    ):
        top_genes[g] = list(subdf["names"].head(N))
        bottom_genes[g] = list(subdf["names"].tail(N))

    positive_n_pth = tmp_path / f"{func.__name__}_positive_n.png"
    top_genes_pth = tmp_path / f"{func.__name__}_top_genes.png"
    negative_n_pth = tmp_path / f"{func.__name__}_negative_n.png"
    bottom_genes_pth = tmp_path / f"{func.__name__}_bottom_genes.png"

    def wrapped(pth, **kwargs):
        func(pbmc, groupby="louvain", dendrogram=False, **kwargs)
        plt.savefig(pth)
        plt.close()

    wrapped(positive_n_pth, n_genes=N)
    wrapped(top_genes_pth, var_names=top_genes)

    check_same_image(positive_n_pth, top_genes_pth, tol=1, root=tmp_path)

    wrapped(negative_n_pth, n_genes=-N)
    wrapped(bottom_genes_pth, var_names=bottom_genes)

    check_same_image(negative_n_pth, bottom_genes_pth, tol=1, root=tmp_path)

    # Shouldn't be able to pass these together
    with pytest.raises(
        ValueError, match="n_genes and var_names are mutually exclusive"
    ):
        wrapped(tmp_path / "not_written.png", n_genes=N, var_names=top_genes)


@pytest.mark.parametrize(
    ("id", "fn"),
    [
        ("heatmap", sc.pl.heatmap),
        ("dotplot", sc.pl.dotplot),
        ("matrixplot", sc.pl.matrixplot),
        ("stacked_violin", sc.pl.stacked_violin),
        ("tracksplot", sc.pl.tracksplot),
    ],
)
def test_genes_symbols(image_comparer, id, fn):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = krumsiek11()

    # add a 'symbols' column
    adata.var["symbols"] = adata.var.index.map(lambda x: f"symbol_{x}")
    symbols = [f"symbol_{x}" for x in adata.var_names]

    fn(adata, symbols, "cell_type", dendrogram=True, gene_symbols="symbols", show=False)
    save_and_compare_images(f"{id}_gene_symbols")


@pytest.fixture(scope="session")
def pbmc_scatterplots_session() -> AnnData:
    # Wrapped in another fixture to avoid mutation
    pbmc = pbmc68k_reduced()
    pbmc.obs["mask"] = pbmc.obs["louvain"].isin(["0", "1", "3"])
    pbmc.layers["sparse"] = pbmc.raw.X / 2
    pbmc.layers["test"] = pbmc.X.copy() + 100
    pbmc.var["numbers"] = [str(x) for x in range(pbmc.shape[1])]
    sc.pp.neighbors(pbmc)
    sc.tl.tsne(pbmc, random_state=0, n_pcs=30)
    sc.tl.diffmap(pbmc)
    return pbmc


@pytest.fixture
def pbmc_scatterplots(pbmc_scatterplots_session) -> AnnData:
    return pbmc_scatterplots_session.copy()


@pytest.mark.parametrize(
    ("id", "fn"),
    [
        ("pca", partial(sc.pl.pca, color="bulk_labels")),
        (
            "pca_with_fonts",
            partial(
                sc.pl.pca,
                color=["bulk_labels", "louvain"],
                legend_loc="on data",
                legend_fontoutline=2,
                legend_fontweight="normal",
                legend_fontsize=10,
            ),
        ),
        pytest.param(
            "3dprojection", partial(sc.pl.pca, color="bulk_labels", projection="3d")
        ),
        (
            "multipanel",
            partial(
                sc.pl.pca,
                color=["CD3D", "CD79A"],
                components=["1,2", "1,3"],
                vmax=5,
                use_raw=False,
                vmin=-5,
                cmap="seismic",
            ),
        ),
        (
            "multipanel_vcenter",
            partial(
                sc.pl.pca,
                color=["CD3D", "CD79A"],
                components=["1,2", "1,3"],
                vmax=5,
                use_raw=False,
                vmin=-5,
                vcenter=1,
                cmap="seismic",
            ),
        ),
        (
            "pca_one_marker",
            partial(sc.pl.pca, color="louvain", marker="^"),
        ),
        (
            "pca_one_marker_multiple_colors",
            partial(sc.pl.pca, color=["louvain", "bulk_labels"], marker="^"),
        ),
        (
            "pca_multiple_markers_multiple_colors",
            partial(sc.pl.pca, color=["louvain", "bulk_labels"], marker=["^", "x"]),
        ),
        (
            "pca_marker_with_dimensions",
            partial(
                sc.pl.pca, color="louvain", marker="^", dimensions=[(0, 1), (1, 2)]
            ),
        ),
        (
            "pca_markers_with_dimensions",
            partial(
                sc.pl.pca,
                color="louvain",
                marker=["^", "x"],
                dimensions=[(0, 1), (1, 2)],
            ),
        ),
        (
            "pca_markers_colors_with_dimensions",
            partial(
                sc.pl.pca,
                color=["louvain", "bulk_labels"],
                marker=["^", "x"],
                dimensions=[(0, 1), (1, 2)],
            ),
        ),
        (
            "pca_sparse_layer",
            partial(sc.pl.pca, color=["CD3D", "CD79A"], layer="sparse", cmap="viridis"),
        ),
        # pytest.param(
        #     "tsne",
        #     partial(sc.pl.tsne, color=["CD3D", "louvain"]),
        #     marks=pytest.mark.xfail(
        #         reason="slight differences even after setting random_state."
        #     ),
        # ),
        ("umap_nocolor", sc.pl.umap),
        (
            "umap",
            partial(
                sc.pl.umap,
                color=["louvain"],
                palette=["b", "grey80", "r", "yellow", "black", "gray", "lightblue"],
                frameon=False,
            ),
        ),
        (
            "umap_gene_expr",
            partial(
                sc.pl.umap,
                color=np.array(["LYZ", "CD79A"]),
                s=20,
                alpha=0.5,
                frameon=False,
                title=["gene1", "gene2"],
            ),
        ),
        (
            "umap_layer",
            partial(
                sc.pl.umap,
                color=np.array(["LYZ", "CD79A"]),
                s=20,
                alpha=0.5,
                frameon=False,
                title=["gene1", "gene2"],
                layer="test",
                vmin=100,
                vcenter=101,
            ),
        ),
        (
            "umap_with_edges",
            partial(sc.pl.umap, color="louvain", edges=True, edges_width=0.1, s=50),
        ),
        # ('diffmap', partial(sc.pl.diffmap, components='all', color=['CD3D'])),
        (
            "umap_symbols",
            partial(sc.pl.umap, color=["1", "2", "3"], gene_symbols="numbers"),
        ),
        (
            "pca_mask",
            partial(
                sc.pl.pca,
                color=["LYZ", "CD79A", "louvain"],
                mask_obs="mask",
            ),
        ),
    ],
)
def test_scatterplots(image_comparer, pbmc_scatterplots, id, fn):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    fn(pbmc_scatterplots, show=False)
    save_and_compare_images(id)


def test_scatter_embedding_groups_and_size(image_comparer):
    # test that the 'groups' parameter sorts
    # cells, such that the cells belonging to the groups are
    # plotted on top. This new ordering requires that the size
    # vector is also ordered (if given).
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.pl.embedding(
        pbmc,
        "umap",
        color=["bulk_labels"],
        groups=["CD14+ Monocyte", "Dendritic"],
        size=(np.arange(pbmc.shape[0]) / 40) ** 1.7,
    )
    save_and_compare_images("embedding_groups_size")


def test_scatter_embedding_add_outline_vmin_vmax_norm(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()

    sc.pl.embedding(
        pbmc,
        "X_umap",
        color=["percent_mito", "n_counts", "bulk_labels", "percent_mito"],
        s=200,
        frameon=False,
        add_outline=True,
        vmax=["p99.0", partial(np.percentile, q=90), None, 0.03],
        vmin=0.01,
        vcenter=[0.015, None, None, 0.025],
        outline_color=("#555555", "0.9"),
        outline_width=(0.5, 0.5),
        cmap="viridis_r",
        alpha=0.9,
        wspace=0.5,
    )
    save_and_compare_images("embedding_outline_vmin_vmax")


def test_scatter_embedding_add_outline_vmin_vmax_norm_ref(tmp_path, check_same_image):
    pbmc = pbmc68k_reduced()

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    norm = mpl.colors.LogNorm()
    with pytest.raises(
        ValueError, match="Passing both norm and vmin/vmax/vcenter is not allowed."
    ):
        sc.pl.embedding(
            pbmc,
            "X_umap",
            color=["percent_mito", "n_counts"],
            norm=norm,
            vmin=0,
            vmax=1,
            vcenter=0.5,
            cmap="RdBu_r",
        )

    try:
        from matplotlib.colors import TwoSlopeNorm as DivNorm
    except ImportError:
        # matplotlib<3.2
        from matplotlib.colors import DivergingNorm as DivNorm

    from matplotlib.colors import Normalize

    norm = Normalize(0, 10000)
    divnorm = DivNorm(200, 150, 6000)

    # allowed
    sc.pl.umap(
        pbmc,
        color=["n_counts", "bulk_labels", "percent_mito"],
        frameon=False,
        vmax=["p99.0", None, None],
        vcenter=[0.015, None, None],
        norm=[None, norm, norm],
        wspace=0.5,
    )

    sc.pl.umap(
        pbmc,
        color=["n_counts", "bulk_labels"],
        frameon=False,
        norm=norm,
        wspace=0.5,
    )
    plt.savefig(tmp_path / "umap_norm_fig0.png")
    plt.close()

    sc.pl.umap(
        pbmc,
        color=["n_counts", "bulk_labels"],
        frameon=False,
        norm=divnorm,
        wspace=0.5,
    )
    plt.savefig(tmp_path / "umap_norm_fig1.png")
    plt.close()

    sc.pl.umap(
        pbmc,
        color=["n_counts", "bulk_labels"],
        frameon=False,
        vcenter=200,
        vmin=150,
        vmax=6000,
        wspace=0.5,
    )
    plt.savefig(tmp_path / "umap_norm_fig2.png")
    plt.close()

    check_same_image(
        tmp_path / "umap_norm_fig1.png",
        tmp_path / "umap_norm_fig2.png",
        tol=1,
        root=tmp_path,
    )

    with pytest.raises(AssertionError):
        check_same_image(
            tmp_path / "umap_norm_fig1.png",
            tmp_path / "umap_norm_fig0.png",
            tol=1,
            root=tmp_path,
            save=False,
        )


def test_timeseries():
    adata = pbmc68k_reduced()
    sc.pp.neighbors(adata, n_neighbors=5, method="gauss", knn=False)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_branchings=1, n_dcs=10)
    sc.pl.dpt_timeseries(adata, as_heatmap=True)


def test_scatter_raw(tmp_path):
    pbmc = pbmc68k_reduced()[:100].copy()
    raw_pth = tmp_path / "raw.png"
    x_pth = tmp_path / "X.png"

    sc.pl.scatter(pbmc, color="HES4", basis="umap", use_raw=True)
    plt.savefig(raw_pth, dpi=60)
    plt.close()

    sc.pl.scatter(pbmc, color="HES4", basis="umap", use_raw=False)
    plt.savefig(x_pth, dpi=60)
    plt.close()

    comp = compare_images(str(raw_pth), str(x_pth), tol=5)
    assert "Error" in comp, "Plots should change depending on use_raw."


def test_binary_scatter(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    data = AnnData(
        np.asarray([[-1, 2, 0], [3, 4, 0], [1, 2, 0]]).T,
        obs=dict(binary=np.asarray([False, True, True])),
    )
    sc.pp.pca(data)
    sc.pl.pca(data, color="binary")
    if pkg_version("scikit-learn") >= Version("1.5.0rc1"):
        save_and_compare_images("binary_pca")
    else:
        save_and_compare_images("binary_pca_old")


def test_scatter_specify_layer_and_raw():
    pbmc = pbmc68k_reduced()
    pbmc.layers["layer"] = pbmc.raw.X.copy()
    with pytest.raises(ValueError, match=r"Cannot use both a layer and.*raw"):
        sc.pl.umap(pbmc, color="HES4", use_raw=True, layer="layer")


@pytest.mark.parametrize(
    "color", ["n_genes", "bulk_labels", ["n_genes", "bulk_labels"]]
)
def test_scatter_no_basis_per_obs(image_comparer, color):
    """Test scatterplot of per-obs points with no basis."""
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.pl.scatter(
        pbmc,
        x="HES4",
        y="percent_mito",
        color=color,
        use_raw=False,
        # palette only applies to categorical, i.e. color=='bulk_labels'
        palette="Set2",
    )
    color_str = color if isinstance(color, str) else "_".join(color)
    save_and_compare_images(f"scatter_HES_percent_mito_{color_str}")


def test_scatter_no_basis_per_var(image_comparer):
    """Test scatterplot of per-var points with no basis."""
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.pl.scatter(pbmc, x="AAAGCCTGGCTAAC-1", y="AAATTCGATGCACA-1", use_raw=False)
    save_and_compare_images("scatter_AAAGCCTGGCTAAC-1_vs_AAATTCGATGCACA-1")


@pytest.fixture
def pbmc_filtered() -> Callable[[], AnnData]:
    pbmc = pbmc68k_reduced()
    sc.pp.filter_genes(pbmc, min_cells=10)
    return pbmc.copy


@pytest.mark.parametrize("use_raw", [True, None])
def test_scatter_no_basis_raw(check_same_image, pbmc_filtered, tmp_path, use_raw):
    """Test scatterplots of raw layer with no basis."""
    adata = pbmc_filtered()

    sc.pl.scatter(adata.raw.to_adata(), x="EGFL7", y="F12", color="FAM185A")
    plt.savefig(path1 := tmp_path / "scatter-raw-to-adata.png")

    sc.pl.scatter(adata, x="EGFL7", y="F12", color="FAM185A", use_raw=use_raw)
    plt.savefig(path2 := tmp_path / f"scatter-{use_raw=}.png")
    plt.close()

    check_same_image(path1, path2, tol=15, root=tmp_path)


@pytest.mark.parametrize(
    ("x", "y", "color", "use_raw"),
    [
        # test that plotting fails with a ValueError if trying to plot
        # var_names only found in raw and use_raw is False
        ("EGFL7", "F12", "FAM185A", False),
        # test that plotting fails if one axis is a per-var value and the
        # other is a per-obs value
        ("HES4", "n_cells", None, None),
        ("percent_mito", "AAAGCCTGGCTAAC-1", None, None),
    ],
)
def test_scatter_no_basis_value_error(pbmc_filtered, x, y, color, use_raw):
    """Test that `scatter()` raises `ValueError` where appropriate.

    If `sc.pl.scatter()` receives variable labels that either cannot be
    found or are incompatible with one another, the function should
    raise a `ValueError`. This test checks that this happens as
    expected.
    """
    with pytest.raises(
        ValueError, match=r"inputs must all come from either `\.obs` or `\.var`"
    ):
        sc.pl.scatter(pbmc_filtered(), x=x, y=y, color=color, use_raw=use_raw)


def test_rankings(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()
    sc.pp.pca(pbmc)
    sc.pl.pca_loadings(pbmc)
    save_and_compare_images("pca_loadings")

    sc.pl.pca_loadings(pbmc, components="1,2,3")
    save_and_compare_images("pca_loadings")

    sc.pl.pca_loadings(pbmc, components=[1, 2, 3])
    save_and_compare_images("pca_loadings")

    sc.pl.pca_loadings(pbmc, include_lowest=False)
    save_and_compare_images("pca_loadings_without_lowest")

    sc.pl.pca_loadings(pbmc, n_points=10)
    save_and_compare_images("pca_loadings_10_points")


# TODO: Make more generic
def test_scatter_rep(tmp_path):
    """Test to make sure I can predict when scatter reps should be the same."""
    rep_args = {
        "raw": {"use_raw": True},
        "layer": {"layer": "layer", "use_raw": False},
        "X": {"use_raw": False},
    }
    states = pd.DataFrame.from_records(
        zip(
            list(chain.from_iterable(repeat(x, 3) for x in ["X", "raw", "layer"])),
            list(chain.from_iterable(repeat("abc", 3))),
            [1, 2, 3, 3, 1, 2, 2, 3, 1],
            strict=True,
        ),
        columns=["rep", "gene", "result"],
    )
    states["outpth"] = [
        tmp_path / f"{state.gene}_{state.rep}_{state.result}.png"
        for state in states.itertuples()
    ]
    pattern = np.array(list(chain.from_iterable(repeat(i, 5) for i in range(3))))
    coords = np.c_[np.arange(15) % 5, pattern]

    adata = AnnData(
        X=np.zeros((15, 3)),
        layers={"layer": np.zeros((15, 3))},
        obsm={"X_pca": coords},
        var=pd.DataFrame(index=list("abc")),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(15)]),
    )
    adata.raw = adata.copy()
    adata.X[np.arange(15), pattern] = 1
    adata.raw.X[np.arange(15), (pattern + 1) % 3] = 1
    adata.layers["layer"][np.arange(15), (pattern + 2) % 3] = 1

    for state in states.itertuples():
        sc.pl.pca(adata, color=state.gene, **rep_args[state.rep], show=False)
        plt.savefig(state.outpth, dpi=60)
        plt.close()

    for s1, s2 in combinations(states.itertuples(), 2):
        comp = compare_images(str(s1.outpth), str(s2.outpth), tol=5)
        if s1.result == s2.result:
            assert comp is None, comp
        else:
            assert "Error" in comp, f"{s1.outpth}, {s2.outpth} aren't supposed to match"


def test_no_copy():
    # https://github.com/scverse/scanpy/issues/1000
    # Tests that plotting functions don't make a copy from a view unless they
    # actually have to
    actual = pbmc68k_reduced()
    sc.pl.umap(actual, color=["bulk_labels", "louvain"], show=False)  # Set colors

    view = actual[np.random.choice(actual.obs_names, size=actual.shape[0] // 5), :]

    sc.pl.umap(view, color=["bulk_labels", "louvain"], show=False)
    assert view.is_view

    rank_genes_groups_plotting_funcs = [
        sc.pl.rank_genes_groups,
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        # TODO: raises ValueError about empty distance matrix – investigate
        # sc.pl.rank_genes_groups_tracksplot,
        sc.pl.rank_genes_groups_violin,
    ]

    # the pbmc68k was generated using rank_genes_groups with method='logreg'
    # which does not generate 'logfoldchanges', although this field is
    # required by `sc.get.rank_genes_groups_df`.
    # After updating rank_genes_groups plots to use the latter function
    # an error appears. Re-running rank_genes_groups with default method
    # solves the problem.
    sc.tl.rank_genes_groups(actual, "bulk_labels")

    # Only plotting one group at a time to avoid generating dendrogram
    # TODO: Generating a dendrogram modifies the object, this should be
    # optional and also maybe not modify the object.
    for plotfunc in rank_genes_groups_plotting_funcs:
        view = actual[actual.obs["bulk_labels"] == "Dendritic"]
        plotfunc(view, ["Dendritic"], show=False)
        assert view.is_view


def test_groupby_index(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    pbmc = pbmc68k_reduced()

    genes = [
        "CD79A",
        "MS4A1",
        "CD8A",
        "CD8B",
        "LYZ",
        "LGALS3",
        "S100A8",
        "GNLY",
        "NKG7",
        "KLRB1",
        "FCGR3A",
        "FCER1A",
        "CST3",
    ]
    pbmc_subset = pbmc[:10].copy()
    sc.pl.dotplot(pbmc_subset, genes, groupby="index")
    save_and_compare_images("dotplot_groupby_index")


# test category order when groupby is a list (#1735)
def test_groupby_list(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=30)

    adata = krumsiek11()

    np.random.seed(1)

    cat_val = adata.obs.cell_type.tolist()
    np.random.shuffle(cat_val)
    cats = adata.obs.cell_type.cat.categories.tolist()
    np.random.shuffle(cats)
    adata.obs["rand_cat"] = pd.Categorical(cat_val, categories=cats)

    with mpl.rc_context({"figure.subplot.bottom": 0.5}):
        sc.pl.dotplot(
            adata, ["Gata1", "Gata2"], groupby=["rand_cat", "cell_type"], swap_axes=True
        )
        save_and_compare_images("dotplot_groupby_list_catorder")


def test_color_cycler(caplog):
    # https://github.com/scverse/scanpy/issues/1885
    import logging

    pbmc = pbmc68k_reduced()
    colors = sns.color_palette("deep")
    cyl = sns.rcmod.cycler("color", sns.color_palette("deep"))

    with (
        caplog.at_level(logging.WARNING),
        plt.rc_context({"axes.prop_cycle": cyl, "patch.facecolor": colors[0]}),
    ):
        sc.pl.umap(pbmc, color="phase")
        plt.show()
        plt.close()

    assert caplog.text == ""


def test_repeated_colors_w_missing_value():
    # https://github.com/scverse/scanpy/issues/2133
    v = pd.Series(np.arange(10).astype(str))
    v[0] = np.nan
    v = v.astype("category")

    ad = sc.AnnData(obs=pd.DataFrame(v, columns=["value"]))
    ad.obsm["X_umap"] = np.random.normal(size=(ad.n_obs, 2))

    sc.pl.umap(ad, color="value")

    ad.uns["value_colors"][1] = ad.uns["value_colors"][0]

    sc.pl.umap(ad, color="value")


@pytest.mark.parametrize(
    "plot",
    [
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        sc.pl.rank_genes_groups_tracksplot,
        # TODO: add other rank_genes_groups plots here once they work
    ],
)
def test_filter_rank_genes_groups_plots(tmp_path, plot, check_same_image):
    N_GENES = 4

    adata = pbmc68k_reduced()

    sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon", pts=True)

    sc.tl.filter_rank_genes_groups(
        adata,
        key_added="rank_genes_groups_filtered",
        min_in_group_fraction=0.25,
        min_fold_change=1,
        max_out_group_fraction=0.5,
    )

    conditions = "logfoldchanges >= 1 & pct_nz_group >= .25 & pct_nz_reference < .5"
    df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups")
    df = df.query(conditions)[["group", "names"]]

    var_names = {
        k: v.head(N_GENES).tolist()
        for k, v in df.groupby("group", observed=True)["names"]
    }

    pth_a = tmp_path / f"{plot.__name__}_filter_a.png"
    pth_b = tmp_path / f"{plot.__name__}_filter_b.png"

    plot(adata, key="rank_genes_groups_filtered", n_genes=N_GENES)
    plt.savefig(pth_a)
    plt.close()

    plot(adata, key="rank_genes_groups", var_names=var_names)
    plt.savefig(pth_b)
    plt.close()

    check_same_image(pth_a, pth_b, tol=1, root=tmp_path)


@needs.skmisc
@pytest.mark.parametrize(
    ("id", "params"),
    [
        pytest.param("scrublet", {}, id="scrublet"),
        pytest.param("scrublet_no_threshold", {}, id="scrublet_no_threshold"),
        pytest.param(
            "scrublet_with_batches", dict(batch_key="batch"), id="scrublet_with_batches"
        ),
    ],
)
def test_scrublet_plots(monkeypatch, image_comparer, id, params):
    save_and_compare_images = partial(image_comparer, ROOT, tol=10)

    adata = pbmc3k()[:200].copy()
    adata.obs["batch"] = 100 * ["a"] + 100 * ["b"]

    with monkeypatch.context() as m:
        if id == "scrublet_no_threshold":
            m.setattr("skimage.filters.threshold_minimum", None)
        sc.pp.scrublet(adata, use_approx_neighbors=False, **params)
    if id == "scrublet_no_threshold":
        assert "threshold" not in adata.uns["scrublet"]

    sc.pl.scrublet_score_distribution(adata, return_fig=True, show=False)
    save_and_compare_images(id)


def test_umap_mask_equal(tmp_path, check_same_image):
    """Check that all desired cells are coloured and masked cells gray."""
    pbmc = pbmc3k_processed()
    mask_obs = pbmc.obs["louvain"].isin(["B cells", "NK cells"])

    ax = sc.pl.umap(pbmc, size=8.0, show=False)
    sc.pl.umap(pbmc[mask_obs], size=8.0, color="LDHB", ax=ax)
    plt.savefig(p1 := tmp_path / "umap_mask_fig1.png")
    plt.close()

    sc.pl.umap(pbmc, size=8.0, color="LDHB", mask_obs=mask_obs)
    plt.savefig(p2 := tmp_path / "umap_mask_fig2.png")
    plt.close()

    check_same_image(p1, p2, tol=1, root=tmp_path)


def test_umap_mask_mult_plots():
    """Check that multiple images are plotted when color is a list."""
    pbmc = pbmc3k_processed()
    color = ["LDHB", "LYZ", "CD79A"]
    mask_obs = pbmc.obs["louvain"].isin(["B cells", "NK cells"])
    axes = sc.pl.umap(pbmc, color=color, mask_obs=mask_obs, show=False)
    assert isinstance(axes, list)
    assert len(axes) == len(color)


def test_umap_categories_dont_change_when_rerun_with_fewer_categories():
    """Check that lowering the categories of interest does not cause a recalculation of colors."""
    pbmc = pbmc3k_processed()
    _ = sc.pl.umap(pbmc, color="louvain", show=False)
    assert len(pbmc.uns["louvain_colors"]) == len(pbmc.obs["louvain"].cat.categories)
    old_colors = pbmc.uns["louvain_colors"].copy()
    pbmc.obs.loc[pbmc.obs["louvain"] == "NK cells", "louvain"] = "B cells"
    pbmc.obs["louvain"] = pbmc.obs["louvain"].cat.remove_unused_categories()
    # see https://github.com/scverse/scanpy/issues/3716 for why this used to fail
    # Recalculation of the UMAP should not cause a re-calculation of colors
    # when there are fewer categories.
    _ = sc.pl.umap(pbmc, color="louvain", show=False)
    assert (old_colors == pbmc.uns["louvain_colors"]).all()


def test_umap_categories_change_when_rerun_with_more_categories():
    """Check that growing the categories of interest causes a recalculation of colors."""
    pbmc = pbmc3k_processed()
    _ = sc.pl.umap(pbmc, color="louvain", show=False)
    assert len(pbmc.uns["louvain_colors"]) == len(pbmc.obs["louvain"].cat.categories)
    pbmc.obs["louvain"] = pbmc.obs["louvain"].cat.add_categories("New Category")
    pbmc.obs.loc[pbmc.obs_names[:5], "louvain"] = "New Category"
    _ = sc.pl.umap(pbmc, color="louvain", show=False)
    assert len(pbmc.obs["louvain"].cat.categories) == len(pbmc.uns["louvain_colors"])


def test_umap_mask_no_modification():
    """Check that mask_obs argument doesn't affect the data being plotted."""
    pbmc = pbmc3k_processed()
    data_copy = pbmc.obs["louvain"].copy()
    sc.pl.umap(
        pbmc, mask_obs=(pbmc.obs["louvain"] == "B cells"), color="louvain", show=False
    )
    pd.testing.assert_series_equal(pbmc.obs["louvain"], data_copy)


def test_string_mask(tmp_path, check_same_image):
    """Check that the same mask given as string or bool array provides the same result."""
    pbmc = pbmc3k_processed()
    pbmc.obs["mask"] = mask_obs = pbmc.obs["louvain"].isin(["B cells", "NK cells"])

    sc.pl.umap(pbmc, mask_obs=mask_obs, color="LDHB")
    plt.savefig(p1 := tmp_path / "umap_mask_fig1.png")
    plt.close()

    sc.pl.umap(pbmc, color="LDHB", mask_obs="mask")
    plt.savefig(p2 := tmp_path / "umap_mask_fig2.png")
    plt.close()

    check_same_image(p1, p2, tol=1, root=tmp_path)


def test_violin_scale_warning(monkeypatch):
    adata = pbmc3k_processed()
    monkeypatch.setattr(sc.pl.StackedViolin, "DEFAULT_SCALE", "count", raising=False)
    with pytest.warns(FutureWarning, match="Don’t set DEFAULT_SCALE"):
        sc.pl.StackedViolin(adata, adata.var_names[:3], groupby="louvain")
