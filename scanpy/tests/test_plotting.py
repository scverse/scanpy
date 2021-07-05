from functools import partial
from pathlib import Path
import sys
from itertools import repeat, chain, combinations

import pytest
from matplotlib.testing import setup
from packaging import version

from scanpy._compat import pkg_version

setup()

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.testing.compare import compare_images
from anndata import AnnData

import scanpy as sc

HERE: Path = Path(__file__).parent
ROOT = HERE / '_images'
FIGS = HERE / 'figures'

sc.pl.set_rcParams_defaults()
sc.set_figure_params(dpi=40, color_map='viridis')

#####
# Test images are saved under the folder ./figures
# if test images need to be updated, simply copy them from
# the ./figures folder to ./_images/


def test_heatmap(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    sc.pl.heatmap(
        adata, adata.var_names, 'cell_type', use_raw=False, show=False, dendrogram=True
    )
    save_and_compare_images('master_heatmap')

    # test swap axes
    sc.pl.heatmap(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        show=False,
        dendrogram=True,
        swap_axes=True,
        figsize=(10, 3),
        cmap='YlGnBu',
    )
    save_and_compare_images('master_heatmap_swap_axes')

    # test heatmap numeric column():

    # set as numeric column the vales for the first gene on the matrix
    adata.obs['numeric_value'] = adata.X[:, 0]
    sc.pl.heatmap(
        adata,
        adata.var_names,
        'numeric_value',
        use_raw=False,
        num_categories=4,
        figsize=(4.5, 5),
        show=False,
    )
    save_and_compare_images('master_heatmap2')

    # test var/obs standardization and layer
    adata.layers['test'] = -1 * adata.X.copy()
    sc.pl.heatmap(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='var',
        layer='test',
    )
    save_and_compare_images('master_heatmap_std_scale_var')

    # test standard_scale_obs
    sc.pl.heatmap(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='obs',
    )
    save_and_compare_images('master_heatmap_std_scale_obs')

    # test var_names as dict
    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.leiden(pbmc, key_added="clusters", resolution=0.5)
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
    save_and_compare_images('master_heatmap_var_as_dict')

    # test that plot elements are well aligned
    # small
    a = AnnData(
        np.array([[0, 0.3, 0.5], [1, 1.3, 1.5], [2, 2.3, 2.5]]),
        obs={"foo": 'a b c'.split()},
        var=pd.DataFrame({"genes": 'g1 g2 g3'.split()}).set_index('genes'),
    )
    a.obs['foo'] = a.obs['foo'].astype('category')
    sc.pl.heatmap(
        a, var_names=a.var_names, groupby='foo', swap_axes=True, figsize=(4, 4)
    )
    save_and_compare_images('master_heatmap_small_swap_alignment')

    sc.pl.heatmap(
        a, var_names=a.var_names, groupby='foo', swap_axes=False, figsize=(4, 4)
    )
    save_and_compare_images('master_heatmap_small_alignment')


@pytest.mark.skipif(
    pkg_version("matplotlib") < version.parse('3.1'),
    reason="https://github.com/mwaskom/seaborn/issues/1953",
)
@pytest.mark.parametrize(
    "obs_keys,name",
    [(None, "master_clustermap"), ("cell_type", "master_clustermap_withcolor")],
)
def test_clustermap(image_comparer, obs_keys, name):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    adata = sc.datasets.krumsiek11()
    sc.pl.clustermap(adata, obs_keys)
    save_and_compare_images(name)


@pytest.mark.parametrize(
    "id,fn",
    [
        (
            "dotplot",
            partial(
                sc.pl.dotplot, groupby='cell_type', title='dotplot', dendrogram=True
            ),
        ),
        (
            "dotplot2",
            partial(
                sc.pl.dotplot,
                groupby='numeric_column',
                use_raw=False,
                num_categories=7,
                title='non categorical obs',
                figsize=(7, 2.5),
            ),
        ),
        (
            "dotplot3",
            partial(
                sc.pl.dotplot,
                groupby='cell_type',
                dot_max=0.7,
                dot_min=0.1,
                cmap='hot_r',
                title='dot_max=0.7 dot_min=0.1, var_groups',
                var_group_positions=[(0, 1), (9, 10)],
                var_group_labels=['A', 'B'],
                dendrogram=True,
            ),
        ),
        (
            "dotplot_std_scale_group",
            partial(
                sc.pl.dotplot,
                groupby='cell_type',
                use_raw=False,
                dendrogram=True,
                layer='test',
                swap_axes=True,
                title='swap_axes, layer=-1*X, scale=group\nsmallest_dot=10',
                standard_scale='group',
                smallest_dot=10,
            ),
        ),
        (
            "dotplot_dict",
            partial(
                sc.pl.dotplot,
                groupby='cell_type',
                dot_max=0.7,
                dot_min=0.1,
                color_map='winter',
                title='var as dict',
                dendrogram=True,
            ),
        ),
        (
            "matrixplot",
            partial(
                sc.pl.matrixplot,
                groupby='cell_type',
                use_raw=False,
                title='matrixplot',
                dendrogram=True,
            ),
        ),
        (
            "matrixplot_std_scale_var_dict",
            partial(
                sc.pl.matrixplot,
                groupby='cell_type',
                dendrogram=True,
                standard_scale='var',
                layer='test',
                cmap='Blues_r',
                title='scale var, custom colorbar_title, layer="test"',
                colorbar_title="Scaled expression",
            ),
        ),
        (
            "matrixplot_std_scale_group",
            partial(
                sc.pl.matrixplot,
                groupby='cell_type',
                use_raw=False,
                standard_scale='group',
                title='scale_group, swap_axes',
                swap_axes=True,
            ),
        ),
        (
            "matrixplot2",
            partial(
                sc.pl.matrixplot,
                groupby='numeric_column',
                use_raw=False,
                num_categories=4,
                title='non-categorical obs, custom figsize',
                figsize=(8, 2.5),
                cmap='RdBu_r',
            ),
        ),
        (
            "stacked_violin",
            partial(
                sc.pl.stacked_violin,
                groupby='cell_type',
                use_raw=False,
                title='stacked_violin',
                dendrogram=True,
            ),
        ),
        (
            "stacked_violin_std_scale_var_dict",
            partial(
                sc.pl.stacked_violin,
                groupby='cell_type',
                dendrogram=True,
                standard_scale='var',
                layer='test',
                title='scale var, layer="test"',
            ),
        ),
        (
            "stacked_violin_std_scale_group",
            partial(
                sc.pl.stacked_violin,
                groupby='cell_type',
                use_raw=False,
                standard_scale='group',
                title='scale_group\nswap_axes',
                swap_axes=True,
                cmap='Blues',
            ),
        ),
        (
            "stacked_violin_no_cat_obs",
            partial(
                sc.pl.stacked_violin,
                groupby='numeric_column',
                use_raw=False,
                num_categories=4,
                title='non-categorical obs, custom figsize',
                figsize=(8, 2.5),
            ),
        ),
    ],
)
def test_dotplot_matrixplot_stacked_violin(image_comparer, id, fn):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    adata.obs['numeric_column'] = adata.X[:, 0]
    adata.layers['test'] = -1 * adata.X.copy()
    genes_dict = {
        'group a': ['Gata2', 'Gata1'],
        'group b': ['Fog1', 'EKLF', 'Fli1', 'SCL'],
        'group c': ['Cebpa', 'Pu.1', 'cJun', 'EgrNab', 'Gfi1'],
    }

    if id.endswith("dict"):
        fn(adata, genes_dict, show=False)
    else:
        fn(adata, adata.var_names, show=False)
    save_and_compare_images(f"master_{id}")


def test_dotplot_obj(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    # test dotplot dot_min, dot_max, color_map, and var_groups
    pbmc = sc.datasets.pbmc68k_reduced()
    genes = [
        'CD79A',
        'MS4A1',
        'CD8A',
        'CD8B',
        'LYZ',
        'LGALS3',
        'S100A8',
        'GNLY',
        'NKG7',
        'KLRB1',
        'FCGR3A',
        'FCER1A',
        'CST3',
    ]
    # test layer, var standardization, smallest_dot,
    # color title, size_title return_fig and dot_edge
    pbmc.layers['test'] = pbmc.X * -1
    plot = sc.pl.dotplot(
        pbmc,
        genes,
        'bulk_labels',
        layer='test',
        dendrogram=True,
        return_fig=True,
        standard_scale='var',
        smallest_dot=40,
        colorbar_title='scaled column max',
        size_title='Fraction of cells',
    )
    plot.style(dot_edge_color='black', dot_edge_lw=0.1, cmap='Reds').show()

    save_and_compare_images('master_dotplot_std_scale_var')


def test_matrixplot_obj(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    adata = sc.datasets.pbmc68k_reduced()
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
        'bulk_labels',
        use_raw=False,
        title='added totals',
        return_fig=True,
    )
    plot.add_totals(sort='descending').style(edge_color='white', edge_lw=0.5).show()
    save_and_compare_images('master_matrixplot_with_totals')

    axes = plot.get_axes()
    assert 'mainplot_ax' in axes, 'mainplot_ax not found in returned axes dict'


def test_stacked_violin_obj(image_comparer, plt):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=26)

    pbmc = sc.datasets.pbmc68k_reduced()
    markers = {
        'T-cell': ['CD3D', 'CD3E', 'IL32'],
        'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
        'myeloid': ['CST3', 'LYZ'],
    }
    plot = sc.pl.stacked_violin(
        pbmc,
        markers,
        'bulk_labels',
        use_raw=False,
        title="return_fig. add_totals",
        return_fig=True,
    )
    plot.add_totals().style(row_palette='tab20').show()
    save_and_compare_images('master_stacked_violin_return_fig')


def test_tracksplot(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    sc.pl.tracksplot(
        adata, adata.var_names, 'cell_type', dendrogram=True, use_raw=False
    )
    save_and_compare_images('master_tracksplot')


def test_multiple_plots(image_comparer):
    # only testing stacked_violin, matrixplot and dotplot
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.pbmc68k_reduced()
    markers = {
        'T-cell': ['CD3D', 'CD3E', 'IL32'],
        'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
        'myeloid': ['CST3', 'LYZ'],
    }
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(20, 5), gridspec_kw={'wspace': 0.7}
    )
    _ = sc.pl.stacked_violin(
        adata,
        markers,
        groupby='bulk_labels',
        ax=ax1,
        title='stacked_violin',
        dendrogram=True,
        show=False,
    )
    _ = sc.pl.dotplot(
        adata,
        markers,
        groupby='bulk_labels',
        ax=ax2,
        title='dotplot',
        dendrogram=True,
        show=False,
    )
    _ = sc.pl.matrixplot(
        adata,
        markers,
        groupby='bulk_labels',
        ax=ax3,
        title='matrixplot',
        dendrogram=True,
        show=False,
    )
    save_and_compare_images('master_multiple_plots')


def test_violin(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=40)

    sc.pl.set_rcParams_defaults()
    sc.set_figure_params(dpi=50, color_map='viridis')

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pl.violin(
        pbmc,
        ['n_genes', 'percent_mito', 'n_counts'],
        stripplot=True,
        multi_panel=True,
        jitter=True,
        show=False,
    )
    save_and_compare_images('master_violin_multi_panel')

    sc.pl.violin(
        pbmc,
        ['n_genes', 'percent_mito', 'n_counts'],
        ylabel=["foo", "bar", "baz"],
        groupby='bulk_labels',
        stripplot=True,
        multi_panel=True,
        jitter=True,
        show=False,
        rotation=90,
    )
    save_and_compare_images('master_violin_multi_panel_with_groupby')

    # test use of layer
    pbmc.layers['negative'] = pbmc.X * -1
    sc.pl.violin(
        pbmc,
        'CST3',
        groupby='bulk_labels',
        stripplot=True,
        multi_panel=True,
        jitter=True,
        show=False,
        layer='negative',
        use_raw=False,
        rotation=90,
    )
    save_and_compare_images('master_violin_multi_panel_with_layer')


# TODO: Generalize test to more plotting types
def test_violin_without_raw(tmpdir):
    # https://github.com/theislab/scanpy/issues/1546
    TESTDIR = Path(tmpdir)

    has_raw_pth = TESTDIR / "has_raw.png"
    no_raw_pth = TESTDIR / "no_raw.png"

    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc_no_raw = pbmc.raw.to_adata().copy()

    sc.pl.violin(pbmc, 'CST3', groupby="bulk_labels", show=False)
    plt.savefig(has_raw_pth)
    plt.close()

    sc.pl.violin(pbmc_no_raw, 'CST3', groupby="bulk_labels", show=False)
    plt.savefig(no_raw_pth)
    plt.close()

    assert compare_images(has_raw_pth, no_raw_pth, tol=5) is None


def test_dendrogram(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=10)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pl.dendrogram(pbmc, 'bulk_labels')
    save_and_compare_images('dendrogram')


def test_correlation(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pl.correlation_matrix(pbmc, 'bulk_labels')
    save_and_compare_images('correlation')


@pytest.mark.parametrize(
    "name,fn",
    [
        (
            "ranked_genes_sharey",
            partial(
                sc.pl.rank_genes_groups, n_genes=12, n_panels_per_row=3, show=False
            ),
        ),
        (
            "ranked_genes",
            partial(
                sc.pl.rank_genes_groups,
                n_genes=12,
                n_panels_per_row=3,
                sharey=False,
                show=False,
            ),
        ),
        (
            "ranked_genes_heatmap",
            partial(
                sc.pl.rank_genes_groups_heatmap, n_genes=4, cmap='YlGnBu', show=False
            ),
        ),
        (
            "ranked_genes_heatmap_swap_axes",
            partial(
                sc.pl.rank_genes_groups_heatmap,
                n_genes=20,
                swap_axes=True,
                use_raw=False,
                show_gene_labels=False,
                show=False,
                vmin=-3,
                vmax=3,
                cmap='bwr',
            ),
        ),
        (
            "ranked_genes_heatmap_swap_axes_vcenter",
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
                cmap='RdBu_r',
            ),
        ),
        (
            "ranked_genes_stacked_violin",
            partial(
                sc.pl.rank_genes_groups_stacked_violin,
                n_genes=3,
                show=False,
                groups=['3', '0', '5'],
            ),
        ),
        (
            "ranked_genes_dotplot",
            partial(sc.pl.rank_genes_groups_dotplot, n_genes=4, show=False),
        ),
        (
            "ranked_genes_dotplot_gene_names",
            partial(
                sc.pl.rank_genes_groups_dotplot,
                var_names={
                    'T-cell': ['CD3D', 'CD3E', 'IL32'],
                    'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
                    'myeloid': ['CST3', 'LYZ'],
                },
                values_to_plot='logfoldchanges',
                cmap='bwr',
                vmin=-3,
                vmax=3,
                show=False,
            ),
        ),
        (
            "ranked_genes_dotplot_logfoldchange",
            partial(
                sc.pl.rank_genes_groups_dotplot,
                n_genes=4,
                values_to_plot="logfoldchanges",
                vmin=-5,
                vmax=5,
                min_logfoldchange=3,
                cmap='RdBu_r',
                swap_axes=True,
                title='log fold changes swap_axes',
                show=False,
            ),
        ),
        (
            "ranked_genes_dotplot_logfoldchange_vcenter",
            partial(
                sc.pl.rank_genes_groups_dotplot,
                n_genes=4,
                values_to_plot="logfoldchanges",
                vmin=-5,
                vcenter=1,
                vmax=5,
                min_logfoldchange=3,
                cmap='RdBu_r',
                swap_axes=True,
                title='log fold changes swap_axes',
                show=False,
            ),
        ),
        (
            "ranked_genes_matrixplot",
            partial(
                sc.pl.rank_genes_groups_matrixplot,
                n_genes=5,
                show=False,
                title='matrixplot',
                gene_symbols='symbol',
                use_raw=False,
            ),
        ),
        (
            "ranked_genes_matrixplot_gene_names_symbol",
            partial(
                sc.pl.rank_genes_groups_matrixplot,
                var_names={
                    'T-cell': ['CD3D__', 'CD3E__', 'IL32__'],
                    'B-cell': ['CD79A__', 'CD79B__', 'MS4A1__'],
                    'myeloid': ['CST3__', 'LYZ__'],
                },
                values_to_plot='logfoldchanges',
                cmap='bwr',
                vmin=-3,
                vmax=3,
                gene_symbols='symbol',
                use_raw=False,
                show=False,
            ),
        ),
        (
            "ranked_genes_matrixplot_n_genes_negative",
            partial(
                sc.pl.rank_genes_groups_matrixplot,
                n_genes=-5,
                show=False,
                title='matrixplot n_genes=-5',
            ),
        ),
        (
            "ranked_genes_matrixplot_swap_axes",
            partial(
                sc.pl.rank_genes_groups_matrixplot,
                n_genes=5,
                show=False,
                swap_axes=True,
                values_to_plot='logfoldchanges',
                vmin=-6,
                vmax=6,
                cmap='bwr',
                title='log fold changes swap_axes',
            ),
        ),
        (
            "ranked_genes_matrixplot_swap_axes_vcenter",
            partial(
                sc.pl.rank_genes_groups_matrixplot,
                n_genes=5,
                show=False,
                swap_axes=True,
                values_to_plot='logfoldchanges',
                vmin=-6,
                vcenter=1,
                vmax=6,
                cmap='bwr',
                title='log fold changes swap_axes',
            ),
        ),
        (
            "ranked_genes_tracksplot",
            partial(
                sc.pl.rank_genes_groups_tracksplot,
                n_genes=3,
                show=False,
                groups=['3', '2', '1'],
            ),
        ),
        (
            "ranked_genes_violin",
            partial(
                sc.pl.rank_genes_groups_violin,
                groups='0',
                n_genes=5,
                use_raw=True,
                jitter=False,
                strip=False,
                show=False,
            ),
        ),
        (
            "ranked_genes_violin_not_raw",
            partial(
                sc.pl.rank_genes_groups_violin,
                groups='0',
                n_genes=5,
                use_raw=False,
                jitter=False,
                strip=False,
                show=False,
            ),
        ),
    ],
)
def test_rank_genes_groups(image_comparer, name, fn):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.rank_genes_groups(pbmc, 'louvain', n_genes=pbmc.raw.shape[1])
    from matplotlib import rcParams

    rcParams['axes.grid'] = True
    rcParams['figure.figsize'] = 4, 4

    # add gene symbol
    pbmc.var['symbol'] = pbmc.var.index + "__"
    fn(pbmc)
    save_and_compare_images(f"master_{name}")
    plt.close()


@pytest.fixture(scope="session")
def gene_symbols_adatas():
    """Create two anndata objects which are equivalent except for var_names

    Both have ensembl ids and hgnc symbols as columns in var. The first has ensembl
    ids as var_names, the second has symbols.
    """
    pbmc = sc.datasets.pbmc3k_processed().raw.to_adata()
    pbmc_counts = sc.datasets.pbmc3k()

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


@pytest.mark.parametrize(
    "func",
    (
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        sc.pl.rank_genes_groups_tracksplot,
        # TODO: add other rank_genes_groups plots here once they work
    ),
)
def test_plot_rank_genes_groups_gene_symbols(
    gene_symbols_adatas, func, check_same_image
):
    a, b = gene_symbols_adatas

    pth_1_a = FIGS / f"{func.__name__}_equivalent_gene_symbols_1_a.png"
    pth_1_b = FIGS / f"{func.__name__}_equivalent_gene_symbols_1_b.png"

    func(a, gene_symbols="gene_symbol")
    plt.savefig(pth_1_a)
    plt.close()

    func(b)
    plt.savefig(pth_1_b)
    pass

    check_same_image(pth_1_a, pth_1_b, tol=1)

    pth_2_a = FIGS / f"{func.__name__}_equivalent_gene_symbols_2_a.png"
    pth_2_b = FIGS / f"{func.__name__}_equivalent_gene_symbols_2_b.png"

    func(a)
    plt.savefig(pth_2_a)
    plt.close()

    func(b, gene_symbols="ensembl_id")
    plt.savefig(pth_2_b)
    plt.close()

    check_same_image(pth_2_a, pth_2_b, tol=1)


@pytest.mark.parametrize(
    "func",
    (
        sc.pl.rank_genes_groups_dotplot,
        sc.pl.rank_genes_groups_heatmap,
        sc.pl.rank_genes_groups_matrixplot,
        sc.pl.rank_genes_groups_stacked_violin,
        sc.pl.rank_genes_groups_tracksplot,
        # TODO: add other rank_genes_groups plots here once they work
    ),
)
def test_rank_genes_groups_plots_n_genes_vs_var_names(tmpdir, func, check_same_image):
    """\
    Checks that passing a negative value for n_genes works, and that passing
    var_names as a dict works.
    """
    N = 3
    pbmc = sc.datasets.pbmc68k_reduced().raw.to_adata()
    groups = pbmc.obs["louvain"].cat.categories[:3]
    pbmc = pbmc[pbmc.obs["louvain"].isin(groups)][::3].copy()

    sc.tl.rank_genes_groups(pbmc, groupby="louvain")

    top_genes = {}
    bottom_genes = {}
    for g, subdf in sc.get.rank_genes_groups_df(pbmc, group=groups).groupby("group"):
        top_genes[g] = list(subdf["names"].head(N))
        bottom_genes[g] = list(subdf["names"].tail(N))

    positive_n_pth = tmpdir / f"{func.__name__}_positive_n.png"
    top_genes_pth = tmpdir / f"{func.__name__}_top_genes.png"
    negative_n_pth = tmpdir / f"{func.__name__}_negative_n.png"
    bottom_genes_pth = tmpdir / f"{func.__name__}_bottom_genes.png"

    def wrapped(pth, **kwargs):
        func(pbmc, groupby="louvain", dendrogram=False, **kwargs)
        plt.savefig(pth)
        plt.close()

    wrapped(positive_n_pth, n_genes=N)
    wrapped(top_genes_pth, var_names=top_genes)

    check_same_image(positive_n_pth, top_genes_pth, tol=1)

    wrapped(negative_n_pth, n_genes=-N)
    wrapped(bottom_genes_pth, var_names=bottom_genes)

    check_same_image(negative_n_pth, bottom_genes_pth, tol=1)

    # Shouldn't be able to pass these together
    with pytest.raises(
        ValueError, match="n_genes and var_names are mutually exclusive"
    ):
        wrapped(tmpdir / "not_written.png", n_genes=N, var_names=top_genes)


@pytest.mark.parametrize(
    "id,fn",
    [
        ("heatmap", sc.pl.heatmap),
        ("dotplot", sc.pl.dotplot),
        ("matrixplot", sc.pl.matrixplot),
        ("stacked_violin", sc.pl.stacked_violin),
        ("tracksplot", sc.pl.tracksplot),
    ],
)
def test_genes_symbols(image_comparer, id, fn):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()

    # add a 'symbols' column
    adata.var['symbols'] = adata.var.index.map(lambda x: "symbol_{}".format(x))
    symbols = ["symbol_{}".format(x) for x in adata.var_names]

    fn(adata, symbols, 'cell_type', dendrogram=True, gene_symbols='symbols', show=False)
    save_and_compare_images(f"master_{id}_gene_symbols")


@pytest.fixture(scope="module")
def pbmc_scatterplots():
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.layers["sparse"] = pbmc.raw.X / 2
    pbmc.layers["test"] = pbmc.X.copy() + 100
    pbmc.var["numbers"] = [str(x) for x in range(pbmc.shape[1])]
    sc.pp.neighbors(pbmc)
    sc.tl.tsne(pbmc, random_state=0, n_pcs=30)
    sc.tl.diffmap(pbmc)
    return pbmc


@pytest.mark.parametrize(
    'id,fn',
    [
        ('pca', partial(sc.pl.pca, color='bulk_labels')),
        (
            'pca_with_fonts',
            partial(
                sc.pl.pca,
                color=['bulk_labels', 'louvain'],
                legend_loc='on data',
                legend_fontoutline=2,
                legend_fontweight='normal',
                legend_fontsize=10,
            ),
        ),
        pytest.param(
            '3dprojection', partial(sc.pl.pca, color='bulk_labels', projection='3d')
        ),
        (
            'multipanel',
            partial(
                sc.pl.pca,
                color=['CD3D', 'CD79A'],
                components=['1,2', '1,3'],
                vmax=5,
                use_raw=False,
                vmin=-5,
                cmap='seismic',
            ),
        ),
        (
            'multipanel_vcenter',
            partial(
                sc.pl.pca,
                color=['CD3D', 'CD79A'],
                components=['1,2', '1,3'],
                vmax=5,
                use_raw=False,
                vmin=-5,
                vcenter=1,
                cmap='seismic',
            ),
        ),
        (
            'pca_sparse_layer',
            partial(sc.pl.pca, color=['CD3D', 'CD79A'], layer='sparse', cmap='viridis'),
        ),
        pytest.param(
            'tsne',
            partial(sc.pl.tsne, color=['CD3D', 'louvain']),
            marks=pytest.mark.xfail(
                reason='slight differences even after setting random_state.'
            ),
        ),
        ('umap_nocolor', sc.pl.umap),
        (
            'umap',
            partial(
                sc.pl.umap,
                color=['louvain'],
                palette=['b', 'grey80', 'r', 'yellow', 'black', 'gray', 'lightblue'],
                frameon=False,
            ),
        ),
        (
            'umap_gene_expr',
            partial(
                sc.pl.umap,
                color=np.array(['LYZ', 'CD79A']),
                s=20,
                alpha=0.5,
                frameon=False,
                title=['gene1', 'gene2'],
            ),
        ),
        (
            'umap_layer',
            partial(
                sc.pl.umap,
                color=np.array(['LYZ', 'CD79A']),
                s=20,
                alpha=0.5,
                frameon=False,
                title=['gene1', 'gene2'],
                layer='test',
                vmin=100,
                vcenter=101,
            ),
        ),
        (
            'umap_with_edges',
            partial(sc.pl.umap, color='louvain', edges=True, edges_width=0.1, s=50),
        ),
        # ('diffmap', partial(sc.pl.diffmap, components='all', color=['CD3D'])),
        (
            'umap_symbols',
            partial(sc.pl.umap, color=['1', '2', '3'], gene_symbols='numbers'),
        ),
    ],
)
def test_scatterplots(image_comparer, pbmc_scatterplots, id, fn):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    # https://github.com/theislab/scanpy/issues/849
    if id == "3dprojection" and version.parse(mpl.__version__) < version.parse("3.3.3"):
        with pytest.raises(ValueError, match=r"known error with matplotlib 3d"):
            fn(pbmc_scatterplots, show=False)
    else:
        fn(pbmc_scatterplots, show=False)
        save_and_compare_images(f"master_{id}")


def test_scatter_embedding_groups_and_size(image_comparer):
    # test that the 'groups' parameter sorts
    # cells, such that the cells belonging to the groups are
    # plotted on top. This new ordering requires that the size
    # vector is also ordered (if given).
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pl.embedding(
        pbmc,
        'umap',
        color=['bulk_labels'],
        groups=['CD14+ Monocyte', 'Dendritic'],
        size=(np.arange(pbmc.shape[0]) / 40) ** 1.7,
    )
    save_and_compare_images('master_embedding_groups_size')


def test_scatter_embedding_add_outline_vmin_vmax_norm(image_comparer, check_same_image):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    pbmc = sc.datasets.pbmc68k_reduced()

    sc.pl.embedding(
        pbmc,
        'X_umap',
        color=['percent_mito', 'n_counts', 'bulk_labels', 'percent_mito'],
        s=200,
        frameon=False,
        add_outline=True,
        vmax=['p99.0', partial(np.percentile, q=90), None, 0.03],
        vmin=0.01,
        vcenter=[0.015, None, None, 0.025],
        outline_color=('#555555', '0.9'),
        outline_width=(0.5, 0.5),
        cmap='viridis_r',
        alpha=0.9,
        wspace=0.5,
    )
    save_and_compare_images('master_embedding_outline_vmin_vmax')

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    norm = mpl.colors.LogNorm()
    with pytest.raises(
        ValueError, match="Passing both norm and vmin/vmax/vcenter is not allowed."
    ):
        sc.pl.embedding(
            pbmc,
            'X_umap',
            color=['percent_mito', 'n_counts'],
            norm=norm,
            vmin=0,
            vmax=1,
            vcenter=0.5,
            cmap='RdBu_r',
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
        color=['n_counts', 'bulk_labels', 'percent_mito'],
        frameon=False,
        vmax=['p99.0', None, None],
        vcenter=[0.015, None, None],
        norm=[None, norm, norm],
        wspace=0.5,
    )

    sc.pl.umap(
        pbmc,
        color=['n_counts', 'bulk_labels'],
        frameon=False,
        norm=norm,
        wspace=0.5,
    )
    plt.savefig(FIGS / 'umap_norm_fig0.png')
    plt.close()

    sc.pl.umap(
        pbmc,
        color=['n_counts', 'bulk_labels'],
        frameon=False,
        norm=divnorm,
        wspace=0.5,
    )
    plt.savefig(FIGS / 'umap_norm_fig1.png')
    plt.close()

    sc.pl.umap(
        pbmc,
        color=['n_counts', 'bulk_labels'],
        frameon=False,
        vcenter=200,
        vmin=150,
        vmax=6000,
        wspace=0.5,
    )
    plt.savefig(FIGS / 'umap_norm_fig2.png')
    plt.close()

    check_same_image(FIGS / 'umap_norm_fig1.png', FIGS / 'umap_norm_fig2.png', tol=1)

    with pytest.raises(AssertionError):
        check_same_image(
            FIGS / 'umap_norm_fig1.png', FIGS / 'umap_norm_fig0.png', tol=1
        )


def test_timeseries():
    adata = sc.datasets.pbmc68k_reduced()
    sc.pp.neighbors(adata, n_neighbors=5, method='gauss', knn=False)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_branchings=1, n_dcs=10)
    sc.pl.dpt_timeseries(adata, as_heatmap=True)


def test_scatter_raw(tmp_path):
    pbmc = sc.datasets.pbmc68k_reduced()[:100].copy()
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


def test_scatter_specify_layer_and_raw():
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.layers["layer"] = pbmc.raw.X.copy()
    with pytest.raises(ValueError):
        sc.pl.umap(pbmc, color="HES4", use_raw=True, layer="layer")


def test_rankings(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.pp.pca(pbmc)
    sc.pl.pca_loadings(pbmc)
    save_and_compare_images('master_pca_loadings')

    sc.pl.pca_loadings(pbmc, components='1,2,3')
    save_and_compare_images('master_pca_loadings')

    sc.pl.pca_loadings(pbmc, components=[1, 2, 3])
    save_and_compare_images('master_pca_loadings')

    sc.pl.pca_loadings(pbmc, include_lowest=False)
    save_and_compare_images('master_pca_loadings_without_lowest')


# TODO: Make more generic
def test_scatter_rep(tmpdir):
    """
    Test to make sure I can predict when scatter reps should be the same
    """
    TESTDIR = Path(tmpdir)
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
        ),
        columns=["rep", "gene", "result"],
    )
    states["outpth"] = [
        TESTDIR / f"{state.gene}_{state.rep}_{state.result}.png"
        for state in states.itertuples()
    ]
    pattern = np.array(list(chain.from_iterable(repeat(i, 5) for i in range(3))))
    coords = np.c_[np.arange(15) % 5, pattern]

    adata = AnnData(
        X=np.zeros((15, 3)),
        layers={"layer": np.zeros((15, 3))},
        obsm={"X_pca": coords},
        var=pd.DataFrame(index=[x for x in list("abc")]),
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


def test_paga(image_comparer):
    # Sometimes things shift a pixel or so, resulting in diffs up to ~27
    # The 1px-edges aren’t that good actually as they’re ignored at this tol …
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=30)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.paga(pbmc, groups='bulk_labels')

    common = dict(threshold=0.5, max_edge_width=1.0, random_state=0, show=False)

    # delete bulk_labels_colors to test the creation of color list by paga
    del pbmc.uns['bulk_labels_colors']
    sc.pl.paga(pbmc, **common)
    save_and_compare_images('master_paga')

    sc.pl.paga(pbmc, color='CST3', **common)
    save_and_compare_images('master_paga_continuous')

    pbmc.obs['cool_feature'] = pbmc[:, 'CST3'].X.squeeze()
    sc.pl.paga(pbmc, color='cool_feature', **common)
    save_and_compare_images('master_paga_continuous_obs')

    sc.pl.paga(pbmc, color=['CST3', 'GATA2'], **common)
    save_and_compare_images('master_paga_continuous_multiple')

    sc.pl.paga_compare(pbmc, legend_fontoutline=2, **common)
    save_and_compare_images('master_paga_compare')

    sc.pl.paga_compare(pbmc, color='CST3', legend_fontsize=5, **common)
    save_and_compare_images('master_paga_compare_continuous')

    sc.pl.paga_compare(pbmc, basis='X_pca', legend_fontweight='normal', **common)
    save_and_compare_images('master_paga_compare_pca')

    colors = {
        c: {cm.Set1(_): 0.33 for _ in range(3)}
        for c in pbmc.obs["bulk_labels"].cat.categories
    }
    colors["Dendritic"] = {cm.Set2(_): 0.25 for _ in range(4)}

    sc.pl.paga(pbmc, color=colors, colorbar=False)
    save_and_compare_images('master_paga_pie')


def test_paga_path(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc68k_reduced()
    sc.tl.paga(pbmc, groups='bulk_labels')

    pbmc.uns['iroot'] = 0
    sc.tl.dpt(pbmc)
    sc.pl.paga_path(
        pbmc,
        nodes=['Dendritic'],
        keys=['HES4', 'SRM', 'CSTB'],
        show=False,
    )
    save_and_compare_images('master_paga_path')


def test_paga_compare(image_comparer):
    # Tests that https://github.com/theislab/scanpy/issues/1887 is fixed
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc3k_processed()
    sc.tl.paga(pbmc, groups="louvain")

    sc.pl.paga_compare(pbmc, basis="umap", show=False)

    save_and_compare_images('master_paga_compare_pbmc3k')


def test_no_copy():
    # https://github.com/theislab/scanpy/issues/1000
    # Tests that plotting functions don't make a copy from a view unless they
    # actually have to
    actual = sc.datasets.pbmc68k_reduced()
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
    sc.tl.rank_genes_groups(actual, 'bulk_labels')

    # Only plotting one group at a time to avoid generating dendrogram
    # TODO: Generating a dendrogram modifies the object, this should be
    # optional and also maybe not modify the object.
    for plotfunc in rank_genes_groups_plotting_funcs:
        view = actual[actual.obs["bulk_labels"] == "Dendritic"]
        plotfunc(view, ["Dendritic"], show=False)
        assert view.is_view


def test_groupby_index(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    pbmc = sc.datasets.pbmc68k_reduced()

    genes = [
        'CD79A',
        'MS4A1',
        'CD8A',
        'CD8B',
        'LYZ',
        'LGALS3',
        'S100A8',
        'GNLY',
        'NKG7',
        'KLRB1',
        'FCGR3A',
        'FCER1A',
        'CST3',
    ]
    pbmc_subset = pbmc[:10].copy()
    sc.pl.dotplot(pbmc_subset, genes, groupby='index')
    save_and_compare_images('master_dotplot_groupby_index')


# test category order when groupby is a list (#1735)
def test_groupby_list(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=30)
    adata = sc.datasets.krumsiek11()

    np.random.seed(1)

    cat_val = adata.obs.cell_type.tolist()
    np.random.shuffle(cat_val)
    cats = adata.obs.cell_type.cat.categories.tolist()
    np.random.shuffle(cats)
    adata.obs['rand_cat'] = pd.Categorical(cat_val, categories=cats)

    with mpl.rc_context({"figure.subplot.bottom": 0.5}):
        sc.pl.dotplot(
            adata, ['Gata1', 'Gata2'], groupby=['rand_cat', 'cell_type'], swap_axes=True
        )
        save_and_compare_images('master_dotplot_groupby_list_catorder')


def test_color_cycler(caplog):
    # https://github.com/theislab/scanpy/issues/1885
    import logging

    pbmc = sc.datasets.pbmc68k_reduced()
    colors = sns.color_palette("deep")
    cyl = sns.rcmod.cycler('color', sns.color_palette("deep"))

    with caplog.at_level(logging.WARNING):
        with plt.rc_context({'axes.prop_cycle': cyl, "patch.facecolor": colors[0]}):
            sc.pl.umap(pbmc, color="phase")
            plt.show()
            plt.close()

    assert caplog.text == ""
