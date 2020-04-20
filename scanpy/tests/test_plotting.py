from functools import partial
from pathlib import Path
from itertools import repeat, chain, combinations

import pytest
from matplotlib.testing import setup
from packaging import version

from scanpy._utils import pkg_version

setup()

import matplotlib.pyplot as plt
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
        adata, adata.var_names, 'cell_type', use_raw=False, show=False, dendrogram=True,
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
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.heatmap(
        adata,
        adata.var_names,
        'Gata2',
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


def test_dotplot(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    sc.pl.dotplot(
        adata, adata.var_names, 'cell_type', use_raw=False, dendrogram=True, show=False,
    )
    save_and_compare_images('master_dotplot')

    # test dotplot numeric column():
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.dotplot(
        adata,
        adata.var_names,
        'Gata2',
        use_raw=False,
        num_categories=7,
        figsize=(7, 2.5),
        show=False,
    )
    save_and_compare_images('master_dotplot2')

    # test dotplot dot_min, dot_max, color_map, and var_groups
    pbmc = sc.datasets.pbmc68k_reduced()
    marker_genes = [
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
    sc.pl.dotplot(
        pbmc,
        marker_genes,
        groupby='louvain',
        dot_max=0.7,
        dot_min=0.1,
        color_map='hot_r',
        var_group_positions=[(0, 1), (11, 12)],
        var_group_labels=['B cells', 'dendritic'],
        figsize=(7, 2.5),
        dendrogram=True,
        show=False,
    )
    save_and_compare_images('master_dotplot3')

    # test dict as markers input
    markers_dict = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
    sc.pl.dotplot(
        pbmc,
        markers_dict,
        groupby='bulk_labels',
        dot_max=0.7,
        dot_min=0.1,
        color_map='winter',
        figsize=(7, 2.5),
        dendrogram=True,
        show=False,
    )
    save_and_compare_images('master_dotplot_dict')

    # test var/group standardization smallest_dot
    sc.pl.dotplot(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='var',
        smallest_dot=40,
    )
    save_and_compare_images('master_dotplot_std_scale_var')

    sc.pl.dotplot(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='group',
        smallest_dot=10,
    )
    save_and_compare_images('master_dotplot_std_scale_group')


def test_matrixplot(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    sc.pl.matrixplot(
        adata, adata.var_names, 'cell_type', use_raw=False, dendrogram=True, show=False,
    )
    save_and_compare_images('master_matrixplot')

    # test swap_axes
    sc.pl.matrixplot(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        swap_axes=True,
    )
    save_and_compare_images('master_matrixplot_swap_axes')

    # test var/group standardization and layer
    adata.layers['test'] = -1 * adata.X.copy()
    sc.pl.matrixplot(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='var',
        layer='test',
        cmap='Blues_r',
    )
    save_and_compare_images('master_matrixplot_std_scale_var')

    sc.pl.matrixplot(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        standard_scale='group',
        swap_axes=True,
    )
    save_and_compare_images('master_matrixplot_std_scale_group')

    # test matrixplot numeric column and alternative cmap
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.matrixplot(
        adata,
        adata.var_names,
        'Gata2',
        use_raw=False,
        num_categories=4,
        figsize=(8, 2.5),
        cmap='viridis',
        show=False,
    )
    save_and_compare_images('master_matrixplot2')


def test_stacked_violin(image_comparer, plt):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=26)

    adata = sc.datasets.krumsiek11()
    sc.pl.stacked_violin(
        adata, adata.var_names, 'cell_type', use_raw=False, color='blue', show=False,
    )

    plt.title("image may have cut labels.\nThis is ok for test")
    save_and_compare_images('master_stacked_violin')

    # test swapped axes
    sc.pl.stacked_violin(
        adata,
        adata.var_names,
        'cell_type',
        use_raw=False,
        swap_axes=True,
        figsize=(3, 5),
        show=False,
    )
    save_and_compare_images('master_stacked_violin_swapped_axes')


def test_tracksplot(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()
    sc.pl.tracksplot(
        adata, adata.var_names, 'cell_type', dendrogram=True, use_raw=False
    )
    save_and_compare_images('master_tracksplot')


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


def test_rank_genes_groups(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    pbmc = sc.datasets.pbmc68k_reduced()

    from matplotlib import rcParams

    rcParams['axes.grid'] = True
    rcParams['figure.figsize'] = 4, 4

    sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3, show=False)
    save_and_compare_images('master_ranked_genes_sharey')

    # test ranked genes panels sharey = False
    sc.pl.rank_genes_groups(
        pbmc, n_genes=12, n_panels_per_row=3, sharey=False, show=False
    )
    save_and_compare_images('master_ranked_genes')

    # test ranked genes using heatmap
    sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=5, cmap='YlGnBu', show=False)
    save_and_compare_images('master_ranked_genes_heatmap')

    # test ranked genes using heatmap (swap_axes=True show_gene_labels=False)
    sc.pl.rank_genes_groups_heatmap(
        pbmc,
        n_genes=20,
        swap_axes=True,
        use_raw=False,
        show_gene_labels=False,
        show=False,
        vmin=-3,
        vmax=3,
        cmap='bwr',
    )
    save_and_compare_images('master_ranked_genes_heatmap_swap_axes')

    # test ranked genes using stacked violin plots
    sc.pl.rank_genes_groups_stacked_violin(pbmc, n_genes=3, show=False)
    save_and_compare_images('master_ranked_genes_stacked_violin', tolerance=20)

    # test ranked genes using dotplot
    sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4, show=False)
    save_and_compare_images('master_ranked_genes_dotplot')

    # test ranked genes using matrixplot
    sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=5, show=False)
    save_and_compare_images('master_ranked_genes_matrixplot')

    # test ranked genes using matrixplot
    sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=5, show=False)
    save_and_compare_images('master_ranked_genes_matrixplot')

    # test ranked genes using matrixplot (swap_axes=True)
    sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=5, swap_axes=True, show=False)
    save_and_compare_images('master_ranked_genes_matrixplot_swap_axes')

    # test ranked genes using tracks_plot
    sc.pl.rank_genes_groups_tracksplot(pbmc, n_genes=5, show=False)
    save_and_compare_images('master_ranked_genes_tracksplot')

    # # test ranked genes using violin plots
    # sc.pl.rank_genes_groups_violin(pbmc, groups=pbmc.obs.bulk_labels.cat.categories[0], n_genes=5,
    #                                jitter=False, strip=False, show=False)
    # save_and_compare_images('master_ranked_genes_stacked_violin', tolerance=tolerance)


def test_rank_genes_symbols(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    adata = sc.datasets.krumsiek11()

    # add a 'symbols' column
    adata.var['symbols'] = adata.var.index.map(lambda x: "symbol_{}".format(x))
    symbols = ["symbol_{}".format(x) for x in adata.var_names]
    sc.pl.heatmap(
        adata,
        symbols,
        'cell_type',
        use_raw=False,
        show=False,
        dendrogram=True,
        gene_symbols='symbols',
    )
    save_and_compare_images('master_heatmap_gene_symbols')

    sc.pl.dotplot(
        adata,
        symbols,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        gene_symbols='symbols',
    )

    save_and_compare_images('master_dotplot_gene_symbols')

    sc.pl.matrixplot(
        adata,
        symbols,
        'cell_type',
        use_raw=False,
        dendrogram=True,
        show=False,
        gene_symbols='symbols',
    )

    save_and_compare_images('master_matrixplot_gene_symbols')

    sc.pl.stacked_violin(
        adata,
        symbols,
        'cell_type',
        use_raw=False,
        color='blue',
        show=False,
        gene_symbols='symbols',
    )
    save_and_compare_images('master_stacked_violin_gene_symbols', tolerance=21)

    sc.pl.tracksplot(
        adata,
        symbols,
        'cell_type',
        dendrogram=True,
        use_raw=False,
        gene_symbols='symbols',
    )
    save_and_compare_images('master_tracksplot_gene_symbols')


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
    "id,fn",
    [
        ("pca", partial(sc.pl.pca, color='bulk_labels')),
        (
            "pca_with_fonts",
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
            "3dprojection", partial(sc.pl.pca, color='bulk_labels', projection='3d'),
        ),
        (
            "multipanel",
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
            "pca_sparse_layer",
            partial(
                sc.pl.pca, color=['CD3D', 'CD79A'], layer="sparse", cmap='viridis',
            ),
        ),
        pytest.param(
            "tsne",
            partial(sc.pl.tsne, color=['CD3D', 'louvain']),
            marks=pytest.mark.xfail(
                reason="slight differences even after setting random_state."
            ),
        ),
        ("umap_nocolor", sc.pl.umap),
        (
            "umap",
            partial(
                sc.pl.umap,
                color=['louvain'],
                palette=['b', 'grey80', 'r', 'yellow', 'black', 'gray', 'lightblue'],
                frameon=False,
            ),
        ),
        (
            "umap_gene_expr",
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
            "umap_layer",
            partial(
                sc.pl.umap,
                color=np.array(['LYZ', 'CD79A']),
                s=20,
                alpha=0.5,
                frameon=False,
                title=['gene1', 'gene2'],
                layer='test',
                vmin=100,
            ),
        ),
        (
            "umap_with_edges",
            partial(sc.pl.umap, color='louvain', edges=True, edges_width=0.1, s=50),
        ),
        # ("diffmap", partial(sc.pl.diffmap, components='all', color=['CD3D'])),
        (
            "umap_symbols",
            partial(sc.pl.umap, color=['1', '2', '3'], gene_symbols="numbers"),
        ),
    ],
)
def test_scatterplots(image_comparer, pbmc_scatterplots, id, fn):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)

    if id == "3dprojection":
        # check if this still happens so we can remove our checks once mpl is fixed
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


def test_scatter_embedding_add_outline_vmin_vmax(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    pbmc = sc.datasets.pbmc68k_reduced()

    sc.pl.embedding(
        pbmc,
        'X_umap',
        color=['percent_mito', 'n_counts', 'bulk_labels'],
        s=200,
        frameon=False,
        add_outline=True,
        vmax=['p99.0', partial(np.percentile, q=90)],
        vmin=0.01,
        outline_color=('#555555', '0.9'),
        outline_width=(0.5, 0.5),
        cmap='viridis_r',
        alpha=0.9,
    )
    save_and_compare_images('master_embedding_outline_vmin_vmax')


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

    # Only plotting one group at a time to avoid generating dendrogram
    # TODO: Generating a dendrogram modifies the object, this should be
    # optional and also maybe not modify the object.
    for plotfunc in rank_genes_groups_plotting_funcs:
        view = actual[actual.obs["bulk_labels"] == "Dendritic"]
        plotfunc(view, ["Dendritic"], show=False)
        assert view.is_view


def test_visium_circles(image_comparer):
    save_and_compare_images = image_comparer(ROOT, FIGS, tol=15)
    adata = sc.read_visium(HERE / '_data' / 'visium_data' / '1.0.0')
    adata.obs = adata.obs.astype({'array_row': 'str'})

    sc.pl.spatial(
        adata,
        color="array_row",
        groups=["24", "33"],
        crop_coord=(100, 400, 400, 100),
        alpha=0.5,
        size=1.3,
    )

    save_and_compare_images('master_spatial_visium')
