# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl
from matplotlib.testing.compare import compare_images

import os.path
from tempfile import NamedTemporaryFile

import scanpy.api as sc

ROOT = os.path.dirname(os.path.abspath(__file__)) + '/_images/'

tolerance = 15

sc.pl.set_rcParams_defaults()
sc.set_figure_params(dpi=80, color_map='viridis')

#####
# To update the test images:
# cd scanpy/test/_images/
# python remake_test_images.py


def test_heatmap():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_heatmap_', delete=False)
    sc.pl.heatmap(adata, adata.var_names, 'cell_type', use_raw=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_heatmap.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)

    # test heatmap numeric column():

    # set as numeric column the vales for the first gene on the matrix
    adata.obs['Gata2'] = adata.X[:, 0]
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_heatmap2_', delete=False)
    sc.pl.heatmap(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=4, figsize=(4.5, 5), show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_heatmap2.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_dotplot():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_dotplot_', delete=False)
    sc.pl.dotplot(adata, adata.var_names, 'cell_type', use_raw=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    tolerance = 15
    res = compare_images(ROOT + '/master_dotplot.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    # test dotplot numeric column():
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=7, figsize=(7, 2.5), show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_dotplot2.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)


def test_matrixplot():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_matrixplot_', delete=False)
    sc.pl.matrixplot(adata, adata.var_names, 'cell_type', use_raw=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    tolerance = 15
    res = compare_images(ROOT + '/master_matrixplot.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    # test matrixplot numeric column and alternative cmap
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.matrixplot(adata, adata.var_names, 'Gata2', use_raw=False,
                     num_categories=4, figsize=(8, 2.5), cmap='viridis', show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_matrixplot2.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)


def test_stacked_violin():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_stacked_violin_', delete=False)
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, color='blue', show=False)

    pl.title("image may have cut labels.\nThis is ok for test")
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_stacked_violin.png', outfile.name, tolerance * 2)

    assert res is None, res

    os.remove(outfile.name)

    # test swapped axes
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False,
                         swap_axes=True, figsize=(3, 5), show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()
    res = compare_images(ROOT + '/master_stacked_violin_swapped_axes.png', outfile.name, tolerance * 2)

    assert res is None, res

    os.remove(outfile.name)


def test_violin():
    sc.pl.set_rcParams_defaults()
    sc.set_figure_params(dpi=80, color_map='viridis')

    pbmc = sc.datasets.pbmc68k_reduced()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_violin_', delete=False)
    sc.pl.violin(pbmc, ['n_genes', 'percent_mito', 'n_counts'],
                 stripplot=True, multi_panel=True, jitter=True, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_violin_multi_panel.png', outfile.name, 40)

    assert res is None, res

    os.remove(outfile.name)


def test_rank_genes_groups():
    pbmc = sc.datasets.pbmc68k_reduced()
    tolerance = 15

    # test ranked genes panels
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_rank_genes_groups_', delete=False)

    sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes_sharey.png', outfile.name, tolerance)
    assert res is None, res

    # test ranked genes panels sharey = False
    sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3, sharey=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes.png', outfile.name, tolerance)
    assert res is None, res

    # test ranked genes using heatmap
    sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=5, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes_heatmap.png', outfile.name, tolerance)
    assert res is None, res

    # test ranked genes using stacked violin plots
    sc.pl.rank_genes_groups_stacked_violin(pbmc, n_genes=3, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes_stacked_violin.png', outfile.name, tolerance)
    assert res is None, res

    # test ranked genes using dotplot
    sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes_dotplot.png', outfile.name, tolerance)
    assert res is None, res

    # test ranked genes using matrixplot
    sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=5, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_ranked_genes_matrixplot.png', outfile.name, tolerance)
    assert res is None, res

    # # test ranked genes using violin plots
    # sc.pl.rank_genes_groups_violin(pbmc, groups=pbmc.obs.bulk_labels.cat.categories[0], n_genes=5,
    #                                jitter=False, strip=False, show=False)
    # pl.savefig(outfile.name, dpi=80)
    # pl.close()
    #
    # res = compare_images(ROOT + '/master_ranked_genes_violin.png', outfile.name, tolerance)
    # assert res is None, res

    os.remove(outfile.name)


def test_scatterplots():

    pbmc = sc.datasets.pbmc68k_reduced()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_scatter_', delete=False)

    # test pca
    sc.pl.pca(pbmc, color='bulk_labels', show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_pca.png', outfile.name, tolerance)
    assert res is None, res

    # test projection='3d'
    sc.pl.pca(pbmc, color='bulk_labels', projection='3d', show=False)
    pl.savefig(outfile.name, dpi=80, show=False)
    pl.close()

    res = compare_images(ROOT + '/master_3dprojection.png', outfile.name, tolerance)
    assert res is None, res

    sc.pl.pca(pbmc, color=['CD3D', 'CD79A'], components=['1,2', '1,3'],
              vmax=5, use_raw=False, vmin=-5, cmap='seismic', show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_multipanel.png', outfile.name, tolerance)
    assert res is None, res

    # test tsne
    # I am removing this test because  slight differences are present even
    # after setting a random_state.
    # sc.tl.tsne(pbmc, random_state=0, n_pcs=30)
    # sc.pl.tsne(pbmc, color=['CD3D', 'louvain'], show=False)
    # pl.savefig(outfile.name, dpi=80)
    # pl.close()

    # res = compare_images(ROOT + '/master_tsne.png', outfile.name, tolerance)
    # assert res is None, res

    # test umap with louvain clusters and palette
    sc.pl.umap(pbmc, color=['louvain'],
               palette=['b', 'g', 'r', 'yellow', 'black', 'gray', 'lightblue'],
               frameon=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_umap.png', outfile.name, tolerance)
    assert res is None, res

    # test umap with gene expression
    sc.pl.umap(pbmc, color=['LYZ', 'CD79A'], s=20, alpha=0.5, frameon=False, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_umap_gene_expr.png', outfile.name, tolerance)
    assert res is None, res

    # test edges = True
    sc.pp.neighbors(pbmc)
    sc.pl.umap(pbmc, color='louvain', edges=True, edges_width=0.1, s=50, show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_umap_with_edges.png', outfile.name, tolerance)
    assert res is None, res

    # test diffmap
    sc.tl.diffmap(pbmc)
    sc.pl.diffmap(pbmc, components='all', color=['CD3D'], show=False)
    pl.savefig(outfile.name, dpi=80)
    pl.close()

    res = compare_images(ROOT + '/master_diffmap.png', outfile.name, tolerance)
    assert res is None, res
