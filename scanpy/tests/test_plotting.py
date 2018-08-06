# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl
from matplotlib.testing.compare import compare_images

import os.path
from tempfile import NamedTemporaryFile

import scanpy.api as sc

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

tolerance = 13  # default matplotlib pixel difference tolerance


def test_heatmap():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_heatmap_', delete=False)
    sc.pl.heatmap(adata, adata.var_names, 'cell_type', use_raw=False)
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_heatmap.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)

    # test heatmap numeric column():

    # set as numeric column the vales for the first gene on the matrix
    adata.obs['Gata2'] = adata.X[:, 0]
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_heatmap2_', delete=False)
    sc.pl.heatmap(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=4, figsize=(4.5, 5))
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_heatmap2.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_dotplot():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_dotplot_', delete=False)
    sc.pl.dotplot(adata, adata.var_names, 'cell_type', use_raw=False)
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_dotplot.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    # test dotplot numeric column():
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=7, figsize=(7, 2.5))
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_dotplot2.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)


def test_stacked_violin():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_stacked_violin_', delete=False)
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False)

    pl.title("image may have cut labels.\nThis is ok for test")
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_stacked_violin.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    # test swapped axes
    sc.pl.stacked_violin(adata, adata.var_names, 'cell_type', use_raw=False, swap_axes=True, figsize=(3, 5))
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_stacked_violin_swapped_axes.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)


def test_violin():
    pbmc = sc.datasets.pbmc68kb_reduced()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_violin_', delete=False)
    sc.pl.violin(pbmc, ['n_genes', 'percent_mito', 'n_counts'], stripplot=True, multi_panel=True, jitter=True)
    pl.savefig(outfile.name, dpi=80)

    res = compare_images(ROOT + '/master_violin_multi_panel.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)


def test_rank_genes_groups():
    pbmc = sc.datasets.pbmc68kb_reduced()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_rank_genes_groups_', delete=False)
    sc.pl.rank_genes_groups(pbmc, n_genes=12, n_panels_per_row=3)
    pl.savefig(outfile.name, dpi=80)

    res = compare_images(ROOT + '/master_ranked_genes.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=5)
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_ranked_genes_heatmap.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    sc.pl.rank_genes_groups_stacked_violin(pbmc, n_genes=3)
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_ranked_genes_stacked_violin.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)

    sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4)
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_ranked_genes_dotplot.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)
