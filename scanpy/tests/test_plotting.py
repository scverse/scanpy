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


def test_heatmap_numeric_column():
    adata = sc.datasets.krumsiek11()
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


def test_dotplot_numeric_column():
    adata = sc.datasets.krumsiek11()
    outfile = NamedTemporaryFile(suffix='.png', prefix='scanpy_test_dotplot2_', delete=False)
    adata.obs['Gata2'] = adata.X[:, 0]
    sc.pl.dotplot(adata, adata.var_names, 'Gata2', use_raw=False,
                  num_categories=7, figsize=(7, 2.5))
    pl.savefig(outfile.name, dpi=80)
    res = compare_images(ROOT + '/master_dotplot2.png', outfile.name, tolerance)

    assert res is None, res

    os.remove(outfile.name)
