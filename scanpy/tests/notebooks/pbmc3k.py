# coding: utf-8
# *First compiled on May 5, 2017. Updated August 14, 2018.*
# # Clustering 3k PBMCs following a Seurat Tutorial
#
# This started out with a demonstration that Scanpy would allow to reproduce most of Seurat's ([Satija *et al.*, 2015](https://doi.org/10.1038/nbt.3192)) clustering tutorial as described on http://satijalab.org/seurat/pbmc3k_tutorial.html (July 26, 2017), which we gratefully acknowledge. In the meanwhile, we have added and removed several pieces.
#
# The data consists in *3k PBMCs from a Healthy Donor* and is freely available from 10x Genomics ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)).


import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images

import numpy as np
import os
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, dpi_save=80, format='png')
sc.logging.print_versions()
sc.settings.autosave = True
sc.settings.autoshow = False
results_file = './write/pmbc3k.h5ad'

ROOT = os.path.dirname(os.path.abspath(__file__)) + '/pbmc3k_images/'

tolerance = 13  # default matplotlib pixel difference tolerance


def test_pbmc3k():

    adata = sc.read('./data/pbmc3k_raw.h5ad', backup_url='http://falexwolf.de/data/pbmc3k_raw.h5ad')

    # Preprocessing

    sc.pl.highest_expr_genes(adata, n_top=20)
    res = compare_images(ROOT + '/highest_expr_genes.png',
                         './figures/highest_expr_genes.png', tolerance * 2)
    assert res is None, res

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse to transform to a dense array after summing
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
