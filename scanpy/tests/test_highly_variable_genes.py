import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path

FILE = Path(__file__).parent / Path('_scripts/seurat_hvg.csv')


def test_highly_variable_genes_basic():
    adata = sc.datasets.blobs()
    sc.pp.highly_variable_genes(adata)

    adata = sc.datasets.blobs()
    np.random.seed(0)
    adata.obs['batch'] = np.random.binomial(3, 0.5, size=(adata.n_obs))
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    sc.pp.highly_variable_genes(adata, batch_key='batch')
    assert 'highly_variable_nbatches' in adata.var.columns
    assert 'highly_variable_intersection' in adata.var.columns

    adata = sc.datasets.blobs()
    adata.obs['batch'] = np.random.binomial(4, 0.5, size=(adata.n_obs))
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    sc.pp.highly_variable_genes(adata, batch_key='batch', n_top_genes=3)
    assert 'highly_variable_nbatches' in adata.var.columns
    assert adata.var['highly_variable'].sum() == 3

def test_higly_variable_genes_compare_to_seurat():
    seurat_hvg_info = pd.read_csv(FILE, sep=' ')

    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    pbmc.var_names_make_unique()

    sc.pp.normalize_per_cell(pbmc, counts_per_cell_after=1e4)
    sc.pp.log1p(pbmc)
    sc.pp.highly_variable_genes(pbmc, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True)

    np.testing.assert_array_equal(seurat_hvg_info['highly_variable'], pbmc.var['highly_variable'])

    # (still) Not equal to tolerance rtol=2e-05, atol=2e-05
    # np.testing.assert_allclose(4, 3.9999, rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(
        seurat_hvg_info['means'],
        pbmc.var['means'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions'],
        pbmc.var['dispersions'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions_norm'],
        pbmc.var['dispersions_norm'],
        rtol=2e-05,
        atol=2e-05,
    )


def test_filter_genes_dispersion_compare_to_seurat():
    seurat_hvg_info = pd.read_csv(FILE, sep=' ')

    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    pbmc.var_names_make_unique()

    sc.pp.normalize_per_cell(pbmc, counts_per_cell_after=1e4)
    sc.pp.filter_genes_dispersion(
        pbmc, flavor='seurat', log=True, subset=False,
        min_mean=0.0125, max_mean=3, min_disp=0.5,
    )

    np.testing.assert_array_equal(seurat_hvg_info['highly_variable'], pbmc.var['highly_variable'])

    # (still) Not equal to tolerance rtol=2e-05, atol=2e-05:
    # np.testing.assert_allclose(4, 3.9999, rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(
        seurat_hvg_info['means'],
        pbmc.var['means'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions'],
        pbmc.var['dispersions'],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info['dispersions_norm'],
        pbmc.var['dispersions_norm'],
        rtol=2e-05,
        atol=2e-05,
    )
