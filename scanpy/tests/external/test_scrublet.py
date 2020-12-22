import pytest

import scanpy as sc
import scanpy.external as sce


def test_scrublet():
    """
    Test that Scrublet run works.

    Check that scrublet runs and detects some doublets.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.pbmc3k()
    sce.pp.scrublet(adata)

    errors = []

    # replace assertions by conditions
    assert "predicted_doublet" in adata.obs.columns
    assert "doublet_score" in adata.obs.columns

    assert adata.obs["predicted_doublet"].any(), "Expect some doublets to be identified"


def test_scrublet_simulate_doublets():
    """
    Test that standalone Scrublet doublet simulation works.

    Check that doublet simulation runs and simulates some doublets..
    """
    pytest.importorskip("scrublet")

    adata_obs = sc.datasets.pbmc3k()
    sc.pp.filter_genes(adata_obs, min_cells=3)
    sc.pp.filter_cells(adata_obs, min_genes=3)
    adata_obs.layers['raw'] = adata_obs.X
    sc.pp.normalize_total(adata_obs)
    logged = sc.pp.log1p(adata_obs, copy=True)

    hvg = sc.pp.highly_variable_genes(logged)
    adata_obs = adata_obs[:, logged.var['highly_variable']]

    adata_sim = sce.pp.scrublet_simulate_doublets(adata_obs, layer='raw')

    assert 'doublet_parents' in adata_sim.obsm.keys()
