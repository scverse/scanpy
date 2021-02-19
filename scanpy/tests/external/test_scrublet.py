import pytest

import scanpy as sc
import scanpy.external as sce
import io
from contextlib import redirect_stdout


def test_scrublet():
    """
    Test that Scrublet run works.

    Check that scrublet runs and detects some doublets.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.pbmc3k()
    sce.pp.scrublet(adata, use_approx_neighbors=False)

    errors = []

    # replace assertions by conditions
    assert "predicted_doublet" in adata.obs.columns
    assert "doublet_score" in adata.obs.columns

    assert adata.obs["predicted_doublet"].any(), "Expect some doublets to be identified"


def test_scrublet_params():
    """
    Test that Scrublet args are passed.

    Check that changes to parameters change scrublet results.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.pbmc3k()

    # Get the default output

    f = io.StringIO()
    with redirect_stdout(f):
        sce.pp.scrublet(adata, use_approx_neighbors=False)

    default = f.getvalue()

    test_params = {
        'expected_doublet_rate': 0.1,
        'synthetic_doublet_umi_subsampling': 0.8,
        'knn_dist_metric': 'manhattan',
        'normalize_variance': False,
        'log_transform': True,
        'mean_center': False,
        'n_prin_comps': 10,
        'n_neighbors': 2,
        'threshold': 0.1,
    }

    # Test each parameter and make sure something changes

    for param in test_params.keys():
        test_args = {
            'adata': adata,
            use_approx_neighbors: False,
            param: test_params[param],
        }
        f = io.StringIO()
        with redirect_stdout(f):
            sc.external.pp.scrublet(**test_args)
        diff_result = f.getvalue()

        assert diff_result != default


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
