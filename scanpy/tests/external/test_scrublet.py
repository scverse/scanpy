import pytest

import scanpy as sc
import scanpy.external as sce
from anndata.tests.helpers import assert_equal
import pandas as pd
import anndata as ad
from scanpy.get import _get_obs_rep
import numpy as np
import scanpy.preprocessing as pp
import scipy.sparse as sparse


def test_scrublet():
    """
    Test that Scrublet run works.

    Check that scrublet runs and detects some doublets.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.pbmc3k()
    sce.pp.scrublet(adata, use_approx_neighbors=False)

    # replace assertions by conditions
    assert "predicted_doublet" in adata.obs.columns
    assert "doublet_score" in adata.obs.columns

    assert adata.obs["predicted_doublet"].any(), "Expect some doublets to be identified"


def test_scrublet_batched():
    """
    Test that Scrublet run works with batched data.

    Check that scrublet runs and detects some doublets.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.pbmc3k()
    adata.obs['batch'] = 1350 * ['a'] + 1350 * ['b']
    split = [adata[adata.obs["batch"] == x].copy() for x in ("a", "b")]

    sce.pp.scrublet(adata, use_approx_neighbors=False, batch_key='batch')

    # replace assertions by conditions
    assert "predicted_doublet" in adata.obs.columns
    assert "doublet_score" in adata.obs.columns

    assert adata.obs["predicted_doublet"].any(), "Expect some doublets to be identified"
    assert (
        'batches' in adata.uns['scrublet'].keys()
    ), "Expect .uns to contain batch info"

    # Check that results are independent
    for s in split:
        sce.pp.scrublet(s, use_approx_neighbors=False)
    merged = sc.concat(split)

    pd.testing.assert_frame_equal(adata.obs[merged.obs.columns], merged.obs)


def test_scrublet_data():
    """
    Test that Scrublet processing is arranged correctly.

    Check that simulations run on raw data.
    """
    pytest.importorskip("scrublet")

    random_state = 1234

    # Run Scrublet and let the main function run simulations
    adata_scrublet_auto_sim = sce.pp.scrublet(
        sc.datasets.pbmc3k(),
        use_approx_neighbors=False,
        copy=True,
        random_state=random_state,
    )

    # Now make our own simulated data so we can check the result from function
    # is the same, and by inference that the processing steps have not been
    # broken

    # Replicate the preprocessing steps used by the main function

    def preprocess_for_scrublet(adata):

        adata_pp = adata.copy()
        pp.filter_genes(adata_pp, min_cells=3)
        pp.filter_cells(adata_pp, min_genes=3)
        adata_pp.layers['raw'] = adata_pp.X.copy()
        pp.normalize_total(adata_pp)
        logged = pp.log1p(adata_pp, copy=True)
        pp.highly_variable_genes(logged)
        adata_pp = adata_pp[:, logged.var['highly_variable']]

        return adata_pp

    # Simulate doublets using the same parents

    def create_sim_from_parents(adata, parents):

        # Now simulate doublets based on the randomly selected parents used
        # previously

        N_sim = parents.shape[0]
        I = sparse.coo_matrix(
            (
                np.ones(2 * N_sim),
                (np.repeat(np.arange(N_sim), 2), parents.flat),
            ),
            (N_sim, adata_obs.n_obs),
        )
        X = I @ adata_obs.layers['raw']
        return ad.AnnData(
            X,
            var=pd.DataFrame(index=adata_obs.var_names),
            obs=pd.DataFrame({"total_counts": np.ravel(X.sum(axis=1))}),
            obsm={"doublet_parents": parents.copy()},
        )

    # Preprocess the data and make the simulated doublets

    adata_obs = preprocess_for_scrublet(sc.datasets.pbmc3k())
    adata_sim = create_sim_from_parents(
        adata_obs, adata_scrublet_auto_sim.uns['scrublet']['doublet_parents']
    )

    # Apply the same post-normalisation the Scrublet function would

    pp.normalize_total(adata_obs, target_sum=1e6)
    pp.normalize_total(adata_sim, target_sum=1e6)

    adata_scrublet_manual_sim = sce.pp.scrublet(
        adata_obs,
        adata_sim=adata_sim,
        use_approx_neighbors=False,
        copy=True,
        random_state=random_state,
    )

    # Require that the doublet scores are the same whether simulation is via
    # the main function or manually provided

    assert (
        adata_scrublet_manual_sim.obs['doublet_score']
        == adata_scrublet_auto_sim.obs['doublet_score']
    ).all()


def test_scrublet_dense():
    """
    Test that Scrublet works for dense matrices.

    Check that scrublet runs and detects some doublets when a dense matrix is supplied.
    """
    pytest.importorskip("scrublet")

    adata = sc.datasets.paul15()[:500].copy()
    sce.pp.scrublet(adata, use_approx_neighbors=False)

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

    # Reduce size of input for faster test
    adata = sc.datasets.pbmc3k()[:500].copy()
    sc.pp.filter_genes(adata, min_counts=100)

    # Get the default output

    default = sce.pp.scrublet(adata, use_approx_neighbors=False, copy=True)

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
            'use_approx_neighbors': False,
            'copy': True,
            param: test_params[param],
        }
        curr = sc.external.pp.scrublet(**test_args)
        with pytest.raises(AssertionError):
            assert_equal(default, curr)


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

    _ = sc.pp.highly_variable_genes(logged)
    adata_obs = adata_obs[:, logged.var['highly_variable']]

    adata_sim = sce.pp.scrublet_simulate_doublets(adata_obs, layer='raw')

    assert 'doublet_parents' in adata_sim.obsm.keys()
