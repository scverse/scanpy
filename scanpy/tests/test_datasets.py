"""
Tests to make sure the example datasets load.
"""

import scanpy as sc
import numpy as np
import pytest
from pathlib import Path
from anndata.tests.helpers import assert_adata_equal


@pytest.fixture(scope="module")
def tmp_dataset_dir(tmpdir_factory):
    new_dir = Path(tmpdir_factory.mktemp("scanpy_data"))
    old_dir = sc.settings.datasetdir
    sc.settings.datasetdir = new_dir  # Set up
    yield sc.settings.datasetdir
    sc.settings.datasetdir = old_dir  # Tear down


@pytest.mark.internet
def test_burczynski06(tmp_dataset_dir):
    adata = sc.datasets.burczynski06()
    assert adata.shape == (127, 22283)
    assert not (adata.X == 0).any()


@pytest.mark.internet
def test_moignard15(tmp_dataset_dir):
    adata = sc.datasets.moignard15()
    assert adata.shape == (3934, 42)


@pytest.mark.internet
def test_paul15(tmp_dataset_dir):
    sc.datasets.paul15()


@pytest.mark.internet
def test_pbmc3k(tmp_dataset_dir):
    adata = sc.datasets.pbmc3k()
    assert adata.shape == (2700, 32738)
    assert "CD8A" in adata.var_names


@pytest.mark.internet
def test_ebi_expression_atlas(tmp_dataset_dir):
    adata = sc.datasets.ebi_expression_atlas("E-MTAB-4888")
    assert adata.shape == (2315, 23852)


def test_krumsiek11(tmp_dataset_dir):
    adata = sc.datasets.krumsiek11()
    assert adata.shape == (640, 11)
    assert all(
        np.unique(adata.obs["cell_type"])
        == np.array(["Ery", "Mk", "Mo", "Neu", "progenitor"])
    )


def test_blobs():
    n_obs = np.random.randint(15, 30)
    n_var = np.random.randint(500, 600)
    adata = sc.datasets.blobs(n_variables=n_var, n_observations=n_obs)
    assert adata.shape == (n_obs, n_var)


def test_toggleswitch():
    sc.datasets.toggleswitch()


def test_pbmc68k_reduced():
    sc.datasets.pbmc68k_reduced()


@pytest.mark.internet
def test_visium_datasets(tmp_dataset_dir, tmpdir):
    # Tests that reading/ downloading works and is does not have global effects
    hheart = sc.datasets.visium_sge("V1_Human_Heart")
    mbrain = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    hheart_again = sc.datasets.visium_sge("V1_Human_Heart")
    assert_adata_equal(hheart, hheart_again)

    # Test that changing the dataset dir doesn't break reading
    sc.settings.datasetdir = Path(tmpdir)
    mbrain_again = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    assert_adata_equal(mbrain, mbrain_again)


def test_download_failure():
    from urllib.error import HTTPError

    with pytest.raises(HTTPError):
        sc.datasets.ebi_expression_atlas("not_a_real_accession")
