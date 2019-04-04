"""
Tests to make sure the example datasets load.
"""

import scanpy as sc
import pytest
from pathlib import Path

@pytest.fixture(scope="module")
def tmp_dataset_dir(tmpdir_factory):
    new_dir = Path(tmpdir_factory.mktemp("scanpy_data"))
    print(new_dir)
    old_dir = sc.settings.datasetdir
    sc.settings.datasetdir = new_dir  # Set up
    yield sc.settings.datasetdir
    sc.settings.datasetdir = old_dir  # Tear down


@pytest.mark.internet
def test_burczynski06(tmp_dataset_dir):
    sc.datasets.burczynski06()


@pytest.mark.internet
def test_krumsiek11(tmp_dataset_dir):
    sc.datasets.krumsiek11()


@pytest.mark.internet
def test_moignard15(tmp_dataset_dir):
    sc.datasets.moignard15()


@pytest.mark.internet
def test_paul15(tmp_dataset_dir):
    sc.datasets.paul15()


@pytest.mark.internet
def test_pbmc3k(tmp_dataset_dir):
    sc.datasets.pbmc3k()


@pytest.mark.internet
def test_ebi_expression_atlas(tmp_dataset_dir):
    adata = sc.datasets.ebi_expression_atlas("E-MTAB-4888")


def test_blobs():
    sc.datasets.blobs()


def test_toggleswitch():
    sc.datasets.toggleswitch()


def test_pbmc68k_reduced():
    sc.datasets.pbmc68k_reduced()
