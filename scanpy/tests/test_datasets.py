"""
Tests to make sure the example datasets load.
"""

import scanpy as sc
import pytest
from pathlib import Path

@pytest.fixture(scope="module")
def dataset_dir(tmpdir_factory):
    new_dir = Path(tmpdir_factory.mktemp("scanpy_data"))
    print(new_dir)
    old_dir = sc.settings.dataset_dir
    sc.settings.dataset_dir = new_dir  # Set up
    yield sc.settings.dataset_dir
    sc.settings.dataset_dir = old_dir  # Tear down


@pytest.mark.internet
def test_burczynski06(dataset_dir):
    sc.datasets.burczynski06()


@pytest.mark.internet
def test_krumsiek11(dataset_dir):
    sc.datasets.krumsiek11()


@pytest.mark.internet
def test_moignard15(dataset_dir):
    sc.datasets.moignard15()


@pytest.mark.internet
def test_paul15(dataset_dir):
    sc.datasets.paul15()


@pytest.mark.internet
def test_pbmc3k(dataset_dir):
    sc.datasets.pbmc3k()


def test_blobs():
    sc.datasets.blobs()


def test_toggleswitch():
    sc.datasets.toggleswitch()


def test_pbmc68k_reduced():
    sc.datasets.pbmc68k_reduced()
