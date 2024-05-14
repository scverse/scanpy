"""
Tests to make sure the example datasets load.
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata.tests.helpers import assert_adata_equal

import scanpy as sc

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    from anndata import AnnData


@pytest.fixture(scope="module")
def _tmp_dataset_dir(
    tmp_path_factory: pytest.TempPathFactory,
) -> Generator[None, None, None]:
    new_dir = tmp_path_factory.mktemp("scanpy_data")
    old_dir = sc.settings.datasetdir
    sc.settings.datasetdir = new_dir  # Set up
    yield
    sc.settings.datasetdir = old_dir  # Tear down


@pytest.mark.internet
def test_burczynski06(_tmp_dataset_dir):
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        adata = sc.datasets.burczynski06()
    assert adata.shape == (127, 22283)
    assert not (adata.X == 0).any()


@pytest.mark.internet
def test_moignard15(_tmp_dataset_dir):
    with warnings.catch_warnings():
        # https://foss.heptapod.net/openpyxl/openpyxl/-/issues/2051
        warnings.filterwarnings(
            "ignore",
            r"datetime\.datetime\.utcnow\(\) is deprecated",
            category=DeprecationWarning,
            module="openpyxl",
        )
        adata = sc.datasets.moignard15()
    assert adata.shape == (3934, 42)


@pytest.mark.internet
def test_paul15(_tmp_dataset_dir):
    sc.datasets.paul15()


@pytest.mark.internet
def test_pbmc3k(_tmp_dataset_dir):
    adata = sc.datasets.pbmc3k()
    assert adata.shape == (2700, 32738)
    assert "CD8A" in adata.var_names


@pytest.mark.internet
def test_pbmc3k_processed(_tmp_dataset_dir):
    with warnings.catch_warnings(record=True) as records:
        adata = sc.datasets.pbmc3k_processed()
    assert adata.shape == (2638, 1838)
    assert adata.raw.shape == (2638, 13714)

    assert len(records) == 0


@pytest.mark.internet
def test_ebi_expression_atlas(_tmp_dataset_dir):
    adata = sc.datasets.ebi_expression_atlas("E-MTAB-4888")
    # The shape changes sometimes
    assert 2261 <= adata.shape[0] <= 2315
    assert 23899 <= adata.shape[1] <= 24051


def test_krumsiek11(_tmp_dataset_dir):
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata = sc.datasets.krumsiek11()
    assert adata.shape == (640, 11)
    assert set(adata.obs["cell_type"]) == {"Ery", "Mk", "Mo", "Neu", "progenitor"}


def test_blobs():
    n_obs = np.random.randint(15, 30)
    n_var = np.random.randint(500, 600)
    adata = sc.datasets.blobs(n_variables=n_var, n_observations=n_obs)
    assert adata.shape == (n_obs, n_var)


def test_toggleswitch():
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        sc.datasets.toggleswitch()


def test_pbmc68k_reduced():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sc.datasets.pbmc68k_reduced()


@pytest.mark.internet
def test_visium_datasets(_tmp_dataset_dir, tmp_path: Path):
    # Tests that reading/ downloading works and is does not have global effects
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        hheart = sc.datasets.visium_sge("V1_Human_Heart")
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        mbrain = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        hheart_again = sc.datasets.visium_sge("V1_Human_Heart")
    assert_adata_equal(hheart, hheart_again)

    # Test that changing the dataset dir doesn't break reading
    sc.settings.datasetdir = tmp_path
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        mbrain_again = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    assert_adata_equal(mbrain, mbrain_again)

    # Test that downloading tissue image works
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        mbrain = sc.datasets.visium_sge("V1_Adult_Mouse_Brain", include_hires_tiff=True)
    expected_image_path = sc.settings.datasetdir / "V1_Adult_Mouse_Brain" / "image.tif"
    image_path = Path(
        mbrain.uns["spatial"]["V1_Adult_Mouse_Brain"]["metadata"]["source_image_path"]
    )
    assert image_path == expected_image_path

    # Test that tissue image exists and is a valid image file
    assert image_path.exists()

    # Test that tissue image is a tif image file (using `file`)
    process = subprocess.run(
        ["file", "--mime-type", image_path], stdout=subprocess.PIPE
    )
    output = process.stdout.strip().decode()  # make process output string
    assert output == str(image_path) + ": image/tiff"


def test_download_failure():
    from urllib.error import HTTPError

    with pytest.raises(HTTPError):
        sc.datasets.ebi_expression_atlas("not_a_real_accession")


# These are tested via doctest
DS_INCLUDED = frozenset({"krumsiek11", "toggleswitch", "pbmc68k_reduced"})
# These have parameters that affect shape and so on
DS_DYNAMIC = frozenset({"blobs", "ebi_expression_atlas"})


@pytest.mark.parametrize(
    "ds_name",
    [
        pytest.param(
            ds, id=ds, marks=[] if ds in DS_INCLUDED else [pytest.mark.internet]
        )
        for ds in set(sc.datasets.__all__) - DS_DYNAMIC
    ],
)
def test_doc_shape(ds_name):
    dataset_fn: Callable[[], AnnData] = getattr(sc.datasets, ds_name)
    assert dataset_fn.__doc__, "No docstring"
    docstring = dedent(dataset_fn.__doc__)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            r"(Observation|Variable) names are not unique",
            category=UserWarning,
        )
        dataset = dataset_fn()
    assert repr(dataset) in docstring
