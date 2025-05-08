"""Tests to make sure the example datasets load."""

from __future__ import annotations

import subprocess
import warnings
from collections import defaultdict
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata.tests.helpers import assert_adata_equal

import scanpy as sc
from testing.scanpy._helpers import data
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata import AnnData


@pytest.fixture(autouse=True)
def _tmp_dataset_dir(tmp_path: Path) -> None:
    """Make sure that datasets are downloaded during the test run.

    The default test environment stores them in a cached location.
    """
    sc.settings.datasetdir = tmp_path / "scanpy_data"


@pytest.mark.internet
def test_burczynski06():
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        adata = sc.datasets.burczynski06()
    assert adata.shape == (127, 22283)
    assert not (adata.X == 0).any()


@pytest.mark.internet
@needs.openpyxl
def test_moignard15():
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
def test_paul15():
    sc.datasets.paul15()


@pytest.mark.internet
def test_pbmc3k():
    adata = sc.datasets.pbmc3k()
    assert adata.shape == (2700, 32738)
    assert "CD8A" in adata.var_names


@pytest.mark.internet
def test_pbmc3k_processed():
    with warnings.catch_warnings(record=True) as records:
        adata = sc.datasets.pbmc3k_processed()
    assert adata.shape == (2638, 1838)
    assert adata.raw.shape == (2638, 13714)

    assert len(records) == 0


@pytest.mark.internet
def test_ebi_expression_atlas(monkeypatch: pytest.MonkeyPatch):
    from scanpy.datasets import _ebi_expression_atlas as ea_mod

    # make sure we use chunks when testing.
    # This dataset has <8M entries, so 4M entries/chunk = 2 chunks
    assert hasattr(ea_mod, "CHUNK_SIZE")
    monkeypatch.setattr(ea_mod, "CHUNK_SIZE", int(4e6))

    adata = sc.datasets.ebi_expression_atlas("E-MTAB-4888")
    # The shape changes sometimes
    assert 2261 <= adata.shape[0] <= 2315
    assert 23899 <= adata.shape[1] <= 24051


def test_krumsiek11():
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


@pytest.mark.filterwarnings("ignore:Use `squidpy.*` instead:FutureWarning")
@pytest.mark.internet
def test_visium_datasets():
    """Tests that reading/ downloading works and is does not have global effects."""
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        hheart = sc.datasets.visium_sge("V1_Human_Heart")
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        hheart_again = sc.datasets.visium_sge("V1_Human_Heart")
    assert_adata_equal(hheart, hheart_again)


@pytest.mark.filterwarnings("ignore:Use `squidpy.*` instead:FutureWarning")
@pytest.mark.internet
def test_visium_datasets_dir_change(tmp_path: Path):
    """Test that changing the dataset dir doesn't break reading."""
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        mbrain = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    sc.settings.datasetdir = tmp_path
    with pytest.warns(UserWarning, match=r"Variable names are not unique"):
        mbrain_again = sc.datasets.visium_sge("V1_Adult_Mouse_Brain")
    assert_adata_equal(mbrain, mbrain_again)


@pytest.mark.filterwarnings("ignore:Use `squidpy.*` instead:FutureWarning")
@pytest.mark.internet
def test_visium_datasets_images():
    """Test that image download works and is does not have global effects."""
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
        ["file", "--mime-type", image_path], stdout=subprocess.PIPE, check=True
    )
    output = process.stdout.strip().decode()  # make process output string
    assert output == f"{image_path}: image/tiff"


def test_download_failure():
    from urllib.error import HTTPError

    with pytest.raises(HTTPError):
        sc.datasets.ebi_expression_atlas("not_a_real_accession")


# These are tested via doctest
DS_INCLUDED = frozenset({"krumsiek11", "toggleswitch", "pbmc68k_reduced"})
# These have parameters that affect shape and so on
DS_DYNAMIC = frozenset({"ebi_expression_atlas"})
# Additional marks for datasets besides “internet”
DS_MARKS = defaultdict(list, moignard15=[needs.openpyxl])


@pytest.mark.parametrize(
    "ds_name",
    [
        pytest.param(
            ds,
            id=ds,
            marks=[
                *(() if ds in DS_INCLUDED else [pytest.mark.internet]),
                *DS_MARKS[ds],
            ],
        )
        for ds in sorted(set(sc.datasets.__all__) - DS_DYNAMIC)
    ],
)
def test_doc_shape(ds_name):
    dataset_fn: Callable[[], AnnData] = getattr(sc.datasets, ds_name)
    assert dataset_fn.__doc__, "No docstring"
    start_line_2 = dataset_fn.__doc__.find("\n") + 1
    docstring = dedent(dataset_fn.__doc__[start_line_2:])
    cached_fn = getattr(data, ds_name, dataset_fn)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            r"(Observation|Variable) names are not unique",
            category=UserWarning,
        )
        dataset = cached_fn()
    assert repr(dataset) in docstring
