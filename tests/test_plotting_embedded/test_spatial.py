from __future__ import annotations

from functools import partial
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.testing.compare import compare_images

import scanpy as sc

HERE: Path = Path(__file__).parent
ROOT = HERE.parent / "_images"
DATA_DIR = HERE.parent / "_data"


pytestmark = [
    pytest.mark.filterwarnings("ignore:Use `squidpy.*` instead:FutureWarning")
]


def check_images(pth1: Path, pth2: Path, *, tol: int) -> None:
    result = compare_images(str(pth1), str(pth2), tol=tol)
    assert result is None, result


def test_visium_circles(image_comparer):  # standard visium data
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(DATA_DIR / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    sc.pl.spatial(
        adata,
        color="array_row",
        groups=["24", "33"],
        crop_coord=(100, 400, 400, 100),
        alpha=0.5,
        size=1.3,
        show=False,
    )

    save_and_compare_images("spatial_visium")


def test_visium_default(image_comparer):  # default values
    from packaging.version import parse as parse_version

    if parse_version(mpl.__version__) < parse_version("3.7.0"):
        pytest.xfail("Matplotlib 3.7.0+ required for this test")

    save_and_compare_images = partial(image_comparer, ROOT, tol=5)

    adata = sc.read_visium(DATA_DIR / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    # Points default to transparent if an image is included
    sc.pl.spatial(adata, show=False)

    save_and_compare_images("spatial_visium_default")


def test_visium_empty_img_key(image_comparer):  # visium coordinates but image empty
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(DATA_DIR / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    sc.pl.spatial(adata, img_key=None, color="array_row", show=False)

    save_and_compare_images("spatial_visium_empty_image")

    sc.pl.embedding(adata, basis="spatial", color="array_row", show=False)
    save_and_compare_images("spatial_visium_embedding")


def test_spatial_general(image_comparer):  # general coordinates
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(DATA_DIR / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})
    spatial_metadata = adata.uns.pop(
        "spatial"
    )  # spatial data don't have imgs, so remove entry from uns
    # Required argument for now
    spot_size = next(iter(spatial_metadata.values()))["scalefactors"][
        "spot_diameter_fullres"
    ]

    sc.pl.spatial(adata, show=False, spot_size=spot_size)
    save_and_compare_images("spatial_general_nocol")

    # category
    sc.pl.spatial(adata, show=False, spot_size=spot_size, color="array_row")
    save_and_compare_images("spatial_general_cat")

    # continuous
    sc.pl.spatial(adata, show=False, spot_size=spot_size, color="array_col")
    save_and_compare_images("spatial_general_cont")


def test_spatial_external_img(image_comparer):  # external image
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(DATA_DIR / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    img = adata.uns["spatial"]["custom"]["images"]["hires"]
    scalef = adata.uns["spatial"]["custom"]["scalefactors"]["tissue_hires_scalef"]
    sc.pl.spatial(
        adata,
        color="array_row",
        scale_factor=scalef,
        img=img,
        basis="spatial",
        show=False,
    )
    save_and_compare_images("spatial_external_img")


@pytest.fixture(scope="module")
def equivalent_spatial_plotters(adata):
    no_spatial = adata.copy()
    del no_spatial.uns["spatial"]

    img_key = "hires"
    library_id = next(iter(adata.uns["spatial"]))
    spatial_data = adata.uns["spatial"][library_id]
    img = spatial_data["images"][img_key]
    scale_factor = spatial_data["scalefactors"][f"tissue_{img_key}_scalef"]
    spot_size = spatial_data["scalefactors"]["spot_diameter_fullres"]

    orig_plotter = partial(sc.pl.spatial, adata, color="1", show=False)
    removed_plotter = partial(
        sc.pl.spatial,
        no_spatial,
        color="1",
        img=img,
        scale_factor=scale_factor,
        spot_size=spot_size,
        show=False,
    )

    return (orig_plotter, removed_plotter)


@pytest.fixture(scope="module")
def equivalent_spatial_plotters_no_img(equivalent_spatial_plotters):
    orig, removed = equivalent_spatial_plotters
    return (partial(orig, img_key=None), partial(removed, img=None, scale_factor=None))


@pytest.fixture(
    params=[
        pytest.param({"crop_coord": (50, 200, 0, 500)}, id="crop"),
        pytest.param({"size": 0.5}, id="size:.5"),
        pytest.param({"size": 2}, id="size:2"),
        pytest.param({"spot_size": 5}, id="spotsize"),
        pytest.param({"bw": True}, id="bw"),
        # Shape of the image for particular fixture, should not be hardcoded like this
        pytest.param({"img": np.ones((774, 1755, 4)), "scale_factor": 1.0}, id="img"),
        pytest.param(
            {"na_color": (0, 0, 0, 0), "color": "1_missing"}, id="na_color.transparent"
        ),
        pytest.param(
            {"na_color": "lightgray", "color": "1_missing"}, id="na_color.lightgray"
        ),
    ]
)
def spatial_kwargs(request):
    return request.param


def test_manual_equivalency(equivalent_spatial_plotters, tmp_path, spatial_kwargs):
    """Tests that manually passing values to sc.pl.spatial is similar to automatic extraction."""
    orig, removed = equivalent_spatial_plotters

    orig_pth = tmp_path / "orig.png"
    removed_pth = tmp_path / "removed.png"

    orig(**spatial_kwargs)
    plt.savefig(orig_pth, dpi=40)
    plt.close()
    removed(**spatial_kwargs)
    plt.savefig(removed_pth, dpi=40)
    plt.close()

    check_images(orig_pth, removed_pth, tol=1)


def test_manual_equivalency_no_img(
    equivalent_spatial_plotters_no_img, tmp_path, spatial_kwargs
):
    if "bw" in spatial_kwargs:
        # Has no meaning when there is no image
        pytest.skip()
    orig, removed = equivalent_spatial_plotters_no_img

    orig_pth = tmp_path / "orig.png"
    removed_pth = tmp_path / "removed.png"

    orig(**spatial_kwargs)
    plt.savefig(orig_pth, dpi=40)
    plt.close()
    removed(**spatial_kwargs)
    plt.savefig(removed_pth, dpi=40)
    plt.close()

    check_images(orig_pth, removed_pth, tol=1)


def test_white_background_vs_no_img(adata, tmp_path, spatial_kwargs):
    if {"bw", "img", "img_key", "na_color"}.intersection(spatial_kwargs):
        # These arguments don't make sense for this check
        pytest.skip()

    white_background = np.ones_like(
        adata.uns["spatial"]["scanpy_img"]["images"]["hires"]
    )
    white_pth = tmp_path / "white_background.png"
    noimg_pth = tmp_path / "no_img.png"

    sc.pl.spatial(
        adata,
        color="2",
        img=white_background,
        scale_factor=1.0,
        show=False,
        **spatial_kwargs,
    )
    plt.savefig(white_pth)
    sc.pl.spatial(adata, color="2", img_key=None, show=False, **spatial_kwargs)
    plt.savefig(noimg_pth)

    check_images(white_pth, noimg_pth, tol=1)


def test_spatial_na_color(adata, tmp_path):
    """Check that na_color defaults to transparent when an image is present, light gray when not."""
    white_background = np.ones_like(
        adata.uns["spatial"]["scanpy_img"]["images"]["hires"]
    )
    lightgray_pth = tmp_path / "lightgray.png"
    transparent_pth = tmp_path / "transparent.png"
    noimg_pth = tmp_path / "noimg.png"
    whiteimg_pth = tmp_path / "whiteimg.png"

    def plot(pth, **kwargs):
        sc.pl.spatial(adata, color="1_missing", show=False, **kwargs)
        plt.savefig(pth, dpi=40)
        plt.close()

    plot(lightgray_pth, na_color="lightgray", img_key=None)
    plot(transparent_pth, na_color=(0.0, 0.0, 0.0, 0.0), img_key=None)
    plot(noimg_pth, img_key=None)
    plot(whiteimg_pth, img=white_background, scale_factor=1.0)

    check_images(lightgray_pth, noimg_pth, tol=1)
    check_images(transparent_pth, whiteimg_pth, tol=1)
    with pytest.raises(AssertionError):
        check_images(lightgray_pth, transparent_pth, tol=1)
