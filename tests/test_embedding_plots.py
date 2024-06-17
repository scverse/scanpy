from __future__ import annotations

from functools import partial
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.testing.compare import compare_images

import scanpy as sc
from testing.scanpy._helpers.data import pbmc3k_processed

HERE: Path = Path(__file__).parent
ROOT = HERE / "_images"

MISSING_VALUES_ROOT = ROOT / "embedding-missing-values"


def check_images(pth1, pth2, *, tol):
    result = compare_images(pth1, pth2, tol=tol)
    assert result is None, result


@pytest.fixture(scope="module")
def adata():
    """A bit cute."""
    from matplotlib.image import imread
    from sklearn.cluster import DBSCAN
    from sklearn.datasets import make_blobs

    empty_pixel = np.array([1.0, 1.0, 1.0, 0]).reshape(1, 1, -1)
    image = imread(HERE.parent / "docs/_static/img/Scanpy_Logo_RGB.png")
    x, y = np.where(np.logical_and.reduce(~np.equal(image, empty_pixel), axis=2))

    # Just using to calculate the hex coords
    hexes = plt.hexbin(x, y, gridsize=(44, 100))
    counts = hexes.get_array()
    pixels = hexes.get_offsets()[counts != 0]
    plt.close()

    labels = DBSCAN(eps=20, min_samples=2).fit(pixels).labels_
    order = np.argsort(labels)
    adata = sc.AnnData(
        make_blobs(
            pd.Series(labels[order]).value_counts().values,
            n_features=20,
            shuffle=False,
            random_state=42,
        )[0],
        obs={"label": pd.Categorical(labels[order].astype(str))},
        obsm={"spatial": pixels[order, ::-1]},
        uns={
            "spatial": {
                "scanpy_img": {
                    "images": {"hires": image},
                    "scalefactors": {
                        "tissue_hires_scalef": 1,
                        "spot_diameter_fullres": 10,
                    },
                }
            }
        },
    )
    sc.pp.pca(adata)

    # Adding some missing values
    adata.obs["label_missing"] = adata.obs["label"].copy()
    adata.obs["label_missing"][::2] = np.nan

    adata.obs["1_missing"] = adata.obs_vector("1")
    adata.obs.loc[
        adata.obsm["spatial"][:, 0] < adata.obsm["spatial"][:, 0].mean(), "1_missing"
    ] = np.nan

    return adata


@pytest.fixture
def fixture_request(request):
    """Returns a Request object.

    Allows you to access names of parameterized tests from within a test.
    """
    return request


@pytest.fixture(
    params=[(0, 0, 0, 1), None],
    ids=["na_color.black_tup", "na_color.default"],
)
def na_color(request):
    return request.param


@pytest.fixture(params=[True, False], ids=["na_in_legend.True", "na_in_legend.False"])
def na_in_legend(request):
    return request.param


@pytest.fixture(
    params=[partial(sc.pl.pca, show=False), partial(sc.pl.spatial, show=False)],
    ids=["pca", "spatial"],
)
def plotfunc(request):
    return request.param


@pytest.fixture(
    params=["on data", "right margin", None],
    ids=["legend.on_data", "legend.on_right", "legend.off"],
)
def legend_loc(request):
    return request.param


@pytest.fixture(
    params=[lambda x: list(x.cat.categories[:3]), lambda x: []],
    ids=["groups.3", "groups.all"],
)
def groupsfunc(request):
    return request.param


@pytest.fixture(
    params=[
        pytest.param(
            {"vmin": None, "vmax": None, "vcenter": None, "norm": None},
            id="vbounds.default",
        ),
        pytest.param(
            {"vmin": 0, "vmax": 5, "vcenter": None, "norm": None}, id="vbounds.numbers"
        ),
        pytest.param(
            {"vmin": "p15", "vmax": "p90", "vcenter": None, "norm": None},
            id="vbounds.percentile",
        ),
        pytest.param(
            {"vmin": 0, "vmax": "p99", "vcenter": 0.1, "norm": None},
            id="vbounds.vcenter",
        ),
        pytest.param(
            {"vmin": None, "vmax": None, "vcenter": None, "norm": Normalize(0, 5)},
            id="vbounds.norm",
        ),
    ]
)
def vbounds(request):
    return request.param


def test_missing_values_categorical(
    *,
    fixture_request,
    image_comparer,
    adata,
    plotfunc,
    na_color,
    na_in_legend,
    legend_loc,
    groupsfunc,
):
    save_and_compare_images = partial(image_comparer, MISSING_VALUES_ROOT, tol=15)

    base_name = fixture_request.node.name

    # Passing through a dict so it's easier to use default values
    kwargs = {}
    kwargs["legend_loc"] = legend_loc
    kwargs["groups"] = groupsfunc(adata.obs["label"])
    if na_color is not None:
        kwargs["na_color"] = na_color
    kwargs["na_in_legend"] = na_in_legend

    plotfunc(adata, color=["label", "label_missing"], **kwargs)

    save_and_compare_images(base_name)


def test_missing_values_continuous(
    *, fixture_request, image_comparer, adata, plotfunc, na_color, legend_loc, vbounds
):
    save_and_compare_images = partial(image_comparer, MISSING_VALUES_ROOT, tol=15)

    base_name = fixture_request.node.name

    # Passing through a dict so it's easier to use default values
    kwargs = {}
    kwargs.update(vbounds)
    kwargs["legend_loc"] = legend_loc
    if na_color is not None:
        kwargs["na_color"] = na_color

    plotfunc(adata, color=["1", "1_missing"], **kwargs)

    save_and_compare_images(base_name)


def test_enumerated_palettes(fixture_request, adata, tmpdir, plotfunc):
    tmpdir = Path(tmpdir)
    base_name = fixture_request.node.name

    categories = adata.obs["label"].cat.categories
    colors_rgb = dict(zip(categories, sns.color_palette(n_colors=12)))

    dict_pth = tmpdir / f"rgbdict_{base_name}.png"
    list_pth = tmpdir / f"rgblist_{base_name}.png"

    # making a copy so colors aren't saved
    plotfunc(adata.copy(), color="label", palette=colors_rgb)
    plt.savefig(dict_pth, dpi=40)
    plt.close()
    plotfunc(adata.copy(), color="label", palette=[colors_rgb[c] for c in categories])
    plt.savefig(list_pth, dpi=40)
    plt.close()

    check_images(dict_pth, list_pth, tol=15)


def test_dimension_broadcasting(adata, tmpdir, check_same_image):
    tmpdir = Path(tmpdir)

    with pytest.raises(ValueError):
        sc.pl.pca(
            adata, color=["label", "1_missing"], dimensions=[(0, 1), (1, 2), (2, 3)]
        )

    dims_pth = tmpdir / "broadcast_dims.png"
    color_pth = tmpdir / "broadcast_colors.png"

    sc.pl.pca(adata, color=["label", "label", "label"], dimensions=(2, 3), show=False)
    plt.savefig(dims_pth, dpi=40)
    plt.close()
    sc.pl.pca(adata, color="label", dimensions=[(2, 3), (2, 3), (2, 3)], show=False)
    plt.savefig(color_pth, dpi=40)
    plt.close()

    check_same_image(dims_pth, color_pth, tol=5)


def test_marker_broadcasting(adata, tmpdir, check_same_image):
    tmpdir = Path(tmpdir)

    with pytest.raises(ValueError):
        sc.pl.pca(adata, color=["label", "1_missing"], marker=[".", "^", "x"])

    dims_pth = tmpdir / "broadcast_markers.png"
    color_pth = tmpdir / "broadcast_colors_for_markers.png"

    sc.pl.pca(adata, color=["label", "label", "label"], marker="^", show=False)
    plt.savefig(dims_pth, dpi=40)
    plt.close()
    sc.pl.pca(adata, color="label", marker=["^", "^", "^"], show=False)
    plt.savefig(color_pth, dpi=40)
    plt.close()

    check_same_image(dims_pth, color_pth, tol=5)


def test_dimensions_same_as_components(adata, tmpdir, check_same_image):
    tmpdir = Path(tmpdir)
    adata = adata.copy()
    adata.obs["mean"] = np.ravel(adata.X.mean(axis=1))

    comp_pth = tmpdir / "components_plot.png"
    dims_pth = tmpdir / "dimension_plot.png"

    # TODO: Deprecate components kwarg
    # with pytest.warns(FutureWarning, match=r"components .* deprecated"):
    sc.pl.pca(
        adata,
        color=["mean", "label"],
        components=["1,2", "2,3"],
        show=False,
    )
    plt.savefig(comp_pth, dpi=40)
    plt.close()

    sc.pl.pca(
        adata,
        color=["mean", "mean", "label", "label"],
        dimensions=[(0, 1), (1, 2), (0, 1), (1, 2)],
        show=False,
    )
    plt.savefig(dims_pth, dpi=40)
    plt.close()

    check_same_image(dims_pth, comp_pth, tol=5)


def test_embedding_colorbar_location(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = pbmc3k_processed().raw.to_adata()

    sc.pl.pca(adata, color="LDHB", colorbar_loc=None)

    save_and_compare_images("no_colorbar")


# Spatial specific


def test_visium_circles(image_comparer):  # standard visium data
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(HERE / "_data" / "visium_data" / "1.0.0")
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

    adata = sc.read_visium(HERE / "_data" / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    # Points default to transparent if an image is included
    sc.pl.spatial(adata, show=False)

    save_and_compare_images("spatial_visium_default")


def test_visium_empty_img_key(image_comparer):  # visium coordinates but image empty
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(HERE / "_data" / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})

    sc.pl.spatial(adata, img_key=None, color="array_row", show=False)

    save_and_compare_images("spatial_visium_empty_image")

    sc.pl.embedding(adata, basis="spatial", color="array_row", show=False)
    save_and_compare_images("spatial_visium_embedding")


def test_spatial_general(image_comparer):  # general coordinates
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = sc.read_visium(HERE / "_data" / "visium_data" / "1.0.0")
    adata.obs = adata.obs.astype({"array_row": "str"})
    spatial_metadata = adata.uns.pop(
        "spatial"
    )  # spatial data don't have imgs, so remove entry from uns
    # Required argument for now
    spot_size = list(spatial_metadata.values())[0]["scalefactors"][
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

    adata = sc.read_visium(HERE / "_data" / "visium_data" / "1.0.0")
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
    library_id = list(adata.uns["spatial"])[0]
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


def test_manual_equivalency(equivalent_spatial_plotters, tmpdir, spatial_kwargs):
    """
    Tests that manually passing values to sc.pl.spatial is similar to automatic extraction.
    """
    orig, removed = equivalent_spatial_plotters

    TESTDIR = Path(tmpdir)
    orig_pth = TESTDIR / "orig.png"
    removed_pth = TESTDIR / "removed.png"

    orig(**spatial_kwargs)
    plt.savefig(orig_pth, dpi=40)
    plt.close()
    removed(**spatial_kwargs)
    plt.savefig(removed_pth, dpi=40)
    plt.close()

    check_images(orig_pth, removed_pth, tol=1)


def test_manual_equivalency_no_img(
    equivalent_spatial_plotters_no_img, tmpdir, spatial_kwargs
):
    if "bw" in spatial_kwargs:
        # Has no meaning when there is no image
        pytest.skip()
    orig, removed = equivalent_spatial_plotters_no_img

    TESTDIR = Path(tmpdir)
    orig_pth = TESTDIR / "orig.png"
    removed_pth = TESTDIR / "removed.png"

    orig(**spatial_kwargs)
    plt.savefig(orig_pth, dpi=40)
    plt.close()
    removed(**spatial_kwargs)
    plt.savefig(removed_pth, dpi=40)
    plt.close()

    check_images(orig_pth, removed_pth, tol=1)


def test_white_background_vs_no_img(adata, tmpdir, spatial_kwargs):
    if {"bw", "img", "img_key", "na_color"}.intersection(spatial_kwargs):
        # These arguments don't make sense for this check
        pytest.skip()

    white_background = np.ones_like(
        adata.uns["spatial"]["scanpy_img"]["images"]["hires"]
    )
    TESTDIR = Path(tmpdir)
    white_pth = TESTDIR / "white_background.png"
    noimg_pth = TESTDIR / "no_img.png"

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


def test_spatial_na_color(adata, tmpdir):
    """
    Check that na_color defaults to transparent when an image is present, light gray when not.
    """
    white_background = np.ones_like(
        adata.uns["spatial"]["scanpy_img"]["images"]["hires"]
    )
    TESTDIR = Path(tmpdir)
    lightgray_pth = TESTDIR / "lightgray.png"
    transparent_pth = TESTDIR / "transparent.png"
    noimg_pth = TESTDIR / "noimg.png"
    whiteimg_pth = TESTDIR / "whiteimg.png"

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
