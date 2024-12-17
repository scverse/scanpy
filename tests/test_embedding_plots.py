from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.testing.compare import compare_images

import scanpy as sc
from testing.scanpy._helpers.data import pbmc3k_processed

if TYPE_CHECKING:
    from scanpy.plotting._utils import _LegendLoc


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
    params=["on data", "right margin", "lower center", None],
    ids=["legend.on_data", "legend.on_right", "legend.on_bottom", "legend.off"],
)
def legend_loc(request) -> _LegendLoc | None:
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
    fixture_request: pytest.FixtureRequest,
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
    *,
    fixture_request: pytest.FixtureRequest,
    image_comparer,
    adata,
    plotfunc,
    na_color,
    vbounds,
):
    save_and_compare_images = partial(image_comparer, MISSING_VALUES_ROOT, tol=15)

    base_name = fixture_request.node.name

    # Passing through a dict so it's easier to use default values
    kwargs = {}
    kwargs.update(vbounds)
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

    with pytest.raises(
        ValueError,
        match=r"Could not broadcast together arguments with shapes: \[2, 3, 1\]",
    ):
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

    with pytest.raises(
        ValueError,
        match=r"Could not broadcast together arguments with shapes: \[2, 1, 3\]",
    ):
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
