from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scanpy as sc

from scanpy.tests.test_plotting import ROOT, FIGS

MISSING_VALUES_ROOT = ROOT / "embedding-missing-values"
MISSING_VALUES_FIGS = FIGS / "embedding-missing-values"


@pytest.fixture(scope="module")
def adata():
    """A bit cute."""
    from matplotlib.image import imread
    from sklearn.datasets import make_blobs
    from sklearn.cluster import DBSCAN

    empty_pixel = np.array([1.0, 1.0, 1.0, 0]).reshape(1, 1, -1)
    image = imread(
        Path(sc.__file__).parent.parent / "docs/_static/img/Scanpy_Logo_RGB.png"
    )
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
            pd.value_counts(labels[order]).values,
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


@pytest.fixture(params=[sc.pl.pca, sc.pl.spatial])
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
    params=[(None, None), (0, 5), ("p15", "p90")],
    ids=["vbounds.default", "vbound.numbers", "vbound.percentile"],
)
def vbounds(request):
    return request.param


def test_missing_values_categorical(
    fixture_request,
    image_comparer,
    adata,
    plotfunc,
    na_color,
    na_in_legend,
    legend_loc,
    groupsfunc,
):
    save_and_compare_images = image_comparer(
        MISSING_VALUES_ROOT, MISSING_VALUES_FIGS, tol=15
    )
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
    fixture_request, image_comparer, adata, plotfunc, na_color, legend_loc, vbounds
):
    save_and_compare_images = image_comparer(
        MISSING_VALUES_ROOT, MISSING_VALUES_FIGS, tol=15
    )
    base_name = fixture_request.node.name

    # Passing through a dict so it's easier to use default values
    kwargs = {}
    kwargs["vmin"], kwargs["vmax"] = vbounds
    kwargs["legend_loc"] = legend_loc
    if na_color is not None:
        kwargs["na_color"] = na_color

    plotfunc(adata, color=["1", "1_missing"], **kwargs)

    save_and_compare_images(base_name)
