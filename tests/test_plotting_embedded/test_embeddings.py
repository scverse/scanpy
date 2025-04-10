from __future__ import annotations

from functools import partial, wraps
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pytest
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.testing.compare import compare_images

import scanpy as sc
from testing.scanpy._helpers.data import pbmc3k_processed

if TYPE_CHECKING:
    from scanpy.plotting._utils import _LegendLoc


HERE: Path = Path(__file__).parent
ROOT = HERE.parent / "_images"

MISSING_VALUES_ROOT = ROOT / "embedding-missing-values"


def check_images(pth1: Path, pth2: Path, *, tol: int) -> None:
    result = compare_images(str(pth1), str(pth2), tol=tol)
    assert result is None, result


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
    if request.param is sc.pl.spatial:

        @wraps(request.param)
        def f(adata, **kwargs):
            with pytest.warns(FutureWarning, match=r"Use `squidpy.*` instead"):
                return sc.pl.spatial(adata, **kwargs)

    else:
        f = request.param
    return partial(f, show=False)


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
    request: pytest.FixtureRequest,
    image_comparer,
    adata,
    plotfunc,
    na_color,
    na_in_legend,
    legend_loc,
    groupsfunc,
):
    save_and_compare_images = partial(image_comparer, MISSING_VALUES_ROOT, tol=15)

    base_name = request.node.name

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
    request: pytest.FixtureRequest,
    image_comparer,
    adata,
    plotfunc,
    na_color,
    vbounds,
):
    save_and_compare_images = partial(image_comparer, MISSING_VALUES_ROOT, tol=15)

    base_name = request.node.name

    # Passing through a dict so it's easier to use default values
    kwargs = {}
    kwargs.update(vbounds)
    if na_color is not None:
        kwargs["na_color"] = na_color

    plotfunc(adata, color=["1", "1_missing"], **kwargs)

    save_and_compare_images(base_name)


def test_enumerated_palettes(request, adata, tmp_path, plotfunc):
    base_name = request.node.name

    categories = adata.obs["label"].cat.categories
    colors_rgb = dict(zip(categories, sns.color_palette(n_colors=12), strict=False))

    dict_pth = tmp_path / f"rgbdict_{base_name}.png"
    list_pth = tmp_path / f"rgblist_{base_name}.png"

    # making a copy so colors aren't saved
    plotfunc(adata.copy(), color="label", palette=colors_rgb)
    plt.savefig(dict_pth, dpi=40)
    plt.close()
    plotfunc(adata.copy(), color="label", palette=[colors_rgb[c] for c in categories])
    plt.savefig(list_pth, dpi=40)
    plt.close()

    check_images(dict_pth, list_pth, tol=15)


def test_dimension_broadcasting(adata, tmp_path, check_same_image):
    with pytest.raises(
        ValueError,
        match=r"Could not broadcast together arguments with shapes: \[2, 3, 1\]",
    ):
        sc.pl.pca(
            adata, color=["label", "1_missing"], dimensions=[(0, 1), (1, 2), (2, 3)]
        )

    dims_pth = tmp_path / "broadcast_dims.png"
    color_pth = tmp_path / "broadcast_colors.png"

    sc.pl.pca(adata, color=["label", "label", "label"], dimensions=(2, 3), show=False)
    plt.savefig(dims_pth, dpi=40)
    plt.close()
    sc.pl.pca(adata, color="label", dimensions=[(2, 3), (2, 3), (2, 3)], show=False)
    plt.savefig(color_pth, dpi=40)
    plt.close()

    check_same_image(dims_pth, color_pth, tol=5, root=tmp_path)


def test_marker_broadcasting(adata, tmp_path, check_same_image):
    with pytest.raises(
        ValueError,
        match=r"Could not broadcast together arguments with shapes: \[2, 1, 3\]",
    ):
        sc.pl.pca(adata, color=["label", "1_missing"], marker=[".", "^", "x"])

    dims_pth = tmp_path / "broadcast_markers.png"
    color_pth = tmp_path / "broadcast_colors_for_markers.png"

    sc.pl.pca(adata, color=["label", "label", "label"], marker="^", show=False)
    plt.savefig(dims_pth, dpi=40)
    plt.close()
    sc.pl.pca(adata, color="label", marker=["^", "^", "^"], show=False)
    plt.savefig(color_pth, dpi=40)
    plt.close()

    check_same_image(dims_pth, color_pth, tol=5, root=tmp_path)


def test_dimensions_same_as_components(adata, tmp_path, check_same_image):
    adata = adata.copy()
    adata.obs["mean"] = np.ravel(adata.X.mean(axis=1))

    comp_pth = tmp_path / "components_plot.png"
    dims_pth = tmp_path / "dimension_plot.png"

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

    check_same_image(dims_pth, comp_pth, tol=5, root=tmp_path)


def test_embedding_colorbar_location(image_comparer):
    save_and_compare_images = partial(image_comparer, ROOT, tol=15)

    adata = pbmc3k_processed().raw.to_adata()

    sc.pl.pca(adata, color="LDHB", colorbar_loc=None)

    save_and_compare_images("no_colorbar")
