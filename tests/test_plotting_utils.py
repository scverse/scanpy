from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, cast

import numpy as np
import pytest
from anndata import AnnData
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

from scanpy.plotting._utils import DefaultProxy, _validate_palette

viridis = cast(ListedColormap, colormaps["viridis"])


@pytest.mark.parametrize(
    "palette",
    [
        pytest.param(viridis.colors, id="viridis"),
        pytest.param(["b", "#cccccc", "r", "yellow", "lightblue"], id="named"),
        pytest.param([(1, 0, 0, 1), (0, 0, 1, 1)], id="rgba"),
    ],
)
@pytest.mark.parametrize("typ", [np.asarray, list])
def test_validate_palette_no_mod(palette, typ):
    palette = typ(palette)
    adata = AnnData(uns=dict(test_colors=palette))
    _validate_palette(adata, "test")
    assert palette is adata.uns["test_colors"], "Palette should not be modified"


@pytest.mark.parametrize(
    "param",
    [
        pytest.param(1, id="direct"),
        pytest.param(field(default=1), id="default"),
        pytest.param(
            field(default_factory=lambda: 1),
            marks=[
                pytest.mark.xfail(
                    "Tries to call factory while class not fully constructed"
                )
            ],
            id="default_factory",
        ),
    ],
)
def test_default_proxy(param):
    @dataclass
    class Test:
        field_: int = param
        DEFAULT: ClassVar[DefaultProxy[int]] = DefaultProxy("field_")

    assert Test(2).field_ == 2
    assert Test(2).DEFAULT == 1
