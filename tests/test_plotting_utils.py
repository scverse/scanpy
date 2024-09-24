from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, cast

import numpy as np
import pytest
from anndata import AnnData
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

from scanpy.plotting._utils import (
    ClassDescriptorEnabled,
    DefaultProxy,
    _validate_palette,
)

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
                    reason="Tries to call factory while class not fully constructed"
                )
            ],
            id="default_factory",
        ),
    ],
)
@pytest.mark.parametrize("set_", ["instance", "field_", "DEFAULT"])
def test_default_proxy(param, set_: str):
    @dataclass
    class Test(metaclass=ClassDescriptorEnabled):
        field_: int = param
        DEFAULT: ClassVar[DefaultProxy[int]] = DefaultProxy("field_")

    instance = Test(2)
    assert instance.field_ == 2
    # instantiating doesn’t update the class
    assert instance.DEFAULT == Test().field_ == 1

    instance.field_ = 3
    # updating the instance doesn’t update the class
    assert Test.field_ == Test.DEFAULT == 1

    if set_ == "instance":
        v = 1
    elif set_ == "field_":
        Test.field_ = v = 4
    elif set_ == "DEFAULT":
        with pytest.warns(FutureWarning):
            Test.DEFAULT = v = 5
    else:
        pytest.fail(f"Unknown {set_=}")

    # updating anything doesn’t update existing instances
    assert instance.field_ == 3
    # setting the fields updates the class, but …
    assert Test.field_ == Test.DEFAULT == v
    # … sadly doesn’t update the __init__ method
    assert Test().field_ == 1
