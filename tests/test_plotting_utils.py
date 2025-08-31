from __future__ import annotations

from string import ascii_lowercase, ascii_uppercase
from typing import TYPE_CHECKING, cast

import numpy as np
import pytest
from anndata import AnnData
from matplotlib import colormaps

from scanpy.plotting._anndata import _check_if_annotations
from scanpy.plotting._utils import _validate_palette

if TYPE_CHECKING:
    from typing import Any, Literal

    from matplotlib.colors import ListedColormap


viridis = cast("ListedColormap", colormaps["viridis"])


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
    ("axis_name", "args", "expected"),
    [
        pytest.param("obs", {}, True, id="valid-nothing"),
        pytest.param("obs", dict(x="B", colors=["obs_a"]), True, id="valid-basic"),
        pytest.param("var", dict(colors=["A", "C", "obs_a"]), False, id="invalid-axis"),
        pytest.param("obs", dict(x="A"), True, id="valid-raw"),
        pytest.param("obs", dict(x="A", use_raw=False), False, id="invalid-noraw"),
        pytest.param("obs", dict(colors=[(0, 0, 0), "red"]), True, id="valid-color"),
    ],
)
def test_check_all_in_axis(
    *, axis_name: Literal["obs", "var"], args: dict[str, Any], expected: bool
):
    raw = AnnData(
        np.random.randn(10, 20),
        dict(obs_a=range(10), obs_names=list(ascii_lowercase[:10])),
        dict(var_a=range(20), var_names=list(ascii_uppercase[:20])),
    )
    adata = raw[:, 1:].copy()
    adata.raw = raw

    assert _check_if_annotations(adata, axis_name, **args) is expected
