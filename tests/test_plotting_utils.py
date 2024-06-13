from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from anndata import AnnData
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

from scanpy.plotting._utils import _validate_palette

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
