from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import scanpy as sc
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import NDArray


@pytest.mark.parametrize(
    ("vec", "expected", "kwargs"),
    [
        pytest.param(
            np.array([0.0, 0.5, 1.0]),
            np.array([0.0, 0.5**1.5 * 20, 20.0]),
            dict(),
            id="default",
        ),
        pytest.param(
            np.array([0.0, 1.0]),
            np.array([2.0, 10.0]),
            dict(smallest_dot=2, largest_dot=10),
            id="smallest_largest",
        ),
        pytest.param(
            np.array([0.1, 0.5, 1.0]),
            # 0.1 is clipped up to dot_min, 1.0 is already at dot_max,
            # 0.5 lands 3/7 of the way in between
            np.array([0.0, (3 / 7) ** 1.5 * 20, 20.0]),
            dict(dot_min=0.2, dot_max=0.9),
            id="clips_and_rescales",
        ),
    ],
)
def test_dot_area(
    vec: NDArray[np.floating], expected: NDArray[np.floating], kwargs: dict[str, Any]
) -> None:
    size = sc.pl.dot_area(np.array(vec), **kwargs)
    np.testing.assert_allclose(size, expected)


@needs.scanpy2
@pytest.mark.parametrize(("dot_min", "dot_max"), [(0, 1), (0.2, 0.9)])
def test_dot_area_matches_hv_dim(dot_min: float, dot_max: float) -> None:
    import holoviews as hv

    vec = np.array([0.1, 0.5, 1.0])
    ds = hv.Dataset(vec, kdims=["x"])

    arr_result = sc.pl.dot_area(vec, dot_min=dot_min, dot_max=dot_max)
    dim_result = sc.pl.dot_area(hv.dim("x"), dot_min=dot_min, dot_max=dot_max).apply(ds)
    np.testing.assert_allclose(arr_result, dim_result)
