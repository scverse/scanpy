from __future__ import annotations

from typing import TYPE_CHECKING, overload

import numpy as np

from scanpy._compat import DaskArray, DaskSeries, ZappyArray

if TYPE_CHECKING:
    import pandas as pd
    from numpy.typing import ArrayLike


@overload
def series_to_array(s: pd.Series, *, dtype: np.dtype | None = None) -> np.ndarray:
    ...


@overload
def series_to_array(s: DaskSeries, *, dtype: np.dtype | None = None) -> DaskArray:
    ...


def series_to_array(
    s: pd.Series | DaskSeries, *, dtype: np.dtype | None = None
) -> np.ndarray | DaskArray:
    """Convert Series to Array, keeping them in-memory or distributed."""
    if isinstance(s, DaskSeries):
        return (
            s.to_dask_array(True)
            if dtype is None
            else s.astype(dtype).to_dask_array(True)
        )
    return s.to_numpy() if dtype is None else s.to_numpy().astype(dtype, copy=False)


@overload
def materialize_as_ndarray(a: ArrayLike) -> np.ndarray:
    ...


@overload
def materialize_as_ndarray(a: tuple[ArrayLike]) -> tuple[np.ndarray]:
    ...


@overload
def materialize_as_ndarray(
    a: tuple[ArrayLike, ArrayLike],
) -> tuple[np.ndarray, np.ndarray]:
    ...


@overload
def materialize_as_ndarray(
    a: tuple[ArrayLike, ArrayLike, ArrayLike],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ...


def materialize_as_ndarray(
    a: ArrayLike | tuple[ArrayLike | ZappyArray | DaskArray, ...],
) -> tuple[np.ndarray] | np.ndarray:
    """Compute distributed arrays and convert them to numpy ndarrays."""
    if not isinstance(a, tuple):
        return np.asarray(a)

    if not any(isinstance(arr, DaskArray) for arr in a):
        return tuple(np.asarray(arr) for arr in a)

    import dask.array as da

    return da.compute(*a, sync=True)
