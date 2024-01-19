from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING, overload

import numpy as np

from scanpy._compat import (
    DaskArray,
    DaskDataFrame,
    DaskDataFrameGroupBy,
    DaskSeries,
    DaskSeriesGroupBy,
    ZappyArray,
)

if TYPE_CHECKING:
    from collections.abc import Generator

    import pandas as pd
    from numpy.typing import ArrayLike
    from pandas.core.groupby.generic import DataFrameGroupBy, SeriesGroupBy


@overload
def series_to_array(s: pd.Series) -> np.ndarray:
    ...


@overload
def series_to_array(s: DaskSeries) -> DaskArray:
    ...


def series_to_array(s: pd.Series | DaskSeries) -> np.ndarray | DaskArray:
    return s.to_dask_array() if isinstance(s, DaskSeries) else s.to_numpy()


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
    """Convert distributed arrays to ndarrays."""
    if not isinstance(a, tuple):
        return np.asarray(a)

    if not any(isinstance(arr, DaskArray) for arr in a):
        return tuple(np.asarray(arr) for arr in a)

    import dask.array as da

    return da.compute(*a, sync=True)


@overload
def dask_compute(value: DaskDataFrame) -> pd.DataFrame:
    ...


@overload
def dask_compute(value: DaskSeries) -> pd.Series:
    ...


@overload
def dask_compute(value: DaskDataFrameGroupBy) -> DataFrameGroupBy:
    ...


@overload
def dask_compute(value: DaskSeriesGroupBy) -> SeriesGroupBy:
    ...


def dask_compute(
    value: DaskDataFrame | DaskSeries | DaskDataFrameGroupBy | DaskSeriesGroupBy,
) -> pd.DataFrame | pd.Series | DataFrameGroupBy | SeriesGroupBy:
    """Compute a dask array or series."""
    if isinstance(
        value, (DaskDataFrame, DaskSeries, DaskDataFrameGroupBy, DaskSeriesGroupBy)
    ):
        with suppress_pandas_warning():
            return value.compute(sync=True)
    return value


@contextmanager
def suppress_pandas_warning() -> Generator[None, None, None]:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"The default of observed=False", category=FutureWarning
        )
        yield
