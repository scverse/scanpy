from __future__ import annotations

from contextlib import contextmanager
from itertools import chain
from typing import TYPE_CHECKING, Literal, overload

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
    from collections.abc import Callable, Generator

    import pandas as pd
    from dask.dataframe import Aggregation
    from numpy.typing import ArrayLike
    from pandas.core.groupby.generic import DataFrameGroupBy, SeriesGroupBy


@overload
def series_to_array(s: pd.Series, *, dtype: np.dtype | None = None) -> np.ndarray:
    ...


@overload
def series_to_array(s: DaskSeries, *, dtype: np.dtype | None = None) -> DaskArray:
    ...


def series_to_array(
    s: pd.Series | DaskSeries, *, dtype: np.dtype | None = None
) -> np.ndarray | DaskArray:
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


try:
    import dask.dataframe as dd
except ImportError:

    def get_mad(dask: Literal[False]) -> Callable[[np.ndarray], np.ndarray]:
        from statsmodels.robust import mad

        return mad
else:

    def _mad1(chunks: DaskSeriesGroupBy):
        return chunks.apply(list)

    def _mad2(grouped: DaskSeriesGroupBy):
        def internal(c):
            if (c != c).all():
                return [np.nan]
            f = [_ for _ in c if _ == _]
            f = [_ if isinstance(_, list) else [_] for _ in f]
            return list(chain.from_iterable(f))

        return grouped.apply(internal)

    def _mad3(grouped: DaskSeriesGroupBy):
        from statsmodels.robust import mad

        return grouped.apply(lambda s: np.nan if len(s) == 0 else mad(s))

    mad_dask = dd.Aggregation("mad", chunk=_mad1, agg=_mad2, finalize=_mad3)

    def get_mad(dask: bool) -> Callable[[np.ndarray], np.ndarray] | Aggregation:
        from statsmodels.robust import mad

        return mad_dask if dask else mad
