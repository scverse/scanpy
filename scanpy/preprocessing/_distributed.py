from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING, Literal, overload

import numpy as np

from scanpy._compat import DaskArray, DaskSeries, DaskSeriesGroupBy, ZappyArray

if TYPE_CHECKING:
    from collections.abc import Callable

    import pandas as pd
    from dask.dataframe import Aggregation
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
