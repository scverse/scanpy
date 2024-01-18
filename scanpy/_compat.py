from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial
from pathlib import Path

from legacy_api_wrap import legacy_api
from packaging import version

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    cache = lru_cache(maxsize=None)

try:
    from dask.array import Array as DaskArray
    from dask.dataframe import DataFrame as DaskDataFrame
    from dask.dataframe import Series as DaskSeries
    from dask.dataframe.groupby import DataFrameGroupBy as DaskDataFrameGroupBy
    from dask.dataframe.groupby import SeriesGroupBy as DaskSeriesGroupBy
except ImportError:

    class DaskArray:
        pass

    class DaskDataFrame:
        pass

    class DaskSeries:
        pass

    class DaskDataFrameGroupBy:
        pass

    class DaskSeriesGroupBy:
        pass


try:
    from zappy.base import ZappyArray
except ImportError:

    class ZappyArray:
        pass


__all__ = [
    "cache",
    "DaskArray",
    "ZappyArray",
    "DaskDataFrame",
    "DaskSeries",
    "DaskDataFrameGroupBy",
    "DaskSeriesGroupBy",
    "fullname",
    "pkg_metadata",
    "pkg_version",
]


def fullname(typ: type) -> str:
    module = typ.__module__
    name = typ.__qualname__
    if module == "builtins" or module is None:
        return name
    return f"{module}.{name}"


try:
    from contextlib import chdir
except ImportError:  # Python < 3.11
    import os
    from contextlib import AbstractContextManager

    @dataclass
    class chdir(AbstractContextManager):
        path: Path
        _old_cwd: list[Path] = field(default_factory=list)

        def __enter__(self) -> None:
            self._old_cwd.append(Path.cwd())
            os.chdir(self.path)

        def __exit__(self, *_excinfo) -> None:
            os.chdir(self._old_cwd.pop())


def pkg_metadata(package):
    from importlib.metadata import metadata as m

    return m(package)


@cache
def pkg_version(package):
    from importlib.metadata import version as v

    return version.parse(v(package))


old_positionals = partial(legacy_api, category=FutureWarning)
