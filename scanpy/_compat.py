from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial, singledispatch
from pathlib import Path

import dask.array as da
import numpy as np
from legacy_api_wrap import legacy_api
from packaging import version
from scipy import sparse as sp

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    cache = lru_cache(maxsize=None)

try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
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


@singledispatch
def sum(X: np.ndarray | sp.spmatrix, axis=None):
    return np.sum(X, axis=axis)


@sum.register
def _(X: da.Array, axis=None):
    def sum_drop_keepdims(*args, **kwargs):
        kwargs.pop("computing_meta", None)
        if isinstance(X._meta, sp.spmatrix):  # bad! why are we getting np matrices??
            kwargs.pop("keepdims", None)
            if isinstance(kwargs["axis"], tuple):
                kwargs["axis"] = kwargs["axis"][0]
        return da.chunk.sum(*args, **kwargs)

    dtype = getattr(np.zeros(1, dtype=X.dtype).sum(), "dtype", object)

    # operates on `np.matrix` for some reason with sparse chunks in dask so need explicit casting
    def aggregate_sum(*args, **kwargs):
        return da.chunk.sum(np.array(args[0]), **kwargs)

    return da.reduction(X, sum_drop_keepdims, aggregate_sum, axis=axis, dtype=dtype)


@singledispatch
def count_nonzero(X: np.ndarray, axis=None):
    return np.count_nonzero(X, axis=axis)


@count_nonzero.register
def _(X: da.Array, axis=None):
    return sum(X > 0, axis=axis)


@count_nonzero.register
def _(X: sp.spmatrix, axis=None):
    return X.getnnz(axis=axis)
