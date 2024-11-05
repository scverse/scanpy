from __future__ import annotations

import sys
from dataclasses import dataclass, field
from functools import cache, partial, wraps
from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING, ParamSpec, TypeVar, overload

from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Callable
    from importlib.metadata import PackageMetadata
    from typing import Literal

P = ParamSpec("P")
R = TypeVar("R")


if TYPE_CHECKING:
    # type checkers are confused and can only see …core.Array
    from dask.array.core import Array as DaskArray
elif find_spec("dask"):
    from dask.array import Array as DaskArray
else:

    class DaskArray:
        pass


if find_spec("zappy") or TYPE_CHECKING:
    from zappy.base import ZappyArray
else:

    class ZappyArray:
        pass


__all__ = [
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


if sys.version_info >= (3, 11):
    from contextlib import chdir
else:
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


def pkg_metadata(package: str) -> PackageMetadata:
    from importlib.metadata import metadata

    return metadata(package)


@cache
def pkg_version(package: str) -> Version:
    from importlib.metadata import version

    return Version(version(package))


if find_spec("legacy_api_wrap") or TYPE_CHECKING:
    from legacy_api_wrap import legacy_api  # noqa: TID251

    old_positionals = partial(legacy_api, category=FutureWarning)
else:
    # legacy_api_wrap is currently a hard dependency,
    # but this code makes it possible to run scanpy without it.
    def old_positionals(*old_positionals: str):
        return lambda func: func


@overload
def njit(fn: Callable[P, R], /) -> Callable[P, R]: ...
@overload
def njit() -> Callable[[Callable[P, R]], Callable[P, R]]: ...
def njit(
    fn: Callable[P, R] | None = None, /
) -> Callable[P, R] | Callable[[Callable[P, R]], Callable[P, R]]:
    """Jit-compile a function using numba.

    On call, this function dispatches to a parallel or sequential numba function,
    depending on if it has been called from a thread pool.

    See <https://github.com/numbagg/numbagg/pull/201/files#r1409374809>
    """

    def decorator(f: Callable[P, R], /) -> Callable[P, R]:
        import numba

        fns: dict[bool, Callable[P, R]] = {
            parallel: numba.njit(f, cache=True, parallel=parallel)  # noqa: TID251
            for parallel in (True, False)
        }

        @wraps(f)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            parallel = not _is_in_unsafe_thread_pool()
            return fns[parallel](*args, **kwargs)

        return wrapper

    return decorator if fn is None else decorator(fn)


def _is_in_unsafe_thread_pool() -> bool:
    import threading

    current_thread = threading.current_thread()
    # ThreadPoolExecutor threads typically have names like 'ThreadPoolExecutor-0_1'
    return (
        current_thread.name.startswith("ThreadPoolExecutor")
        and _thread_backend() == "workqueue"
    )


@cache
def _thread_backend() -> Literal["tbb", "omp", "workqueue"]:
    import importlib

    import numba

    if numba.config.THREADING_LAYER != "default":
        return numba.config.THREADING_LAYER

    for layer in numba.config.THREADING_LAYER_PRIORITY:
        if layer not in {"tbb", "omp"}:
            continue
        # Note that `importlib.util.find_spec` doesn't work for these
        # it will falsely return True
        try:
            importlib.import_module(f"numba.np.ufunc.{layer}pool")
        except ImportError:
            pass
        else:
            return layer

    return "workqueue"
