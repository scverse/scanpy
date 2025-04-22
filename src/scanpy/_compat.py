from __future__ import annotations

import os
import sys
import warnings
from dataclasses import dataclass, field
from functools import cache, partial, wraps
from importlib.util import find_spec
from pathlib import Path
from types import UnionType
from typing import TYPE_CHECKING, Literal, ParamSpec, TypeVar, cast, overload

from packaging.version import Version
from scipy import sparse

if TYPE_CHECKING:
    from collections.abc import Callable
    from functools import _SingleDispatchCallable
    from importlib.metadata import PackageMetadata


P = ParamSpec("P")
R = TypeVar("R")


_CSMatrix = sparse.csr_matrix | sparse.csc_matrix  # noqa: TID251
"""Only use if you want to specially handle matrices as opposed to arrays"""

CSRBase = sparse.csr_matrix  # noqa: TID251
CSCBase = sparse.csc_matrix  # noqa: TID251
SpBase = sparse.spmatrix  # noqa: TID251
CSBase = _CSMatrix


if TYPE_CHECKING:
    # type checkers are confused and can only see …core.Array
    from dask.array.core import Array as DaskArray
elif find_spec("dask"):
    from dask.array import Array as DaskArray
else:
    DaskArray = type("Array", (), {})
    DaskArray.__module__ = "dask.array"


if find_spec("zappy") or TYPE_CHECKING:
    from zappy.base import ZappyArray
else:
    ZappyArray = type("ZappyArray", (), {})
    ZappyArray.__module__ = "zappy.base"


__all__ = [
    "DaskArray",
    "ZappyArray",
    "_numba_threading_layer",
    "deprecated",
    "fullname",
    "njit",
    "old_positionals",
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


if sys.version_info >= (3, 11):

    def add_note(exc: BaseException, note: str) -> None:
        exc.add_note(note)
else:

    def add_note(exc: BaseException, note: str) -> None:
        if not hasattr(exc, "__notes__"):
            exc.__notes__ = []
        exc.__notes__.append(note)


if sys.version_info >= (3, 13):
    from warnings import deprecated as _deprecated
else:
    from typing_extensions import deprecated as _deprecated


deprecated = partial(_deprecated, category=FutureWarning)


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
            if not parallel:
                msg = (
                    "Detected unsupported threading environment. "
                    f"Trying to run {f.__name__} in serial mode. "
                    "In case of problems, install `tbb`."
                )
                warnings.warn(msg, stacklevel=2)
            return fns[parallel](*args, **kwargs)

        return wrapper

    return decorator if fn is None else decorator(fn)


LayerType = Literal["default", "safe", "threadsafe", "forksafe"]
Layer = Literal["tbb", "omp", "workqueue"]


LAYERS: dict[LayerType, set[Layer]] = {
    "default": {"tbb", "omp", "workqueue"},
    "safe": {"tbb"},
    "threadsafe": {"tbb", "omp"},
    "forksafe": {"tbb", "workqueue", *(() if sys.platform == "linux" else {"omp"})},
}


def _is_in_unsafe_thread_pool() -> bool:
    import threading

    current_thread = threading.current_thread()
    # ThreadPoolExecutor threads typically have names like 'ThreadPoolExecutor-0_1'
    return (
        current_thread.name.startswith("ThreadPoolExecutor")
        and _numba_threading_layer() not in LAYERS["threadsafe"]
    )


@cache
def _numba_threading_layer() -> Layer:
    """Get numba’s threading layer.

    This function implements the algorithm as described in
    <https://numba.readthedocs.io/en/stable/user/threading-layer.html>
    """
    import importlib

    import numba

    if (available := LAYERS.get(numba.config.THREADING_LAYER)) is None:
        # given by direct name
        return numba.config.THREADING_LAYER

    # given by layer type (safe, …)
    for layer in cast("list[Layer]", numba.config.THREADING_LAYER_PRIORITY):
        if layer not in available:
            continue
        if layer != "workqueue":
            try:  # `importlib.util.find_spec` doesn’t work here
                importlib.import_module(f"numba.np.ufunc.{layer}pool")
            except ImportError:
                continue
        # the layer has been found
        return layer
    msg = (
        f"No loadable threading layer: {numba.config.THREADING_LAYER=} "
        f" ({available=}, {numba.config.THREADING_LAYER_PRIORITY=})"
    )
    raise ValueError(msg)


def _register_union(
    sdc: _SingleDispatchCallable, typ: type | UnionType
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    if sys.version_info >= (3, 11) or not isinstance(typ, UnionType):
        return sdc.register(typ)

    def decorator(f: Callable[P, R]) -> Callable[P, R]:
        for subtype in typ.__args__:
            sdc.register(subtype)(f)
        return f

    return decorator
