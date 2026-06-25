from __future__ import annotations

import warnings
from functools import cache
from importlib.util import find_spec
from pathlib import Path
from types import FunctionType
from typing import TYPE_CHECKING

from packaging.version import Version
from scipy import sparse

if TYPE_CHECKING:
    from collections.abc import Callable
    from importlib.metadata import PackageMetadata


__all__ = [
    "CSBase",
    "CSCBase",
    "CSRBase",
    "DaskArray",
    "SpBase",
    "fullname",
    "get_namespace",
    "is_array_api",
    "pkg_metadata",
    "pkg_version",
    "set_module",
    "warn",
]


SpBase = sparse.spmatrix | sparse.sparray  # noqa: TID251
"""Only use when you directly convert it to a known subclass."""

_CSArray = sparse.csr_array | sparse.csc_array  # noqa: TID251
"""Only use if you want to specially handle arrays as opposed to matrices."""

_CSMatrix = sparse.csr_matrix | sparse.csc_matrix  # noqa: TID251
"""Only use if you want to specially handle matrices as opposed to arrays."""

CSRBase = sparse.csr_matrix | sparse.csr_array  # noqa: TID251
CSCBase = sparse.csc_matrix | sparse.csc_array  # noqa: TID251
CSBase = _CSArray | _CSMatrix


if TYPE_CHECKING:
    # type checkers are confused and can only see …core.Array
    from dask.array.core import Array as DaskArray
elif find_spec("dask"):
    from dask.array import Array as DaskArray
else:
    DaskArray = type("Array", (), dict(__module__="dask.array"))


def fullname(typ: type) -> str:
    module = typ.__module__
    name = typ.__qualname__
    if module == "builtins" or module is None:
        return name
    return f"{module}.{name}"


def pkg_metadata(package: str) -> PackageMetadata:
    from importlib.metadata import metadata

    return metadata(package)


def is_array_api(x: object) -> bool:
    # returns true if x is array api compatible
    # exclusing the ones that are already handled by the script
    from array_api_compat import is_array_api_obj

    # excluding packages that are handled by both array-api-compat and script
    if isinstance(x, DaskArray):
        return False
    if isinstance(x, DaskArray):
        return False
    return is_array_api_obj(x)


def get_namespace(x):
    # get array-api namespace for x
    from array_api_compat import get_namespace

    return get_namespace(x)


@cache
def pkg_version(package: str) -> Version:
    from importlib.metadata import version

    return Version(version(package))


def set_module[T: FunctionType | type](module: str) -> Callable[[T], T]:
    def decorator(obj: T) -> T:
        obj.__module__ = module
        return obj

    return decorator


# File prefixes for us and decorators we use
_FILE_PREFIXES: tuple[str, ...] = (str(Path(__file__).parent),)


def warn(
    message: str,
    category: type[Warning],
    *,
    source: str | None = None,
    skip_file_prefixes: tuple[str, ...] = (),
    more_file_prefixes: tuple[str, ...] = (),
) -> None:
    """Issue a warning, skipping frames from certain file prefixes."""
    __tracebackhide__ = True

    if not skip_file_prefixes:
        skip_file_prefixes = (*_FILE_PREFIXES, *more_file_prefixes)
    elif more_file_prefixes:
        msg = "Cannot use both `skip_file_prefixes` and `more_file_prefixes`."
        raise TypeError(msg)
    warnings.warn(  # noqa: TID251
        message, category, source=source, skip_file_prefixes=skip_file_prefixes
    )
