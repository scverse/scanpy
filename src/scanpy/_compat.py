from __future__ import annotations

import sys
import warnings
from functools import cache, partial
from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

import legacy_api_wrap
from packaging.version import Version
from scipy import sparse

if TYPE_CHECKING:
    from importlib.metadata import PackageMetadata


__all__ = [
    "CSBase",
    "CSCBase",
    "CSRBase",
    "DaskArray",
    "SpBase",
    "deprecated",
    "fullname",
    "old_positionals",
    "pkg_metadata",
    "pkg_version",
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
    DaskArray = type("Array", (), {})
    DaskArray.__module__ = "dask.array"


def fullname(typ: type) -> str:
    module = typ.__module__
    name = typ.__qualname__
    if module == "builtins" or module is None:
        return name
    return f"{module}.{name}"


def pkg_metadata(package: str) -> PackageMetadata:
    from importlib.metadata import metadata

    return metadata(package)


@cache
def pkg_version(package: str) -> Version:
    from importlib.metadata import version

    return Version(version(package))


# File prefixes for us and decorators we use
_FILE_PREFIXES: tuple[str, ...] = (
    str(Path(__file__).parent),
    str(Path(legacy_api_wrap.__file__).parent),
)


old_positionals = partial(
    legacy_api_wrap.legacy_api,  # noqa: TID251
    category=FutureWarning,
    skip_file_prefixes=_FILE_PREFIXES,
)


# we’re not using _FILE_PREFIXES here,
# since a wholesale deprecated function shouldn’t be used internally anyway
if TYPE_CHECKING:
    from warnings import deprecated
else:
    if sys.version_info >= (3, 13):
        from warnings import deprecated as _deprecated
    else:
        from typing_extensions import deprecated as _deprecated
    deprecated = partial(_deprecated, category=FutureWarning)


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
