from __future__ import annotations

import sys
from dataclasses import dataclass, field
from functools import cache, partial
from pathlib import Path

from legacy_api_wrap import legacy_api
from packaging.version import Version

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


def pkg_metadata(package):
    from importlib.metadata import metadata

    return metadata(package)


@cache
def pkg_version(package):
    from importlib.metadata import version

    return Version(version(package))


old_positionals = partial(legacy_api, category=FutureWarning)
