from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial
from pathlib import Path

from legacy_api_wrap import legacy_api
from packaging.version import Version

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
    from importlib.metadata import metadata

    return metadata(package)


@cache
def pkg_version(package):
    from importlib.metadata import version

    return Version(version(package))


old_positionals = partial(legacy_api, category=FutureWarning)
