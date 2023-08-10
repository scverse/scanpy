from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from packaging import version


try:
    from contextlib import chdir
except ImportError:
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


def pkg_version(package):
    from importlib.metadata import version as v

    return version.parse(v(package))
