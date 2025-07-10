from __future__ import annotations

from contextlib import contextmanager
from enum import EnumMeta, IntEnum
from logging import getLevelName
from typing import TYPE_CHECKING

from .._compat import deprecated

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal

    _VerbosityName = Literal["error", "warning", "info", "hint", "debug"]
    _LoggingLevelName = Literal["CRITICAL", "ERROR", "WARNING", "INFO", "HINT", "DEBUG"]


_VERBOSITY_TO_LOGLEVEL: dict[int | _VerbosityName, _LoggingLevelName] = {
    "error": "ERROR",
    "warning": "WARNING",
    "info": "INFO",
    "hint": "HINT",
    "debug": "DEBUG",
}
_VERBOSITY_TO_LOGLEVEL.update(dict(enumerate(list(_VERBOSITY_TO_LOGLEVEL.values()))))


class VerbosityMeta(EnumMeta):
    @property
    @deprecated("Use `Verbosity.warning` instead")
    def warn(cls) -> Verbosity:
        return Verbosity.warning


class Verbosity(IntEnum, metaclass=VerbosityMeta):
    """Logging verbosity levels for :attr:`scanpy.settings.verbosity`."""

    error = 0
    """Error (`0`)"""
    warning = 1
    """Warning (`1`)"""
    info = 2
    """Info (`2`)"""
    hint = 3
    """Hint (`3`)"""
    debug = 4
    """Debug (`4`)"""

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Verbosity):
            return self is other
        if isinstance(other, int):
            return self.value == other
        if isinstance(other, str):
            return self.name == other
        return NotImplemented

    def __hash__(self) -> int:
        # See https://docs.astral.sh/ruff/rules/eq-without-hash/
        return super().__hash__()

    @property
    def level(self) -> int:
        """The :ref:`logging level <levels>` corresponding to this verbosity level."""
        # getLevelName(str) returns the int level…
        return getLevelName(_VERBOSITY_TO_LOGLEVEL[self.name])

    @contextmanager
    def override(
        self, verbosity: Verbosity | _VerbosityName | int
    ) -> Generator[Verbosity, None, None]:
        """Temporarily override verbosity.

        >>> import scanpy as sc
        >>> sc.settings.verbosity = sc.Verbosity.info
        >>> with sc.settings.verbosity.override(sc.settings.verbosity.debug):
        ...     sc.settings.verbosity
        <Verbosity.debug: 4>
        >>> sc.settings.verbosity
        <Verbosity.info: 2>
        """
        from scanpy import settings

        settings.verbosity = verbosity
        try:
            yield self
        finally:
            settings.verbosity = self
