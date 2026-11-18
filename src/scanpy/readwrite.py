"""Deprecated private module, use :mod:`scanpy.io` instead."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from . import io
from ._compat import warn
from .io._download import download

if TYPE_CHECKING:
    from os import PathLike
    from typing import Any


def _present_or_download(
    path: PathLike[str] | str, backup_url: str | None = None
) -> bool:
    path = Path(path)
    if path.is_file():
        return True
    if backup_url is None:
        return False
    download(backup_url, path)
    return True


def __getattr__(name: str) -> Any:
    if name == "_check_datafile_present_and_download":
        msg = (
            "`scanpy.readwrite` is private; download files yourself, e.g. with `pooch`"
        )
        warn(msg, FutureWarning)
        return _present_or_download
    if name in io.__all__:
        msg = f"`scanpy.readwrite` is private, import `{name}` from `scanpy.io` instead"
        warn(msg, FutureWarning)
        return getattr(io, name)
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
