"""Helpers for reading and writing."""

from __future__ import annotations

from pathlib import PurePath
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

__all__ = ["download", "slugify"]


def slugify(path: str | PurePath) -> str:
    """Make a path into a filename."""
    if not isinstance(path, PurePath):
        path = PurePath(path)
    parts = list(path.parts)
    if parts[0] == "/":
        parts.pop(0)
    elif len(parts[0]) == 3 and parts[0][1:] == ":\\":
        parts[0] = parts[0][0]  # C:\ → C
    filename = "-".join(parts)
    assert "/" not in filename, filename
    assert not filename[1:].startswith(":"), filename
    return filename


def download(url: str, path: Path) -> None:
    """Download `url` to `path`, unless it is already there.

    pooch writes to a temp file first, so `path` never holds a partial download
    (#4097), and it retries a few times when the transfer fails.
    """
    import pooch

    # A one-file `Pooch`: `pooch.retrieve` would do, but has no retry knob.
    pooch.create(
        path=path.parent,
        base_url="",
        registry={path.name: None},
        urls={path.name: url},
        retry_if_failed=3,
    ).fetch(path.name, progressbar=True)
