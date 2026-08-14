"""Wrap prose Returns sections so scverse_misc’s pydocstring parser can handle them.

`scverse_misc.sphinx_ext` (priority 100) parses and re-emits docstrings via pydocstring,
which doesn’t understand `scanpydoc`’s (priority 500) prose Returns format
and corrupts rST roles like :class:`...` by splitting on `:`.

We bracket scverse_misc with two hooks:
- priority 50: wrap prose Returns content under a dummy type entry
- priority 200: remove the dummy entry and restore the original prose
"""

from __future__ import annotations

import re
import textwrap
from itertools import count
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator

    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata


# Matches a numpy-style section underline
_UNDERLINE_RE = re.compile(r"^-+$")
# Matches a numpy identifier entry (optionally with " : type"), but NOT prose
_NUMPY_ITEM_RE = re.compile(r"^\w+( : .*)?\s*$")

_DUMMY = "dummy"


def _find_returns_sections(lines: list[str]) -> Generator[tuple[int, int]]:
    """Yield (content_start, content_end) index pairs for each Returns section.

    content_start: first line after the `-------` underline
    content_end: one past the last content line (exclusive, trailing blanks excluded)
    """
    i = 0
    while i < len(lines) - 1:
        if lines[i] != "Returns" or not _UNDERLINE_RE.match(lines[i + 1]):
            i += 1
            continue

        content_start = i + 2
        # Find the start of the next section (non-empty line followed by underline)
        content_end = next(
            k
            for k in count(content_start)
            if k >= len(lines)
            or (
                k + 1 < len(lines)
                and lines[k].strip()
                and _UNDERLINE_RE.match(lines[k + 1])
            )
        )
        i = content_end
        # Trim trailing blank lines
        while content_end > content_start and not lines[content_end - 1].strip():
            content_end -= 1
        yield (content_start, content_end)


def _first_content(lines: list[str], start: int, end: int) -> int | None:
    for i in range(start, end):
        if lines[i].strip():
            return i
    return None


def _wrap(  # noqa: PLR0917
    app: Sphinx, objtype: str, name: str, obj: object, options: object, lines: list[str]
) -> None:
    """Wrap prose Returns content under a dummy type entry (priority 50)."""
    for start, end in reversed(list(_find_returns_sections(lines))):
        first = _first_content(lines, start, end)
        if first is None or lines[first] == _DUMMY:
            continue
        # Proper numpy def-list entry — skip
        if (
            _NUMPY_ITEM_RE.match(lines[first])
            and first + 1 < end
            and lines[first + 1].startswith("    ")
        ):
            continue
        # Prose: indent all content and prepend dummy
        content = lines[start:end]
        lines[start:end] = [
            _DUMMY,
            *textwrap.indent("\n".join(content), "    ").splitlines(),
        ]


def _unwrap(  # noqa: PLR0917
    app: Sphinx, objtype: str, name: str, obj: object, options: object, lines: list[str]
) -> None:
    """Remove the dummy type entry and restore prose (priority 200)."""
    for start, end in reversed(list(_find_returns_sections(lines))):
        first = _first_content(lines, start, end)
        if first is None or lines[first] != _DUMMY:
            continue
        # Collect everything after the dummy line and dedent
        dedented = textwrap.dedent("\n".join(lines[first + 1 : end])).splitlines()
        # Strip leading blank lines introduced by the wrap
        while dedented and not dedented[0].strip():
            dedented.pop(0)
        lines[start:end] = dedented


def setup(app: Sphinx) -> ExtensionMetadata:
    """App setup hook."""
    app.connect("autodoc-process-docstring", _wrap, priority=50)
    app.connect("autodoc-process-docstring", _unwrap, priority=200)
    return {"version": "1", "parallel_read_safe": True}
