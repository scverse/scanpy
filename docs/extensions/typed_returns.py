# code from https://github.com/theislab/scanpy/blob/master/docs/extensions/typed_returns.py
# with some minor adjustment
from __future__ import annotations

import re
from collections.abc import Generator, Iterable

from sphinx.application import Sphinx
from sphinx.ext.napoleon import NumpyDocstring


def _process_return(lines: Iterable[str]) -> Generator[str, None, None]:
    for line in lines:
        if m := re.fullmatch(r"(?P<param>\w+)\s+:\s+(?P<type>[\w.]+)", line):
            yield f"-{m['param']} (:class:`~{m['type']}`)"
        else:
            yield line


def _parse_returns_section(self: NumpyDocstring, section: str) -> list[str]:
    lines_raw = self._dedent(self._consume_to_next_section())
    if lines_raw[0] == ":":
        del lines_raw[0]
    lines = self._format_block(":returns: ", list(_process_return(lines_raw)))
    if lines and lines[-1]:
        lines.append("")
    return lines


def setup(app: Sphinx):
    """Set app."""
    NumpyDocstring._parse_returns_section = _parse_returns_section
