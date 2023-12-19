from __future__ import annotations

import re
from typing import TYPE_CHECKING

from sphinx.ext.napoleon import NumpyDocstring

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def process_return(lines):
    for line in lines:
        m = re.fullmatch(r"(?P<param>\w+)\s+:\s+(?P<type>[\w.]+)", line)
        if m:
            # Once this is in scanpydoc, we can use the fancy hover stuff
            yield f'**{m["param"]}** : :class:`~{m["type"]}`'
        else:
            yield line


def scanpy_parse_returns_section(self, section):
    lines_raw = list(process_return(self._dedent(self._consume_to_next_section())))
    if lines_raw[0] == ":":
        # Remove the “:” inserted by sphinx-autodoc-typehints
        # https://github.com/tox-dev/sphinx-autodoc-typehints/blob/a5c091f725da8374347802d54c16c3d38833d41c/src/sphinx_autodoc_typehints/patches.py#L66
        lines_raw.pop(0)
    lines = self._format_block(":returns: ", lines_raw)
    if lines and lines[-1]:
        lines.append("")
    return lines


def setup(app: Sphinx):
    NumpyDocstring._parse_returns_section = scanpy_parse_returns_section
