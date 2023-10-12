"""Override MyST’s cite role with one that works."""

from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

from docutils import nodes, utils

if TYPE_CHECKING:
    from typing import Any
    from collections.abc import Mapping, Sequence

    from sphinx.application import Sphinx
    from docutils.parsers.rst.states import Inliner


def cite_role(
    name: str,
    rawsource: str,
    text: str,
    lineno: int,
    inliner: Inliner,
    options: Mapping[str, Any] = MappingProxyType({}),
    content: Sequence[str] = (),
) -> tuple[list[nodes.Node], list[nodes.system_message]]:
    key = utils.unescape(text)
    node = nodes.citation_reference(f'[{key}]_', key)
    return [node], []


def setup(app: Sphinx):
    app.add_role('cite', cite_role, override=True)
