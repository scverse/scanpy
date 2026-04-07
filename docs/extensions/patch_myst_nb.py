"""Extension to patch https://github.com/executablebooks/MyST-NB/pull/599."""

# TODO once MyST-NB 1.1.1/1.2.0 is out, this can be removed.

from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING

from myst_nb.core.render import MditRenderMixin, NbElementRenderer

if TYPE_CHECKING:
    from docutils import nodes
    from myst_nb.core.render import MimeData, SelfType
    from sphinx.application import Sphinx


get_orig = MditRenderMixin.get_cell_level_config
ru_orig = NbElementRenderer.render_unhandled
ru_i_orig = NbElementRenderer.render_unhandled_inline


def get_cell_level_config(
    self: SelfType,
    field: str,
    cell_metadata: dict[str, object],
    line: int | None = None,
):
    """Correct version of ``MditRenderMixin.get_cell_level_config``."""
    rv = get_orig(self, field, cell_metadata, line)
    return copy(rv)


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_config_value("myst_ignore_mime_types", [], "env")
    ignore = set(app.config.myst_ignore_mime_types)

    def render_unhandled(
        self: NbElementRenderer, data: MimeData
    ) -> list[nodes.Element]:
        if data.mime_type in ignore:
            return []
        return ru_orig(self, data)

    def render_unh_i(self: NbElementRenderer, data: MimeData) -> list[nodes.Element]:
        if data.mime_type in ignore:
            return []
        return ru_i_orig(self, data)

    MditRenderMixin.get_cell_level_config = get_cell_level_config
    NbElementRenderer.render_unhandled = render_unhandled
    NbElementRenderer.render_unhandled_inline = render_unh_i
