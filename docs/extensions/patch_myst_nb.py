"""Extension to patch https://github.com/executablebooks/MyST-NB/pull/599."""

# TODO once MyST-NB 1.1.1/1.2.0 is out, this can be removed.

from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING

from myst_nb.core.render import MditRenderMixin

if TYPE_CHECKING:
    from sphinx.application import Sphinx


get_orig = MditRenderMixin.get_cell_level_config


def get_cell_level_config(
    self: MditRenderMixin,
    field: str,
    cell_metadata: dict[str, object],
    line: int | None = None,
):
    """Correct version of ``MditRenderMixin.get_cell_level_config``."""
    rv = get_orig(self, field, cell_metadata, line)
    return copy(rv)


def setup(app: Sphinx):
    """App setup hook."""
    MditRenderMixin.get_cell_level_config = get_cell_level_config
