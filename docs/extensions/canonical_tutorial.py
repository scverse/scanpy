from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from typing import ClassVar

    from docutils import nodes
    from sphinx.application import Sphinx


class CanonicalTutorial(SphinxDirective):
    """In the scanpy-tutorials repo, this links to the canonical location (here!)."""

    required_arguments: ClassVar = 1

    def run(self) -> list[nodes.Node]:
        return []


def setup(app: Sphinx) -> None:
    app.add_directive("canonical-tutorial", CanonicalTutorial)
