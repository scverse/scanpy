"""Add `array-support` directive."""

from __future__ import annotations

from typing import TYPE_CHECKING

from docutils import nodes
from sphinx.util.docutils import SphinxDirective

from scanpy._utils import _docs

if TYPE_CHECKING:
    from typing import ClassVar

    from sphinx.application import Sphinx


class ArraySupport(SphinxDirective):
    """In the scanpy-tutorials repo, this links to the canonical location (here!)."""

    required_arguments: ClassVar = 1
    optional_arguments: ClassVar = 999

    option_spec: ClassVar = {
        "except": lambda arg: arg.split(" "),
    }

    def run(self) -> list[nodes.Node]:  # noqa: D102
        array_types = list(_docs.parse(self.arguments, self.options.get("except", ())))
        colspecs = [nodes.colspec(stub=True), *(nodes.colspec() for _ in range(2))]
        thead = nodes.thead(
            "",
            nodes.row(
                "",
                *(
                    nodes.entry("", nodes.paragraph("", t))
                    for t in ["", "direct", "dask"]
                ),
            ),
        )
        tbody = nodes.tbody()
        for array_type in _docs.parse(["np", "sp"], inner=True):
            dask_array_type = _docs.DaskArray(array_type)
            tbody += [
                nodes.row(
                    "",
                    nodes.entry("", nodes.paragraph("", str(array_type))),
                    *(
                        nodes.entry("", nodes.paragraph("", str(at in array_types)))
                        for at in [array_type, dask_array_type]
                    ),
                )
            ]
        return [
            nodes.table(
                "",
                nodes.title("", "Array support"),
                nodes.tgroup("", *colspecs, thead, tbody, cols=3),
            )
        ]


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_directive("array-support", ArraySupport)
