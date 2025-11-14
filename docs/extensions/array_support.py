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
        headers = ("Array type", "supported", "… in dask")
        data: list[tuple[_docs.ArrayType, bool, bool]] = []
        for array_type in _docs.parse(["np", "sp"], inner=True):
            dask_array_type = _docs.DaskArray(array_type)
            data.append((
                array_type,
                array_type in array_types,
                dask_array_type in array_types,
            ))

        return self._render_table(headers, data)

    def _render_table(
        self,
        headers: tuple[str, str, str],
        data: list[tuple[_docs.ArrayType, bool, bool]],
    ) -> list[nodes.Node]:
        colspecs = [nodes.colspec(stub=True), *(nodes.colspec() for _ in range(2))]
        thead = nodes.thead(
            "",
            nodes.row(
                "",
                *(nodes.entry("", nodes.paragraph("", t)) for t in headers),
            ),
        )
        tbody = nodes.tbody()
        for array_type, support, in_dask in data:
            cells = [
                self._render_array_type(array_type),
                self._render_support(support),
                self._render_support(in_dask),
            ]
            children = (
                nodes.entry("", nodes.paragraph("", "", *cell)) for cell in cells
            )
            tbody += [nodes.row("", *children)]
        return [
            nodes.table(
                "",
                nodes.title("", "Array support"),
                nodes.tgroup("", *colspecs, thead, tbody, cols=3),
            )
        ]

    def _render_array_type(self, array_type: _docs.ArrayType, /) -> list[nodes.Node]:
        nodes_, msgs = self.parse_inline(array_type.rst())
        assert not msgs, msgs
        return nodes_

    def _render_support(self, support: bool, /) -> list[nodes.Node]:  # noqa: FBT001
        return [nodes.Text("✅" if support else "❌")]


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_directive("array-support", ArraySupport)
