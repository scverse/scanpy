"""Add `array-support` directive."""

from __future__ import annotations

from itertools import groupby
from typing import TYPE_CHECKING

from docutils import nodes
from sphinx.util.docutils import SphinxDirective

from scanpy._utils import _docs

if TYPE_CHECKING:
    from collections.abc import Collection, Generator, Iterable, Sequence
    from typing import ClassVar

    from sphinx.application import Sphinx


ALL_INNER = list(_docs.parse(["np", "sp"], inner=True))


class ArraySupport(SphinxDirective):
    """Document array support."""

    required_arguments: ClassVar = 1

    @property
    def _array_support(self) -> dict[str, tuple[list[str], list[str]]]:
        return self.config.array_support

    def run(self) -> list[nodes.Node]:  # noqa: D102
        if self.arguments[0] == "all":
            return self._render_overview()

        if not self.arguments[0] not in self._array_support:
            self.error(
                f"API not in `array_support`, add it in `docs/conf.py`: {self.arguments[0]}"
            )
        array_types = list(_docs.parse(*self._array_support[self.arguments[0]]))
        headers = ("Array type", "supported", "… experimentally in dask :class:`~dask.array.Array`")
        data: list[tuple[_docs.Inner, bool, bool]] = []
        for array_type in ALL_INNER:
            dask_array_type = _docs.DaskArray(array_type)
            data.append((
                array_type,
                array_type in array_types,
                dask_array_type in array_types,
            ))

        rows = self._render_support_data(data)
        return self._render_table(headers, rows)

    def _render_overview(self) -> list[nodes.Node]:
        headers = ["Function", *(at.rst(short=True) for at in ALL_INNER)]
        rows: list[nodes.row] = []
        for fn, (include, exclude) in self._array_support.items():
            row_header, _ = self.parse_inline(f":func:`scanpy.{fn}`")
            ats = frozenset(_docs.parse(include, exclude))
            cells: list[list[nodes.Node]] = [
                row_header,
                *(
                    [
                        nodes.Text(
                            ("✅" if at in ats else "❌") + ("⚡" if dt in ats else "")
                        )
                    ]
                    for at, dt in zip(
                        ALL_INNER, map(_docs.DaskArray, ALL_INNER), strict=True
                    )
                ),
            ]
            rows.append(
                nodes.row(
                    "",
                    *(
                        nodes.entry("", nodes.paragraph("", "", *cell))
                        for cell in cells
                    ),
                )
            )
        return self._render_table(headers, rows)

    def _render_support_data(
        self,
        data: list[tuple[_docs.Inner, bool, bool]],
    ) -> Generator[nodes.row, None, None]:
        for t, group in groupby(data, key=lambda r: type(r[0])):
            group = list(group)  # noqa: PLW2901
            if (  # if all sparse types have the same support, just one row
                t is _docs.ScipySparse
                and (support := one({s for _, s, _ in group})) is not None
                and (in_dask := one({d for _, _, d in group})) is not None
            ):
                refs: list[nodes.Node] = [
                    nodes.inline("", "scipy.sparse.{"),
                    *self.parse_inline(":class:`csr <scipy.sparse.csr_array>`")[0],
                    nodes.inline("", ","),
                    *self.parse_inline(":class:`csc <scipy.sparse.csc_matrix>`")[0],
                    nodes.inline("", "}_{"),
                    *self.parse_inline(":class:`array <scipy.sparse.csc_array>`")[0],
                    nodes.inline("", ","),
                    *self.parse_inline(":class:`matrix <scipy.sparse.csr_matrix>`")[0],
                    nodes.inline("", "}"),
                ]
                header = [nodes.literal("", "", *refs)]
                yield self._render_row(header, support=support, in_dask=in_dask)
            else:  # otherwise, show them individually
                for array_type, support, in_dask in group:
                    yield self._render_row(
                        self._render_array_type(array_type),
                        support=support,
                        in_dask=in_dask,
                    )

    def _render_row(
        self, header: Sequence[nodes.Node], *, support: bool, in_dask: bool
    ) -> nodes.row:
        cells: list[Sequence[nodes.Node]] = [
            header,
            [nodes.Text("✅" if support else "❌")],
            [nodes.Text("✅" if in_dask else "❌")],
        ]
        children = (nodes.entry("", nodes.paragraph("", "", *cell)) for cell in cells)
        return nodes.row("", *children)

    def _render_table(
        self, headers: Collection[str], rows: Iterable[nodes.row]
    ) -> list[nodes.Node]:
        colspecs = [
            nodes.colspec(stub=True),
            *(nodes.colspec() for _ in range(len(headers) - 1)),
        ]
        header_nodes = [
            nodes.entry("", nodes.paragraph("", "", *self.parse_inline(t)[0]))
            for t in headers
        ]
        thead = nodes.thead("", nodes.row("", *header_nodes))
        tbody = nodes.tbody("", *rows)
        return [
            nodes.table(
                "",
                nodes.title("", "Array type support"),
                nodes.tgroup("", *colspecs, thead, tbody, cols=len(colspecs)),
                ids=["array-support"],
            )
        ]

    def _render_array_type(self, array_type: _docs.ArrayType, /) -> list[nodes.Node]:
        nodes_, msgs = self.parse_inline(array_type.rst())
        assert not msgs, msgs
        return nodes_


def one[T](arg: Collection[T]) -> T | None:
    """Return the only item in `arg` or None if `arg` is not of length 1."""
    try:
        [item] = arg
    except ValueError:
        return None
    return item


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_directive("array-support", ArraySupport)
    app.add_config_value("array_support", {}, "env")
