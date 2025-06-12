from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from typing import Literal


__all__ = ["FilterCellsCutoffs", "FilterGenesCutoffs", "HVGFlavor"]


class HVGFlavor(NamedTuple):
    flavor: Literal["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]


class FilterCellsCutoffs(NamedTuple):
    min_genes: int | None = None
    min_counts: int | None = None
    max_genes: int | None = None
    max_counts: int | None = None

    @property
    def n(self) -> int:
        return sum([i is not None for i in self])


class FilterGenesCutoffs(NamedTuple):
    min_cells: int | None = None
    min_counts: int | None = None
    max_cells: int | None = None
    max_counts: int | None = None

    @property
    def n(self) -> int:
        return sum([i is not None for i in self])
