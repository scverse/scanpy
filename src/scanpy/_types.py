from __future__ import annotations

from typing import Literal, NamedTuple

__all__ = ["FilterCellsCutoffs", "FilterGenesCutoffs", "HVGFlavor"]

HVGFlavor = Literal["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]


class FilterCellsCutoffs(NamedTuple):
    min_genes: int | None
    min_counts: int | None
    max_genes: int | None
    max_counts: int | None

    @property
    def n(self) -> int:
        return sum([i is not None for i in self])


class FilterGenesCutoffs(NamedTuple):
    min_cells: int | None
    min_counts: int | None
    max_cells: int | None
    max_counts: int | None

    @property
    def n(self) -> int:
        return sum([i is not None for i in self])
