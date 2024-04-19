from __future__ import annotations

from typing import TYPE_CHECKING

import scanpy as sc

if TYPE_CHECKING:
    from anndata import AnnData

_pbmc68k_reduced: AnnData | None = None
_pbmc3k: AnnData | None = None


def pbmc68k_reduced():
    global _pbmc68k_reduced
    if _pbmc68k_reduced is None:
        _pbmc68k_reduced = sc.datasets.pbmc68k_reduced()
    return _pbmc68k_reduced.copy()


def pbmc3k():
    global _pbmc3k
    if _pbmc3k is None:
        _pbmc3k = sc.datasets.pbmc3k()
    return _pbmc3k.copy()
