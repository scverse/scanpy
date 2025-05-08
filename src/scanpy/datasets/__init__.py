"""Builtin Datasets."""

from __future__ import annotations

from ._datasets import (
    blobs,
    burczynski06,
    krumsiek11,
    moignard15,
    paul15,
    pbmc3k,
    pbmc3k_processed,
    pbmc68k_reduced,
    toggleswitch,
    visium_sge,
)
from ._ebi_expression_atlas import ebi_expression_atlas

__all__ = [
    "blobs",
    "burczynski06",
    "ebi_expression_atlas",
    "krumsiek11",
    "moignard15",
    "paul15",
    "pbmc3k",
    "pbmc3k_processed",
    "pbmc68k_reduced",
    "toggleswitch",
    "visium_sge",
]
