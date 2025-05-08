"""Biomart queries."""

from __future__ import annotations

from ._queries import (
    biomart_annotations,
    enrich,  # gprofiler queries
    gene_coordinates,
    mitochondrial_genes,
)

__all__ = [
    "biomart_annotations",
    "enrich",
    "gene_coordinates",
    "mitochondrial_genes",
]
