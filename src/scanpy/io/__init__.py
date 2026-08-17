"""Reading and Writing."""

from __future__ import annotations

from anndata.io import (
    read_csv,
    read_excel,
    read_h5ad,
    read_hdf,
    read_mtx,
    read_text,
    read_umi_tools,
)

from ._read import read, read_10x_h5, read_10x_mtx
from ._write import write, write_cellbrowser, write_spring_project

__all__ = [
    "read",
    "read_10x_h5",
    "read_10x_mtx",
    "read_csv",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "write",
    "write_cellbrowser",
    "write_spring_project",
]
