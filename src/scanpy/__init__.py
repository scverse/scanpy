"""Single-Cell Analysis in Python."""

from __future__ import annotations

import sys

from packaging.version import Version

from ._utils import check_versions
from ._version import __version__

check_versions()
del check_versions

# the actual API
# (start with settings as several tools are using it)

from ._settings import Verbosity, settings

set_figure_params = settings._set_figure_params

import anndata

if Version(anndata.__version__) >= Version("0.11.0rc2"):
    from anndata.io import (
        read_csv,
        read_excel,
        read_h5ad,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
        read_umi_tools,
    )
else:
    from anndata import (
        read_csv,
        read_excel,
        read_h5ad,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
        read_umi_tools,
    )
from anndata import AnnData, concat

from . import datasets, experimental, external, get, logging, metrics, queries
from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
from .neighbors import Neighbors
from .readwrite import read, read_10x_h5, read_10x_mtx, read_visium, write

# has to be done at the end, after everything has been imported
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})
from ._utils import annotate_doc_types

annotate_doc_types(sys.modules[__name__], "scanpy")
del sys, annotate_doc_types

__all__ = [
    "AnnData",
    "Neighbors",
    "Verbosity",
    "__version__",
    "concat",
    "datasets",
    "experimental",
    "external",
    "get",
    "logging",
    "metrics",
    "pl",
    "pp",
    "queries",
    "read",
    "read_10x_h5",
    "read_10x_mtx",
    "read_csv",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_loom",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "read_visium",
    "set_figure_params",
    "settings",
    "tl",
    "write",
]
