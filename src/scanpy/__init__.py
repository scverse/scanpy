"""Single-Cell Analysis in Python."""

from __future__ import annotations

import sys

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="../..", relative_to=__file__)
    del get_version
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ModuleNotFoundError:
        raise RuntimeError(
            "scanpy is not correctly installed. Please install it, e.g. with pip."
        )

from ._utils import check_versions

check_versions()
del check_versions

# the actual API
# (start with settings as several tools are using it)

from ._settings import Verbosity, settings

set_figure_params = settings.set_figure_params

from anndata import (
    AnnData,
    concat,
    read_csv,
    read_excel,
    read_h5ad,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
)

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
    "__version__",
    "AnnData",
    "concat",
    "read_csv",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_loom",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "read",
    "read_10x_h5",
    "read_10x_mtx",
    "read_visium",
    "write",
    "datasets",
    "experimental",
    "external",
    "get",
    "logging",
    "metrics",
    "queries",
    "pl",
    "pp",
    "tl",
    "Verbosity",
    "settings",
    "Neighbors",
    "set_figure_params",
]
