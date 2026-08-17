"""Single-Cell Analysis in Python."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

# start with settings as several tools are using it
from ._settings import Preset, Verbosity, settings  # isort: skip

from anndata import AnnData, concat
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

from . import datasets, experimental, external, get, logging, metrics, queries
from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
from ._utils import annotate_doc_types
from .neighbors import Neighbors
from .readwrite import read, read_10x_h5, read_10x_mtx, read_visium, write

if TYPE_CHECKING:
    from typing import Any

__all__ = [
    "AnnData",
    "Neighbors",
    "Preset",
    "Verbosity",
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


from .plotting.legacy.mpl_settings import set_figure_params

annotate_doc_types(sys.modules[__name__], "scanpy")

# has to be done at the end, after everything has been imported
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})


def __getattr__(name: str) -> Any:
    if name == "__version__":
        from importlib.metadata import version

        from ._compat import warn

        msg = "`__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead"
        warn(msg, FutureWarning)
        return version("scanpy")

    raise AttributeError
