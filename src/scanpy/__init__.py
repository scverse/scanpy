"""Single-Cell Analysis in Python."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

# start with settings as several tools are using it
from ._settings import Preset, Verbosity, settings  # isort: skip

from anndata import AnnData, concat

from . import datasets, experimental, get, io, logging, metrics, queries
from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
from ._utils import annotate_doc_types
from .neighbors import Neighbors

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
    "get",
    "io",
    "logging",
    "metrics",
    "pl",
    "pp",
    "queries",
    "set_figure_params",
    "settings",
    "tl",
]


from .plotting.legacy.mpl_settings import set_figure_params

annotate_doc_types(sys.modules[__name__], "scanpy")

# has to be done at the end, after everything has been imported
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})


def __getattr__(name: str) -> Any:
    import anndata

    from ._compat import warn

    if name == "__version__":
        from importlib.metadata import version

        msg = "`__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead"
        warn(msg, FutureWarning)
        return version("scanpy")

    if name == "external":
        import scanpy.external

        return scanpy.external

    if name in {"read", "read_10x_h5", "read_10x_mtx", "write"} or (
        name in io.__all__ and name in anndata.io.__all__
    ):
        msg = "Import from `scanpy.io` instead"
        warn(msg, FutureWarning)
        return getattr(io, name)
    if name == "read_loom":  # deprecated in anndata, hence missing from its `__all__`
        msg = "`read_loom` is deprecated and will be removed, use `read_h5ad` instead"
        warn(msg, FutureWarning)
        return anndata.io.read_loom
    if name == "read_visium":
        from .io._read import read_visium

        return read_visium  # warns on call already

    raise AttributeError
