"""Single-Cell Analysis in Python."""

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
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
from ._settings import settings, Verbosity
from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets, logging, queries, external, get, metrics, experimental

from anndata import AnnData, concat
from anndata import (
    read_h5ad,
    read_csv,
    read_excel,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
)
from .readwrite import read, read_10x_h5, read_10x_mtx, write, read_visium
from .neighbors import Neighbors

set_figure_params = settings.set_figure_params

# has to be done at the end, after everything has been imported
import sys

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})
from ._utils import annotate_doc_types

annotate_doc_types(sys.modules[__name__], "scanpy")
del sys, annotate_doc_types
