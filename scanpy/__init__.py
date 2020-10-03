"""Single-Cell Analysis in Python."""

from ._metadata import __version__, __author__, __email__

from ._utils import check_versions
check_versions()
del check_versions

# the actual API
from ._settings import settings, Verbosity  # start with settings as several tools are using it
from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets, logging, queries, external, get

from anndata import AnnData
try:
    from anndata import concat
except ImportError:
    from anndata._core.merge import concat
from anndata import read_h5ad, read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools
from .readwrite import read, read_10x_h5, read_10x_mtx, write, read_visium
from .neighbors import Neighbors

set_figure_params = settings.set_figure_params

# has to be done at the end, after everything has been imported
import sys
from ._utils import annotate_doc_types
annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, annotate_doc_types
