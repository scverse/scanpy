# some technical stuff
import sys
from ._utils import pkg_version, check_versions, annotate_doc_types

__author__ = ', '.join([
    'Alex Wolf',
    'Philipp Angerer',
    'Fidel Ramirez',
    'Isaac Virshup',
    'Sergei Rybakov',
    'Gokcen Eraslan',
    'Tom White',
    'Malte Luecken',
    'Davide Cittaro',
    'Tobias Callies',
    'Marius Lange',
    'Andrés R. Muñoz-Rojas',
])
__email__ = ', '.join([
    'f.alex.wolf@gmx.de',
    'philipp.angerer@helmholtz-muenchen.de',
    # We don’t need all, the main authors are sufficient.
])
try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    __version__ = str(pkg_version(__name__))

check_versions()
del pkg_version, check_versions

# the actual API
from ._settings import settings, Verbosity  # start with settings as several tools are using it
from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from . import datasets, logging, queries, external, get

from anndata import AnnData
from anndata import read_h5ad, read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools
from .readwrite import read, read_10x_h5, read_10x_mtx, write, read_visium
from .neighbors import Neighbors

set_figure_params = settings.set_figure_params

# has to be done at the end, after everything has been imported
annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, annotate_doc_types
