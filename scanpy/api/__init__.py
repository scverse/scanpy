from anndata import AnnData
from ..neighbors import Neighbors

from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

from .. import __version__

from . import tl
from . import pl
from . import pp
from ..readwrite import read, read_10x_h5, write, read_params, write_params
from . import datasets
from . import export_to
from . import logging

# unfortunately, we cannot put this here as long as we have simple global
# variables in settings... they couldn't be set in this case...
# the main drawback is that we have to import set_figure_params
# to show in the docs for that reason...
# it would be nice to make the simple data types "properties of the
# module"... putting setters and getters for all of them wouldn't be very nice
from .. import settings
# for now - or maybe as the permanently favored solution - put the single function here
from ..settings import set_figure_params

# some stuff that is not actually documented...
from .. import utils
from .. import rtools
