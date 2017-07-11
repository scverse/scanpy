# Author: F. Alex Wolf (http://falexwolf.de)
"""Scanpy - Single-Cell Analysis in Python

This is the API.
"""

from .. import __version__

from .. import settings
"""Settings"""
from .. import logging
"""Logging"""
from . import tools
tl = tools       # abbreviation
"""Tools"""
from .. import plotting
pl = plotting    # abbreviation
"""Plotting functions"""
from .. import preprocessing
pp = preprocessing  # abbreviation
"""Preprocessing functions"""
from ..readwrite import read, read_10x_h5, write, read_params, write_params
"""Reading and writing."""
from .. import examples
from ..examples import init_run, read_run, write_run
"""Manage runs and builtin examples."""
from ..data_structs import AnnData
"""Main class for storing an annotated data matrix."""
from .. import utils
"""Utils."""
