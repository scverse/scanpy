# Author: F. Alex Wolf (http://falexwolf.de)
"""Just a namespace to be used within examples.
"""

from .. import settings
sett = settings  # abbreviation
"""Settings"""
from .. import logging
logg = logging   # abbreviation
"""Logging"""
from .. import tools
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
from ..data_structs import AnnData
"""Main class for storing an annotated data matrix."""
