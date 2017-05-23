# Author: F. Alex Wolf (http://falexwolf.de)
"""Scanpy - Single-Cell Analysis in Python

This is the API.
"""

from . import settings
sett = settings  # abbreviation
"""Settings"""
from . import tools
tl = tools       # abbreviation
"""Tools"""
from . import plotting
pl = plotting    # abbreviation
"""Plotting functions"""
from . import preprocessing
pp = preprocessing  # abbreviation
"""Preprocessing functions"""
from .readwrite import read, read_10x_h5, write, read_params, write_params
"""Reading and writing."""
from .examples import show_exdata, show_exparams, get_example
"""Builtin examples."""
from .data_structs import AnnData
"""Main class for storing an annotated data matrix."""
