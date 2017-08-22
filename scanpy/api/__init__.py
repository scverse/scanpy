"""A high-level API for analyzing single-cell data.

Import this as::

   import scanpy.api as sc   

Consider using the builtin abbreviations (``sc.pp``, ``sc.tl``, ``sc.pl``) for
the main modules (``sc.preprocessing``, ``sc.tools``, ``sc.plotting``) Scanpy's
API.

.. raw:: html

   <h4>Main modules</h4>

.. autosummary::
   :toctree: generated/

   preprocessing                
   tools
   plotting
   
.. raw:: html

   <h4>Reading and Writing</h4>

.. autosummary::
   :toctree: generated/

   read
   write
   read_10x_h5
   
.. raw:: html

   <h4>Data Structures</h4>

.. autosummary::
   :toctree: generated/

   AnnData
   DataGraph
"""

from .. import __version__

from .. import settings
"""Settings"""
from .. import logging
"""Logging"""
from . import tools
tl = tools       # abbreviation
"""Tools"""
from . import plotting
pl = plotting    # abbreviation
"""Plotting functions"""
from . import preprocessing
pp = preprocessing  # abbreviation
"""Preprocessing functions"""
from ..readwrite import read, read_10x_h5, write, read_params, write_params
"""Reading and writing."""
from .. import examples
from ..examples import init_run, read_run, write_run
"""Manage runs and builtin examples."""
from ..data_structs import AnnData, DataGraph
"""Main class for storing an annotated data matrix."""
from .. import utils
"""Utils."""
