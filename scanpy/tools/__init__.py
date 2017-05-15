"""Tools

The user only sees convenience functions.
"""

# order alphabetically
from .dbscan import dbscan
from .diffmap import diffmap
from .diffrank import diffrank
from .dpt import dpt
from .pca import pca
from .sim import sim
from .spring import spring
from .tsne import tsne

try:
    # development tools
    from .paths_ import paths
    from .tgdyn import tgdyn
    from .tgdyn_simple import tgdyn_simple
except ImportError:
    pass
