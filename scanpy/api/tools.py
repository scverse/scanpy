"""Tools

The user only sees convenience functions.
"""

# order alphabetically
from ..tools.dbscan import dbscan
from ..tools.diffmap import diffmap
from ..tools.diffrank import diffrank
from ..tools.dpt import dpt
from ..tools.pca import pca
from ..tools.sim import sim
from ..tools.spring import spring
from ..tools.tsne import tsne

try:
    # development tools
    from ..tools.aga import aga
    from ..tools.louvain import louvain
except ImportError:
    pass
