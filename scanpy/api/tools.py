"""Tools

The user only sees convenience functions.
"""

# order alphabetically
from ..tools.dbscan import dbscan
from ..tools.diffmap import diffmap
from ..tools.rank_genes_groups import rank_genes_groups
from ..tools.dpt import dpt
from ..tools.louvain import louvain
from ..tools.pca import pca
from ..tools.sim import sim
from ..tools.spring import spring
from ..tools.tsne import tsne

try:
    # development tools
    from ..tools.draw_graph import draw_graph
    from ..tools.aga import aga
    from ..tools.aga import aga_contract_graph
except ImportError:
    pass
