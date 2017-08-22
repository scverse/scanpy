"""Machine Learning and Statistics tools.

.. raw:: html

   <h4>Visualization</h4>

.. autosummary::
   :toctree: generated/
   
   pca
   tsne
   diffmap
   draw_graph

.. raw:: html

   <h4>Branching trajectories and pseudotime, clustering, differential expression</h4>

.. autosummary::
   :toctree: generated/
   
   dpt
   louvain
   rank_genes_groups

.. raw:: html

   <h4>Simulations</h4>

.. autosummary::
   :toctree: generated/
   
   sim
"""

from ..tools.pca import pca
from ..tools.tsne import tsne
from ..tools.diffmap import diffmap
from ..tools.draw_graph import draw_graph

from ..tools.aga import aga
from ..tools.aga import aga_contract_graph
from ..tools.dbscan import dbscan
from ..tools.rank_genes_groups import rank_genes_groups
from ..tools.dpt import dpt
from ..tools.louvain import louvain
from ..tools.sim import sim
