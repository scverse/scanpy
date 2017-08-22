"""Plotting with AnnData; visualizing annotation added by `tools`.

.. raw:: html

   <h4>Generic plotting for AnnData</h4>

.. autosummary::
   :toctree: generated/

   scatter
   violin
   ranking

.. raw:: html

   <h4>Plotting tool results</h4>

Methods that extract and visualize tool-specific annotation in an AnnData object.

.. raw:: html

   <h5>Visualization</h5>

.. autosummary::
   :toctree: generated/
   
   pca
   pca_loadings
   pca_scatter
   pca_variance_ratio
   draw_graph
   diffmap
   draw_graph

.. raw:: html

   <h5>Branching trajectories and pseudotime, clustering, differential expression</h5>

.. autosummary::
   :toctree: generated/
   
   dpt
   dpt_scatter
   dpt_groups_pseudotime
   dpt_timeseries
   louvain
   rank_genes_groups
   rank_genes_groups_violin

.. raw:: html

   <h5>Simulations</h5>

.. autosummary::
   :toctree: generated/
   
   sim
"""


from ..plotting.rcmod import reset_rcParams

from ..plotting import palettes

from ..plotting.ann_data import scatter, violin
from ..plotting.ann_data import ranking

from ..plotting.utils import matrix
from ..plotting.utils import timeseries, timeseries_subplot, timeseries_as_heatmap

from ..plotting.preprocessing import filter_genes_dispersion

from ..plotting.tools import pca, pca_loadings, pca_scatter, pca_variance_ratio
from ..plotting.tools import diffmap
from ..plotting.tools import draw_graph
from ..plotting.tools import tsne
from ..plotting.tools import aga, aga_attachedness, aga_graph, aga_path, aga_scatter
from ..plotting.tools import dpt, dpt_scatter, dpt_groups_pseudotime, dpt_timeseries
from ..plotting.tools import louvain
from ..plotting.tools import rank_genes_groups, rank_genes_groups_violin
from ..plotting.tools import sim

# reset matplotlib.rcParams
reset_rcParams()
