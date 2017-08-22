"""Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization, preprocessing recipes.

.. raw:: html

   <h4>Recipes</h4>

.. autosummary::
   :toctree: generated/
   
   recipe_zheng17
   recipe_weinreb16

.. raw:: html

   <h4>Basic Preprocessing</h4>

.. autosummary::
   :toctree: generated/

   filter_cells   
   filter_genes
   filter_genes_dispersion

   log1p
   pca
   normalize_per_cell
   regress_out
   scale
   subsample
"""

from ..preprocessing.simple import *
from ..preprocessing.recipes import *
