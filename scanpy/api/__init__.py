"""Scanpy's high-level API provides an overview of all features relevant to pratical use::

   import scanpy.api as sc   



.. raw:: html

   <h3>Preprocessing tools</h3>

Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization, preprocessing recipes.

.. raw:: html

   <h4>Basic Preprocessing</h4>

.. autosummary::
   :toctree: generated/

   pp.filter_cells   
   pp.filter_genes
   pp.filter_genes_dispersion

   pp.log1p
   pp.pca
   pp.normalize_per_cell
   pp.regress_out
   pp.scale
   pp.subsample

.. raw:: html

   <h4>Recipes</h4>

.. autosummary::
   :toctree: generated/
   
   pp.recipe_zheng17
   pp.recipe_weinreb16


.. raw:: html

   <h3>Machine Learning and Statistics tools<h3>

.. raw:: html

   <h4>Visualization</h4>

.. autosummary::
   :toctree: generated/
   
   tl.pca
   tl.tsne
   tl.diffmap
   tl.draw_graph

.. raw:: html

   <h4>Branching trajectories and pseudotime, clustering, differential expression</h4>

.. autosummary::
   :toctree: generated/
   
   tl.dpt
   tl.louvain
   tl.rank_genes_groups

.. raw:: html

   <h4>Simulations</h4>

.. autosummary::
   :toctree: generated/
   
   tl.sim


.. raw:: html

   <h3>Generic methods</h3>

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

.. raw:: html

   <h3>Plotting</h3>

.. raw:: html

   <h4>Generic plotting with AnnData</h4>

.. autosummary::
   :toctree: generated/

   pl.scatter
   pl.violin
   pl.ranking

.. raw:: html

   <h4>Plotting tool results</h4>

Methods that extract and visualize tool-specific annotation in an AnnData object.

.. raw:: html

   <h5>Visualization</h5>

.. autosummary::
   :toctree: generated/
   
   pl.pca
   pl.pca_loadings
   pl.pca_scatter
   pl.pca_variance_ratio
   pl.tsne
   pl.diffmap
   pl.draw_graph

.. raw:: html

   <h5>Branching trajectories and pseudotime, clustering, differential expression</h5>

.. autosummary::
   :toctree: generated/
   
   pl.dpt
   pl.dpt_scatter
   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.louvain
   pl.rank_genes_groups
   pl.rank_genes_groups_violin

.. raw:: html

   <h5>Simulations</h5>

.. autosummary::
   :toctree: generated/
   
   pl.sim   
"""

from .. import __version__

from .. import settings
from .. import logging
from . import tl
tools = tl
from . import pl
plotting = pl
from . import pp
preprocessing = pp
from ..readwrite import read, read_10x_h5, write, read_params, write_params
from .. import examples
from ..examples import init_run, read_run, write_run
from ..data_structs import AnnData, DataGraph
from .. import utils
