.. automodule:: scanpy.api

Plotting: PL
============

The following describes the plotting submodule of ``scanpy.api``.

.. note::

    See the :ref:`settings` section for all important plotting configurations.


Generic
-------

.. autosummary::
   :toctree: .

   pl.scatter
   pl.violin
   pl.heatmap
   pl.clustermap
   pl.ranking


Preprocessing
-------------

Methods for visualizing quality control and results of preprocessing functions.

.. autosummary::
   :toctree: .

   pl.highest_expr_genes
   pl.filter_genes_dispersion


Tools
-----

Methods that extract and visualize tool-specific annotation in an
:class:`~anndata.AnnData` object.  For any method in module ``tl``, there is
a method with the same name in ``pl``.

**PCA**

.. autosummary::
   :toctree: .

   pl.pca
   pl.pca_loadings
   pl.pca_scatter
   pl.pca_variance_ratio
   
**Embeddings**

.. autosummary::
   :toctree: .
             
   pl.tsne
   pl.umap
   pl.diffmap
   pl.draw_graph
   pl.phate

**Branching trajectories and pseudotime, clustering**

Visualize clusters using one of the embedding methods passing ``color='louvain'``.

.. autosummary::
   :toctree: .

   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.paga
   pl.paga_path
   pl.paga_compare

**Marker genes**

.. autosummary::
   :toctree: .

   pl.rank_genes_groups
   pl.rank_genes_groups_violin

**Simulations**

.. autosummary::
   :toctree: .

   pl.sim
