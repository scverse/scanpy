.. automodule:: scanpy.api

Plotting: PL
============

The following describes the plotting submodule of ``scanpy.api``.

.. note::

    See the :ref:`settings` section for some important plotting configurations.


Generic
-------

.. autosummary::
   :toctree: .

   pl.scatter
   pl.ranking

Thin wrappers for Seaborn [Waskom16]_ functions.

.. autosummary::
   :toctree: .

   pl.violin
   pl.clustermap

   
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
:class:`~scanpy.api.AnnData` object.  For any method in module ``tl``, there is
a method with the same name in ``pl``.

**Embeddings**

.. autosummary::
   :toctree: .

   pl.pca
   pl.pca_loadings
   pl.pca_scatter
   pl.pca_variance_ratio
   pl.tsne
   pl.umap
   pl.diffmap
   pl.draw_graph
   pl.phate

**Branching trajectories and pseudotime, clustering**

.. autosummary::
   :toctree: .

   pl.louvain
   pl.dpt
   pl.dpt_scatter
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
