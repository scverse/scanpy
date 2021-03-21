Plotting API
============

.. currentmodule:: scanpy

.. note::
   See the :ref:`settings` section for all important plotting configurations.

.. _pl-generic:

Generic
-------

.. autosummary::
   :toctree: generated/

   pl.scatter
   pl.heatmap
   pl.dotplot
   pl.tracksplot
   pl.violin
   pl.stacked_violin
   pl.matrixplot
   pl.clustermap
   pl.ranking
   pl.dendrogram


Classes
-------

These classes allow fine tuning of visual parameters.

.. autosummary::
   :toctree: generated/classes

    pl.DotPlot
    pl.MatrixPlot
    pl.StackedViolin


Preprocessing
-------------

Methods for visualizing quality control and results of preprocessing functions.

.. autosummary::
   :toctree: generated/

   pl.highest_expr_genes
   pl.filter_genes_dispersion
   pl.highly_variable_genes


Tools
-----

Methods that extract and visualize tool-specific annotation in an
:class:`~anndata.AnnData` object.  For any method in module ``tl``, there is
a method with the same name in ``pl``.

PCA
~~~
.. autosummary::
   :toctree: generated/

   pl.pca
   pl.pca_loadings
   pl.pca_variance_ratio
   pl.pca_overview

Embeddings
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   pl.tsne
   pl.umap
   pl.diffmap
   pl.draw_graph
   pl.spatial
   pl.embedding

Compute densities on embeddings.

.. autosummary::
   :toctree: generated/

   pl.embedding_density

Branching trajectories and pseudotime, clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Visualize clusters using one of the embedding methods passing ``color='louvain'``.

.. autosummary::
   :toctree: generated/

   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.paga
   pl.paga_path
   pl.paga_compare

Marker genes
~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   pl.rank_genes_groups
   pl.rank_genes_groups_violin
   pl.rank_genes_groups_stacked_violin
   pl.rank_genes_groups_heatmap
   pl.rank_genes_groups_dotplot
   pl.rank_genes_groups_matrixplot
   pl.rank_genes_groups_tracksplot

Simulations
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   pl.sim
