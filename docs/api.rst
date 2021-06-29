.. module:: scanpy
.. automodule:: scanpy
   :noindex:

API
===


Import Scanpy as::

   import scanpy as sc

.. note::
   Additional functionality is available in the broader :doc:`ecosystem <../ecosystem>`, with some tools being wrapped in the :mod:`scanpy.external` module.

Preprocessing: `pp`
-------------------

.. module:: scanpy.pp
.. currentmodule:: scanpy

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Any transformation of the data matrix that is not a *tool*. Other than *tools*, preprocessing steps usually don't return an easily interpretable annotation, but perform a basic transformation on the data matrix.

Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For visual quality control, see :func:`~scanpy.pl.highest_expr_genes` and
:func:`~scanpy.pl.filter_genes_dispersion` in :mod:`scanpy.pl`.

.. autosummary::
   :toctree: generated/

   pp.calculate_qc_metrics
   pp.filter_cells
   pp.filter_genes
   pp.highly_variable_genes
   pp.log1p
   pp.pca
   pp.normalize_total
   pp.regress_out
   pp.scale
   pp.subsample
   pp.downsample_counts

Recipes
~~~~~~~

.. autosummary::
   :toctree: generated/

   pp.recipe_zheng17
   pp.recipe_weinreb17
   pp.recipe_seurat

Batch effect correction
~~~~~~~~~~~~~~~~~~~~~~~

Also see `Data integration`_. Note that a simple batch correction method is available via :func:`pp.regress_out`. Checkout :mod:`scanpy.external` for more.

.. autosummary::
   :toctree: generated/

   pp.combat

Neighbors
~~~~~~~~~

.. autosummary::
   :toctree: generated/

   pp.neighbors


Tools: `tl`
-----------

.. module:: scanpy.tl
.. currentmodule:: scanpy

Any transformation of the data matrix that is not *preprocessing*. In contrast to a *preprocessing* function, a *tool* usually adds an easily interpretable annotation to the data matrix, which can then be visualized with a corresponding plotting function.

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.pca
   tl.tsne
   tl.umap
   tl.draw_graph
   tl.diffmap

Compute densities on embeddings.

.. autosummary::
   :toctree: generated/

   tl.embedding_density

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.leiden
   tl.louvain
   tl.dendrogram
   tl.dpt
   tl.paga

Data integration
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.ingest

Marker genes
~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.rank_genes_groups
   tl.filter_rank_genes_groups
   tl.marker_gene_overlap

Gene scores, Cell cycle
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.score_genes
   tl.score_genes_cell_cycle

Simulations
~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   tl.sim


Plotting: `pl`
--------------

.. module:: scanpy.pl
.. currentmodule:: scanpy

The plotting module :mod:`scanpy.pl` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

See :tutorial:`plotting/core` for an overview of how to use these functions.

.. note::
   See the :ref:`settings` section for all important plotting configurations.

.. _pl-generic:


Generic
~~~~~~~

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
~~~~~~~

These classes allow fine tuning of visual parameters.

.. autosummary::
   :toctree: generated/classes

    pl.DotPlot
    pl.MatrixPlot
    pl.StackedViolin


Preprocessing
~~~~~~~~~~~~~

Methods for visualizing quality control and results of preprocessing functions.

.. autosummary::
   :toctree: generated/

   pl.highest_expr_genes
   pl.filter_genes_dispersion
   pl.highly_variable_genes


Tools
~~~~~

Methods that extract and visualize tool-specific annotation in an
:class:`~anndata.AnnData` object.  For any method in module ``tl``, there is
a method with the same name in ``pl``.

PCA
^^^
.. autosummary::
   :toctree: generated/

   pl.pca
   pl.pca_loadings
   pl.pca_variance_ratio
   pl.pca_overview

Embeddings
^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Visualize clusters using one of the embedding methods passing ``color='louvain'``.

.. autosummary::
   :toctree: generated/

   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.paga
   pl.paga_path
   pl.paga_compare

Marker genes
^^^^^^^^^^^^
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
^^^^^^^^^^^
.. autosummary::
   :toctree: generated/

   pl.sim


Reading
-------

.. note::
   For reading annotation use :ref:`pandas.read_… <pandas:io>`
   and add it to your :class:`anndata.AnnData` object. The following read functions are
   intended for the numeric data in the data matrix `X`.

Read common file formats using

.. autosummary::
   :toctree: generated/

   read

Read 10x formatted hdf5 files and directories containing `.mtx` files using

.. autosummary::
   :toctree: generated/

   read_10x_h5
   read_10x_mtx
   read_visium

Read other formats using functions borrowed from :mod:`anndata`

.. autosummary::
   :toctree: generated/

   read_h5ad
   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools


Get object from `AnnData`: `get`
--------------------------------

.. module:: scanpy.get
.. currentmodule:: scanpy

The module `sc.get` provides convenience functions for getting values back in
useful formats.

.. autosummary::
   :toctree: generated/

   get.obs_df
   get.var_df
   get.rank_genes_groups_df


Queries
-------

.. module:: scanpy.queries
.. currentmodule:: scanpy

This module provides useful queries for annotation and enrichment.

.. autosummary::
   :toctree: generated/

   queries.biomart_annotations
   queries.gene_coordinates
   queries.mitochondrial_genes
   queries.enrich


Metrics
-------

.. module:: scanpy.metrics
.. currentmodule:: scanpy

Collections of useful measurements for evaluating results.

.. autosummary::
   :toctree: generated/

   metrics.confusion_matrix
   metrics.gearys_c
   metrics.morans_i


Classes
-------

:class:`~anndata.AnnData` is reexported from :mod:`anndata`.

Represent data as a neighborhood structure, usually a knn graph.

.. autosummary::
   :toctree: generated/

   Neighbors


.. _settings:

Settings
--------

A convenience function for setting some default :obj:`matplotlib.rcParams` and a
high-resolution jupyter display backend useful for use in notebooks.

.. autosummary::
   :toctree: generated/

   set_figure_params

An instance of the :class:`~scanpy._settings.ScanpyConfig` is available as `scanpy.settings` and allows configuring Scanpy.

.. autosummary::
   :toctree: generated/

   _settings.ScanpyConfig

Some selected settings are discussed in the following.

Influence the global behavior of plotting functions. In non-interactive scripts,
you'd usually want to set `settings.autoshow` to ``False``.

.. no :toctree: here because they are linked under the class
.. autosummary::

   ~_settings.ScanpyConfig.autoshow
   ~_settings.ScanpyConfig.autosave

The default directories for saving figures, caching files and storing datasets.

.. autosummary::

   ~_settings.ScanpyConfig.figdir
   ~_settings.ScanpyConfig.cachedir
   ~_settings.ScanpyConfig.datasetdir

The verbosity of logging output, where verbosity levels have the following
meaning: 0='error', 1='warning', 2='info', 3='hint', 4=more details, 5=even more
details, etc.

.. autosummary::

   ~_settings.ScanpyConfig.verbosity

Print versions of packages that might influence numerical results.

.. autosummary::
   :toctree: generated/

   logging.print_header
   logging.print_versions


Datasets
--------

.. module:: scanpy.datasets
.. currentmodule:: scanpy

.. autosummary::
   :toctree: generated/

   datasets.blobs
   datasets.ebi_expression_atlas
   datasets.krumsiek11
   datasets.moignard15
   datasets.pbmc3k
   datasets.pbmc3k_processed
   datasets.pbmc68k_reduced
   datasets.paul15
   datasets.toggleswitch
   datasets.visium_sge


Deprecated functions
--------------------

.. autosummary::
   :toctree: generated/

   pp.filter_genes_dispersion
   pp.normalize_per_cell
