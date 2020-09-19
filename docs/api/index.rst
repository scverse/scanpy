.. module:: scanpy
.. automodule:: scanpy
   :noindex:

API
===


Import Scanpy as::

   import scanpy as sc

.. note::
   Wrappers to external functionality are found in :mod:`scanpy.external`.

Preprocessing: `pp`
-------------------

.. module:: scanpy.pp
.. currentmodule:: scanpy

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Any transformation of the data matrix that is not a *tool*. Other than *tools*, preprocessing steps usually don't return an easily interpretable annotation, but perform a basic transformation on the data matrix.

Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For visual quality control, see :func:`~scanpy.pl.highest_expr_genes` and
:func:`~scanpy.pl.filter_genes_dispersion` in :mod:`scanpy.plotting`.

.. autosummary::
   :toctree: .

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
   :toctree: .

   pp.recipe_zheng17
   pp.recipe_weinreb17
   pp.recipe_seurat

Batch effect correction
~~~~~~~~~~~~~~~~~~~~~~~

Also see `Data integration`_. Note that a simple batch correction method is available via :func:`pp.regress_out`. Checkout :mod:`scanpy.external` for more.

.. autosummary::
   :toctree: .

   pp.combat

Neighbors
~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.neighbors


Tools: `tl`
-----------

.. module:: scanpy.tl
.. currentmodule:: scanpy

Any transformation of the data matrix that is not *preprocessing*. In contrast to a *preprocessing* function, a *tool* usually adds an easily interpretable annotation to the data matrix, which can then be visualized with a corresponding plotting function.

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.pca
   tl.tsne
   tl.umap
   tl.draw_graph
   tl.diffmap

Compute densities on embeddings.

.. autosummary::
   :toctree: .

   tl.embedding_density

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.leiden
   tl.louvain
   tl.dendrogram
   tl.dpt
   tl.paga

Data integration
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.ingest

Marker genes
~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.rank_genes_groups
   tl.filter_rank_genes_groups
   tl.marker_gene_overlap

Gene scores, Cell cycle
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.score_genes
   tl.score_genes_cell_cycle

Simulations
~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.sim


Plotting: `pl`
--------------

.. module:: scanpy.pl
.. currentmodule:: scanpy

The plotting module :mod:`scanpy.plotting` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. autosummary::
   :toctree: .

   plotting


Reading
-------

.. note::
   For reading annotation use :ref:`pandas.read_… <pandas:io>`
   and add it to your :class:`anndata.AnnData` object. The following read functions are
   intended for the numeric data in the data matrix `X`.

Read common file formats using

.. autosummary::
   :toctree: .

   read

Read 10x formatted hdf5 files and directories containing `.mtx` files using

.. autosummary::
   :toctree: .

   read_10x_h5
   read_10x_mtx
   read_visium

Read other formats using functions borrowed from :mod:`anndata`

.. autosummary::
   :toctree: .

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
   :toctree:

   get.obs_df
   get.var_df
   get.rank_genes_groups_df


Queries
-------

.. module:: scanpy.queries
.. currentmodule:: scanpy

This module provides useful queries for annotation and enrichment.

.. autosummary::
   :toctree: .

   queries.biomart_annotations
   queries.gene_coordinates
   queries.mitochondrial_genes
   queries.enrich


Classes
-------

:class:`~anndata.AnnData` is reexported from :mod:`anndata`.

Represent data as a neighborhood structure, usually a knn graph.

.. autosummary::
   :toctree: .

   Neighbors


.. _settings:

Settings
--------

A convenience function for setting some default :obj:`matplotlib.rcParams` and a
high-resolution jupyter display backend useful for use in notebooks.

.. autosummary::
   :toctree: .

   set_figure_params

An instance of the :class:`~scanpy._settings.ScanpyConfig` is available as `scanpy.settings` and allows configuring Scanpy.

.. autosummary::
   :toctree: .

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
   :toctree: .

   logging.print_header
   logging.print_versions


Datasets
--------

.. module:: scanpy.datasets
.. currentmodule:: scanpy

.. autosummary::
   :toctree: .

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


Further modules
---------------

.. autosummary::
   :toctree: .

   plotting


Deprecated functions
--------------------

.. autosummary::
   :toctree: .

   pp.filter_genes_dispersion
   pp.normalize_per_cell
