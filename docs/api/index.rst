.. automodule:: scanpy.api

API
===


Import Scanpy's high-level API as::

   import scanpy.api as sc

Preprocessing: PP
------------------

Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization, preprocessing recipes.

**Basic Preprocessing**

.. autosummary::
   :toctree: .

   pp.filter_cells
   pp.filter_genes
   pp.filter_genes_dispersion
   pp.log1p
   pp.pca
   pp.normalize_per_cell
   pp.regress_out
   pp.scale
   pp.subsample
   pp.downsample_counts
   
For visual quality control, also see :func:`~scanpy.api.pl.highest_expr_gens` and :func:`pl.filter_genes_dispersion` :doc:`plotting API <plotting>`.

**Recipes**

.. autosummary::
   :toctree: .

   pp.recipe_zheng17
   pp.recipe_weinreb17

**Batch effect correction** - beyond :func:`pp.regress_out`

.. autosummary::
   :toctree: .

   pp.mnn_correct

**Neighbors**

.. autosummary::
   :toctree: .

   pp.neighbors


Tools: TL
----------

**Embeddings**

.. autosummary::
   :toctree: .

   tl.pca
   tl.tsne
   tl.umap
   tl.draw_graph
   tl.diffmap
   tl.phate

**Clustering, branching trajectories and pseudotime based on single-cell graph**

.. autosummary::
   :toctree: .

   tl.louvain
   tl.dpt
   tl.paga

**Marker genes**

.. autosummary::
   :toctree: .

   tl.rank_genes_groups

**Gene scores, Cell cycle**

.. autosummary::
   :toctree: .

   tl.score_genes
   tl.score_genes_cell_cycle
   tl.sandbag
   tl.cyclone

**Simulations**

.. autosummary::
   :toctree: .

   tl.sim


Plotting: PL
------------

The plotting :doc:`plotting API <plotting>` largely parallels the ``tl.*`` and
``pp.*`` functions. For most tools and for some preprocessing functions, you'll
find a plotting function with the same name.

.. toctree::
   :hidden:
   :maxdepth: 1

   plotting


Reading
-------

*Note:* For reading annotation use
`pandas.read_… <http://pandas.pydata.org/pandas-docs/stable/io.html>`_ and add
it to your `AnnData` object. The following read functions are intended for
the numeric data in the data matrix `X`.

Read common file formats using

.. autosummary::
   :toctree: .

   read

Read 10x formatted hdf5 files using

.. autosummary::
   :toctree: .

   read_10x_h5

Read other formats using functions borrowed from `anndata
<http://anndata.readthedocs.io>`_

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


Queries
-------

.. autosummary::
   :toctree: .

   queries.mitochondrial_genes


Classes
-------

:class:`~scanpy.api.AnnData` is borrowed from `anndata <http://anndata.readthedocs.io>`_.

.. autosummary::
   :toctree: .

   AnnData

Represent data as a neighborhood structure, usually a knn graph.

.. autosummary::
   :toctree: .

   Neighbors


.. _settings:

Settings
--------

A convenience function for setting some default ``matplotlib.rcParams`` and a
high-resolution jupyter display backend useful for use in notebooks.

.. autosummary::
   :toctree: .

   set_figure_params

Influence the global behavior of plotting functions. In non-interactive scripts,
you'd usually want to set :class:`settings.autoshow` to ``False``.

==============================================  ===================================
:class:`settings.autoshow`                      Automatically show figures (default: ``True``).
:class:`settings.autosave`                      Automatically save figures (default: ``False``).
==============================================  ===================================

The default directories for saving figures and caching files.

==============================================  ===================================
:class:`settings.figdir`                        Directory for saving figures (default: ``'./figures/'``).
:class:`settings.cachedir`                      Directory for cache files (default: ``'./cache/'``).
==============================================  ===================================

The verbosity of logging output, where verbosity levels have the following
meaning: 0='error', 1='warning', 2='info', 3='hint', 4=more details, 5=even more
details, etc.

==============================================  ===================================
:class:`settings.verbosity`                     Verbosity level (default: 1).
==============================================  ===================================

Print versions of packages that might influence numerical results.

.. autosummary::
   :toctree: .

   logging.print_versions


Datasets
--------

.. autosummary::
   :toctree: .

   datasets.blobs
   datasets.krumsiek11
   datasets.moignard15
   datasets.paul15
   datasets.toggleswitch


Exporting
---------

.. autosummary::
   :toctree: .

   export_to.spring_project
