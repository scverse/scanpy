from anndata import AnnData
from ..neighbors import Neighbors

from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

from .. import __version__

from . import tl
from . import pl
from . import pp
from ..readwrite import read, read_10x_h5, read_10x_mtx, write, read_params, write_params
from . import datasets
from . import export_to
from . import logging
from . import queries

from .. import plotting

# unfortunately, we cannot put this here as long as we have simple global
# variables in settings... they couldn't be set in this case...
# the main drawback is that we have to import set_figure_params
# to show in the docs for that reason...
# it would be nice to make the simple data types "properties of the
# module"... putting setters and getters for all of them wouldn't be very nice
from .._settings import settings
# for now - or maybe as the permanently favored solution - put the single function here
# from ..settings import set_figure_params
set_figure_params = settings.set_figure_params

# some stuff that is not actually documented...
from .. import utils

import sys
utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys


__doc__ = """\
Global API (deprecated)
=======================

.. warning::

    .. deprecated:: 1.3.7

       Use the top level module instead: `import scanpy as sc`.

For the deprecated high-level API documented on this page, use `import scanpy.api as sc`.

Preprocessing: PP
------------------

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For visual quality control, see :func:`~scanpy.api.pl.highest_expr_gens` and
:func:`~scanpy.api.pl.filter_genes_dispersion` in the :doc:`plotting API
<plotting>`.

.. autosummary::
   :toctree: .

   pp.calculate_qc_metrics
   pp.filter_cells
   pp.filter_genes
   pp.highly_variable_genes
   pp.filter_genes_dispersion
   pp.log1p
   pp.pca
   pp.normalize_per_cell
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

Note that a simple batch correction method is available via :func:`pp.regress_out`.

``pp.bbknn`` is just an alias for :func:`bbknn.bbknn`. Refer to it for the documentation.

.. autosummary::
   :toctree: .

   pp.bbknn
   pp.mnn_correct

Imputation
~~~~~~~~~~

Note that the fundamental limitations of imputation are still under `debate
<https://github.com/theislab/scanpy/issues/189>`__.

.. autosummary::
   :toctree: .

   pp.dca
   pp.magic

Neighbors
~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.neighbors


Tools: TL
----------

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.pca
   tl.tsne
   tl.umap
   tl.draw_graph
   tl.diffmap
   tl.phate

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.leiden
   tl.louvain
   tl.dpt
   tl.paga

Marker genes
~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.rank_genes_groups

Gene scores, Cell cycle
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.score_genes
   tl.score_genes_cell_cycle
   tl.sandbag
   tl.cyclone

Simulations
~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.sim


Plotting: PL
------------

The plotting :doc:`plotting API <plotting>` largely parallels the ``tl.*`` and
``pp.*`` functions. For most tools and for some preprocessing functions, you'll
find a plotting function with the same name.


Reading
-------

*Note:* For reading annotation use
:ref:`pandas.read_… <pandas:/io.rst#io-tools-text-csv-hdf5>` and add
it to your `AnnData` object. The following read functions are intended for
the numeric data in the data matrix `X`.

Read common file formats using

.. autosummary::
   :toctree: .

   read

Read 10x formatted hdf5 files and directories containing `.mtx` files using

.. autosummary::
   :toctree: .

    read_10x_h5
    read_10x_mtx

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


Queries
-------

.. autosummary::
   :toctree: .

   queries.mitochondrial_genes


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
:class:`settings.datasetdir`                    Directory for example datasets (default: ``'~/scanpy_datasets/'``).
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
   datasets.pbmc3k
   datasets.pbmc68k_reduced
   datasets.paul15
   datasets.toggleswitch

Exporting
---------

.. autosummary::
   :toctree: .

   export_to.spring_project
"""
