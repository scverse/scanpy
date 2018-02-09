"""\
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

**Recipes**

.. autosummary::
   :toctree: .

   pp.recipe_zheng17
   pp.recipe_weinreb16


Tools: TL
----------

**Embeddings**

.. autosummary::
   :toctree: .

   tl.pca
   tl.tsne
   tl.diffmap
   tl.draw_graph

**Branching trajectories and pseudotime, clustering, differential expression**

.. autosummary::
   :toctree: .

   tl.aga
   tl.louvain
   tl.dpt
   tl.rank_genes_groups

**Simulations**

.. autosummary::
   :toctree: .

   tl.sim


Reading and Writing
-------------------

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

For writing, use `AnnData.write_… <http://anndata.readthedocs.io/en/latest/api.html>`_.


Exporting
---------

.. autosummary::
   :toctree: .

   export_to.spring_project


Datasets
--------

.. autosummary::
   :toctree: .

   datasets.blobs
   datasets.krumsiek11
   datasets.moignard15
   datasets.paul15
   datasets.paul15_raw
   datasets.toggleswitch


Plotting: PL
-------------

.. autosummary::
   :toctree: .

   pl.set_rcParams_Scanpy
   pl.set_rcParams_Defaults

.. raw:: html

**Generic plotting with AnnData**

.. autosummary::
   :toctree: .

   pl.scatter
   pl.ranking

Thin wrappers for Seaborn functions.

.. autosummary::
   :toctree: .

   pl.violin
   pl.clustermap

**Plotting tool results**

Methods that extract and visualize tool-specific annotation in an AnnData object.
For any method in module `tl`, there is a method with the same name in `pl`.

*Embeddings*

.. autosummary::
   :toctree: .

   pl.pca
   pl.pca_loadings
   pl.pca_scatter
   pl.pca_variance_ratio
   pl.tsne
   pl.diffmap
   pl.draw_graph

*Branching trajectories and pseudotime, clustering, differential expression*

.. autosummary::
   :toctree: .

   pl.aga
   pl.aga_graph
   pl.aga_path
   pl.louvain
   pl.dpt
   pl.dpt_scatter
   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.rank_genes_groups
   pl.rank_genes_groups_violin

*Simulations*

.. autosummary::
   :toctree: .

   pl.sim


Settings
--------

Global settings.

==============================================  ===================================
`settings.verbosity`                            Verbosity level (default: 1).
`settings.file_format_figs`                     Format for saving figures (default: 'png').
`settings.figdir`                               Default directory for saving figures (default: './figures').
`settings.autoshow`                             Automatically show figures (default: `True`).
`settings.autoshow`                             Automatically show figures (default: `True`).
==============================================  ===================================

The verbosity levels have the following meaning:

===  ======================================
 0   Only show 'error' messages.
 1   Also show 'warning' messages.
 2   Also show 'info' messages.
 3   Also show 'hint' messages.
 4   Show very detailed progress.
...  Show even more detailed progress.
===  ======================================


Logging
-------

================================================  ===================================
`logging.print_version_and_date()`                Print the version and the date.
`logging.print_versions_dependencies_numerics()`  Print the versions of dependencies.
================================================  ===================================

AnnData
-------

This only temporarily replicates the docs of `anndata
<http://anndata.readthedocs.io>`_. Occurances of :class:`~scanpy.api.AnnData`
will directly link to `anndata
<http://anndata.readthedocs.io>`_ in the future.

.. autosummary::
   :toctree: .

   AnnData
"""

from anndata import AnnData
from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

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
from . import datasets
from ..data_structs import DataGraph
from .. import utils
from . import export_to
