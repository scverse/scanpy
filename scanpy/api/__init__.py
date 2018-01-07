"""\
API
---

Scanpy's high-level API provides an overview of all features relevant to pratical use::

   import scanpy.api as sc



.. raw:: html

   <h3>Preprocessing tools</h3>

Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization, preprocessing recipes.

.. raw:: html

   <h5>Basic Preprocessing</h5>

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

.. raw:: html

   <h5>Recipes</h5>

.. autosummary::
   :toctree: .

   pp.recipe_zheng17
   pp.recipe_weinreb16


.. raw:: html

   <h3>Machine Learning and Statistics tools<h3>

.. raw:: html

   <h5>Visualization</h5>

.. autosummary::
   :toctree: .

   tl.pca
   tl.tsne
   tl.diffmap
   tl.draw_graph

.. raw:: html

   <h5>Branching trajectories and pseudotime, clustering, differential expression</h5>

.. autosummary::
   :toctree: .

   tl.louvain
   tl.dpt
   tl.aga
   tl.rank_genes_groups

.. raw:: html

   <h5>Simulations</h5>

.. autosummary::
   :toctree: .

   tl.sim


.. raw:: html

   <h3>Data Structures</h3>

.. autosummary::
   :toctree: .

   AnnData
   DataGraph


.. raw:: html

   <h3>Reading and Writing</h3>

Note: For reading annotation use
`pandas.read_… <http://pandas.pydata.org/pandas-docs/stable/io.html>`_ and add
it to your `AnnData` object. The following read functions are intended for
numeric data.

A powerful multi-format file-reader for common file formats. The format is
inferred from the extension of the filename.

.. autosummary::
   :toctree: .

   read

Reading a variety of other formats, borrowed from `anndata
<http://anndata.readthedocs.io>`_.

.. autosummary::
   :toctree: .

   read_h5ad
   read_10x_h5
   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools

For writing, use `AnnData.write_… <http://anndata.readthedocs.io/en/latest/api.html>`_.

.. raw:: html

   <h3>Exporting</h3>

.. autosummary::
   :toctree: .

   export_to.spring_project


.. raw:: html

   <h3>Plotting</h3>

.. autosummary::
   :toctree: .

   pl.set_rcParams_Scanpy
   pl.set_rcParams_Defaults

.. raw:: html

   <h5>Generic plotting with AnnData</h5>

.. autosummary::
   :toctree: .

   pl.scatter
   pl.violin
   pl.ranking

.. raw:: html

   <h5>Plotting tool results</h5>

Methods that extract and visualize tool-specific annotation in an AnnData object.

.. raw:: html

   <h6>Visualization</h6>

.. autosummary::
   :toctree: .

   pl.pca
   pl.pca_loadings
   pl.pca_scatter
   pl.pca_variance_ratio
   pl.tsne
   pl.diffmap
   pl.draw_graph

.. raw:: html

   <h6>Branching trajectories and pseudotime, clustering, differential expression</h6>

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

.. raw:: html

   <h6>Simulations</h6>

.. autosummary::
   :toctree: .

   pl.sim


.. raw:: html

   <h3>Datasets</h3>

Simple functions that provide annotated datasets for benchmarking. See
`here <https://scanpy.readthedocs.io/en/latest/examples.html>`_ for extensive
documented tutorials and use cases.

All of these functions return an Annotated Data object.

.. autosummary::
   :toctree: .

   datasets.blobs
   datasets.krumsiek11
   datasets.moignard15
   datasets.paul15
   datasets.paul15_raw
   datasets.toggleswitch
"""

# .. raw:: html

#    <h3>Logging</h3>

# .. autosummary::
#    :toctree: .

#    logging.print_version_and_date
#    logging.print_versions_dependencies_numerics


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
