from . import tl
from . import pl
from . import pp
from . import exporting

import sys
from .. import _utils

_utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, _utils


__doc__ = """\
External API
============

.. note::
   More tools that integrate well with scanpy and anndata can be found on the :doc:`ecosystem page <../ecosystem>`.

Import Scanpy's wrappers to external tools as::

   import scanpy.external as sce

If you'd like to include a tool here, consider making a pull request (:doc:`instructions <../dev/external-tools>`).
If the tool already uses `scanpy` or `anndata`, it may fit better in the :doc:`ecosystem page <../ecosystem>`.

Preprocessing: PP
------------------

Data integration
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.bbknn
   pp.harmony_integrate
   pp.mnn_correct
   pp.scanorama_integrate


Sample demultiplexing, Doublet detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.scrublet
   pp.scrublet_simulate_doublets
   pl.scrublet_score_distribution
   pp.hashsolo

Imputation
~~~~~~~~~~

Note that the fundamental limitations of imputation are still under `debate
<https://github.com/theislab/scanpy/issues/189>`__.

.. autosummary::
   :toctree: .

   pp.dca
   pp.magic


Tools: TL
----------

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phate
   tl.palantir
   tl.trimap
   tl.sam

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phenograph
   tl.harmony_timeseries
   tl.wishbone
   tl.palantir
   tl.palantir_results

Gene scores, Cell cycle
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.sandbag
   tl.cyclone


Plotting: PL
------------

.. autosummary::
   :toctree: .

   pl.phate
   pl.trimap
   pl.sam
   pl.wishbone_marker_trajectory

Exporting
---------

.. autosummary::
   :toctree: .

   exporting.spring_project
   exporting.cellbrowser
"""
