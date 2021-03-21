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


Import Scanpy's wrappers to external tools as::

   import scanpy.external as sce

If you'd like to see your tool included here, please open a `pull request <https://github.com/theislab/scanpy>`_!

Preprocessing: PP
------------------

Data integration
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: ../generated/

   pp.bbknn
   pp.harmony_integrate
   pp.mnn_correct
   pp.scanorama_integrate


Sample demultiplexing, Doublet detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: ../generated/

   pp.scrublet
   pp.scrublet_simulate_doublets
   pl.scrublet_score_distribution
   pp.hashsolo

Imputation
~~~~~~~~~~

Note that the fundamental limitations of imputation are still under `debate
<https://github.com/theislab/scanpy/issues/189>`__.

.. autosummary::
   :toctree: ../generated/

   pp.dca
   pp.magic


Tools: TL
----------

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: ../generated/

   tl.phate
   tl.palantir
   tl.trimap
   tl.sam

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: ../generated/

   tl.phenograph
   tl.harmony_timeseries
   tl.wishbone
   tl.palantir
   tl.palantir_results

Gene scores, Cell cycle
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: ../generated/

   tl.sandbag
   tl.cyclone


Plotting: PL
------------

.. autosummary::
   :toctree: ../generated/

   pl.phate
   pl.trimap
   pl.sam
   pl.wishbone_marker_trajectory

Exporting
---------

.. autosummary::
   :toctree: ../generated/

   exporting.spring_project
   exporting.cellbrowser
"""
