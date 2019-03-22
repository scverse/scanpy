from . import tl
from . import pl
from . import pp

from .. import _exporting as exporting

import sys
from .. import utils
utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, utils


__doc__ = """\
External API
============


Import Scanpy's wrappers to external tools as::

   import scanpy.external as sce

Preprocessing: PP
------------------

Batch effect correction
~~~~~~~~~~~~~~~~~~~~~~~

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


Tools: TL
----------

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phate
   tl.palantir

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phenograph

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
   tl.palantir


Exporting
---------

.. autosummary::
   :toctree: .

   exporting.spring_project
   exporting.cellbrowser
"""
