"""\
External API
============


Import Scanpy's wrappers to external tools as::

   import scanpy.external as sce

Preprocessing: PP
------------------

Batch effect correction
~~~~~~~~~~~~~~~~~~~~~~~

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


Tools: TL
----------

Embeddings
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phate

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


Exporting
---------

.. autosummary::
   :toctree: .

   exporting.spring_project
   exporting.cellbrowser
"""

from . import tl
from . import pl
from . import pp

from .. import _exporting as exporting

import sys
from .. import utils
utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, utils
