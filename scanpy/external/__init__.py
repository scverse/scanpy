from . import tl
from . import pl
from . import pp
<<<<<<< HEAD

from .. import _exporting as exporting

import sys
from .. import utils
utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, utils
=======
from . import exporting

import sys
from .. import _utils

_utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys, _utils
>>>>>>> upstream/master


__doc__ = """\
External API
============


Import Scanpy's wrappers to external tools as::

   import scanpy.external as sce

<<<<<<< HEAD
=======
If you'd like to see your tool included here, please open a `pull request <https://github.com/theislab/scanpy>`_!

>>>>>>> upstream/master
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
<<<<<<< HEAD
=======
   tl.trimap
>>>>>>> upstream/master

Clustering and trajectory inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.phenograph
<<<<<<< HEAD
   tl.harmony_timeseries
=======
>>>>>>> upstream/master

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
<<<<<<< HEAD
   tl.palantir
   tl.harmony_timeseries
=======
   pl.trimap
   tl.palantir
>>>>>>> upstream/master


Exporting
---------

.. autosummary::
   :toctree: .

   exporting.spring_project
   exporting.cellbrowser
"""
