.. note::
   Also see the `release notes`__ of :mod:`anndata`.

.. __: https://anndata.readthedocs.io

.. role:: small
.. role:: smaller

.. sidebar:: Key Contributors

   `anndata graph`_ & `scanpy graph`_;
   the ☀ signifies current maintainers

   * `Isaac Virshup`_: anndata overhaul, diverse contributions ☀
   * Gökcen Eraslan: diverse contributions ☀
   * Sergei Rybakov: diverse contributions ☀
   * Fidel Ramirez: plotting ☀
   * `Tom White`_: distributed computing
   * Philipp Angerer: initial anndata conception/development, software quality ☀
   * `Alex Wolf`_: initial anndata & scanpy conception/development ☀
   * `Fabian Theis`_ & lab: enabling guidance, support and environment

.. _anndata graph: https://github.com/theislab/anndata/graphs/contributors
.. _scanpy graph: https://github.com/theislab/scanpy/graphs/contributors
.. _Isaac Virshup: https://twitter.com/ivirshup
.. _Tom White: https://twitter.com/tom_e_white
.. _Alex Wolf: https://twitter.com/falexwolf
.. _Fabian Theis: https://twitter.com/fabian_theis


Version 1.4.5
-------------

New functionality:

- :func:`~scanpy.tl.ingest` integrates embeddings and annotations of an `adata` with a reference dataset, see the `ingest tutorial`_ :pr:`651` :smaller:`thanks to S Rybakov`
- :mod:`scanpy.queries` recieved many updates. This includes enrichment through gprofiler_ and more advanced biomart queries :pr:`467` :smaller:`thanks to I Virshup`

.. _gprofiler: https://biit.cs.ut.ee/gprofiler/
.. _ingest tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/integrating-pbmcs-using-ingest.html

Code design:

- :mod:`scanpy.pp.downsample_counts` now always preserves the dtype of it's input, instead of converting to floats to int :noteversion:`1.4.5` :pr:`865` :smaller:`thanks to I Virshup`
- allow specifying a base for :func:`~scanpy.pp.log1p` :pr:`931` :smaller:`thanks to G Eraslan`
- run neighbors on a GPU using rapids :pr:`850` :smaller:`thanks to T White`

And many further improvements and bug fixes.
