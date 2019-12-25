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


Version 1.4.5 :small:`December 25, 2019`
----------------------------------------

New functionality:

- :func:`~scanpy.tl.ingest` integrating new data with embeddings and annotations of reference data, see the `ingest tutorial`_ :pr:`651` :smaller:`S Rybakov`
- :mod:`~scanpy.queries` recieved many updates. This includes enrichment through gprofiler_ and more advanced biomart queries :pr:`467` :smaller:`I Virshup`

.. _gprofiler: https://biit.cs.ut.ee/gprofiler/
.. _ingest tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/integrating-pbmcs-using-ingest.html

Code design:

- :mod:`~scanpy.pp.downsample_counts` now always preserves the dtype of it's input, instead of converting to floats to int :noteversion:`1.4.5` :pr:`865` :smaller:`I Virshup`
- allow specifying a base for :func:`~scanpy.pp.log1p` :pr:`931` :smaller:`G Eraslan`
- run neighbors on a GPU using rapids :pr:`850` :smaller:`T White`

.. warning::

   * Changed default solver of :func:`~scanpy.tl.pca` from `auto` to `arpack`.
   * Changed default behavior of `use_raw` of :func:`~scanpy.tl.score_genes` from `False` to `None`.
