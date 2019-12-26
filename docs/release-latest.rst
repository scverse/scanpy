.. role:: small
.. role:: smaller

New functionality:

- :func:`~scanpy.tl.ingest` integrating new data with embeddings and annotations of reference data, see the `ingest tutorial`_ :pr:`651` :smaller:`S Rybakov`
- :mod:`~scanpy.queries` recieved many updates. This includes enrichment through gprofiler_ and more advanced biomart queries :pr:`467` :smaller:`I Virshup`

.. _gprofiler: https://biit.cs.ut.ee/gprofiler/
.. _ingest tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/integrating-pbmcs-using-ingest.html

Code design:

- :mod:`~scanpy.pp.downsample_counts` now always preserves the dtype of it's input, instead of converting floats to ints :pr:`865` :smaller:`I Virshup`
- allow specifying a base for :func:`~scanpy.pp.log1p` :pr:`931` :smaller:`G Eraslan`
- run neighbors on a GPU using rapids :pr:`850` :smaller:`T White`

.. warning::

   * Changed default solver of :func:`~scanpy.tl.pca` from `auto` to `arpack`.
   * Changed default behavior of `use_raw` of :func:`~scanpy.tl.score_genes` from `False` to `None`.
