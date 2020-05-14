.. role:: small
.. role:: smaller

On master
~~~~~~~~~~

.. rubric:: Performance

.. warning::

   The new :func:`~scanpy.pp.pca` implementation can result in slightly different results than previous releases when passed a sparse matrix. See the pr (:pr:`1066`) and documentation for more info.

- :func:`~scanpy.pp.pca` now uses efficient implicit centering for sparse matrices. This can lead to signifigantly improved performance for large datasets :pr:`1066` :smaller:`A Tarashansky`
- :func:`~scanpy.tl.score_genes` now has an efficient implementation for sparse matrices with missing values. :pr:`1196` :smaller:`redst4r`.


.. rubric:: New functionality

- Many functions, like :func:`~scanpy.pp.neighbors` and :func:`~scanpy.tl.umap`, have been updated to store cell by cell graphs in :attr:`~anndata.AnnData.obsp` :pr:`1118` :smaller:`S Rybakov`
- :func:`~scanpy.pp.scale` now saves mean and standard deviation in the :attr:`~anndata.AnnData.var` :pr:`1173` :smaller:`A Wolf`
- :func:`~scanpy.pp.scale` and :func:`~scanpy.pp.log1p` can now be used on any element in :attr:`~anndata.AnnData.layers` or :attr:`~anndata.AnnData.obsm` :pr:`1173` :smaller:`I Virshup`
- :func:`~scanpy.tl.score_genes` has improved logging :pr:`1119` :smaller:`G Eraslan`
- :func:`~scanpy.pl.stacked_violin` can now be used as a subplot :pr:`1084` :smaller:`P Angerer`

.. rubric:: External tools

- Added :func:`~scanpy.external.pp.scvi` for fitting scVI model :pr:`1085` :smaller:`G Xing`
- Added a guide for using :ref:`Scanpy in R <conversion_to_r>` :pr:`1186` :smaller:`L Zappia`
- Updates to :func:`~scanpy.external.tl.harmony_timeseries` :pr:`#1091` :smaller:`A Mousa`

.. rubric:: Bug fixes

- :func:`~scanpy.pp.combat` now works when `obs_names` aren't unique. :pr:`1215` :smaller:`I Virshup`
- :func:`~scanpy.pp.scale` can now be used on dense arrays without centering :pr:`1160` :smaller:`simonwm`
- :func:`~scanpy.pp.regress_out` now works when some features are constant :pr:`1194` :smaller:`simonwm`
- Fixed bug in :func:`~scanpy.pp.normalize_total`, which would error if the passed object was a view :pr:`1200` :smaller:`I Virshup`
- Fixed bug in :func:`~scanpy.pp.neighbors` which could cause the `n_pcs` argument to not work :pr:`1124` :smaller:`V Bergen`
- Fixed out of date urls in :func:`~scanpy.datasets.ebi_expression_atlas` :pr:`1102` :smaller:`I Virshup`
- Fix :func:`~scanpy.tl.ingest` for UMAP `v0.4+` :pr:`1165` :smaller:`S Rybakov`
- Fix :func:`~scanpy.tl.louvain` for louvain `v0.6+` :pr:`1197` :smaller:`I Virshup`


1.4.6 :small:`2020-03-17`
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: Functionality in `external`

- :func:`~scanpy.external.tl.sam` self-assembling manifolds [Tarashansky19]_ :pr:`903` :smaller:`A Tarashansky`
- :func:`~scanpy.external.tl.harmony_timeseries` for trajectory inference on discrete time points :pr:`994` :smaller:`A Mousa`
- :func:`~scanpy.external.tl.wishbone` for trajectory inference (bifurcations) :pr:`1063` :smaller:`A Mousa`

.. rubric:: Code design

- :mod:`~scanpy.pl.violin` now reads `.uns['colors_...']` :pr:`1029` :smaller:`michalk8`

.. rubric:: Bug fixes

- adapt :func:`~scanpy.tl.ingest` for UMAP 0.4 :pr:`1038` :pr:`1106` :smaller:`S Rybakov`
- compat with matplotlib 3.1 and 3.2 :pr:`1090` :smaller:`I Virshup, P Angerer`
- fix PAGA for new igraph :pr:`1037` :smaller:`P Angerer`
- fix rapids compat of louvain :pr:`1079` :smaller:`LouisFaure`

1.4.5 :small:`2019-12-30`
~~~~~~~~~~~~~~~~~~~~~~~~~

Please install `scanpy==1.4.5.post3` instead of `scanpy==1.4.5`.

.. rubric:: New functionality

- :func:`~scanpy.tl.ingest` maps labels and embeddings of reference data to new data :tutorial:`integrating-data-using-ingest` :pr:`651` :smaller:`S Rybakov, A Wolf`
- :mod:`~scanpy.queries` recieved many updates including enrichment through gprofiler_ and more advanced biomart queries :pr:`467` :smaller:`I Virshup`
- :func:`~scanpy.set_figure_params` allows setting `figsize` and accepts `facecolor='white'`, useful for working in dark mode  :smaller:`A Wolf`

.. _gprofiler: https://biit.cs.ut.ee/gprofiler/

.. rubric:: Code design

- :mod:`~scanpy.pp.downsample_counts` now always preserves the dtype of it's input, instead of converting floats to ints :pr:`865` :smaller:`I Virshup`
- allow specifying a base for :func:`~scanpy.pp.log1p` :pr:`931` :smaller:`G Eraslan`
- run neighbors on a GPU using rapids :pr:`850` :smaller:`T White`
- param docs from typed params :smaller:`P Angerer`
- :func:`~scanpy.tl.embedding_density` now only takes one positional argument; similar for :func:`~scanpy.pl.embedding_density`, which gains a param `groupby` :pr:`965` :smaller:`A Wolf`
- webpage overhaul, ecosystem page, release notes, tutorials overhaul :pr:`960` :pr:`966` :smaller:`A Wolf`

.. warning::

   * changed default `solver` in :func:`~scanpy.tl.pca` from `auto` to `arpack`
   * changed default `use_raw` in :func:`~scanpy.tl.score_genes` from `False` to `None`
