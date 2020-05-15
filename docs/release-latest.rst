.. role:: small
.. role:: smaller

1.5.0 :small:`2020-05-15`
~~~~~~~~~~~~~~~~~~~~~~~~~~

The `1.5.0` release adds a lot of new functionality, much of which takes advantage updates made to :mod:`anndata`. A few highlights of this release include support for 10X's visium platform, a more efficient PCA, and a generalization of how we handle cell by cell graphs. Thanks to all the contributors who made this possible!

.. rubric:: New functionality – 10X Visium support

- :func:`~scanpy.read_visium` read function for Visium data :pr:`1034` :smaller:`G Palla, P Angerer, I Virshup`
- :func:`~scanpy.datasets.visium_sge` download and import Visium datasets from 10x genomics website :pr:`1013` :smaller:`M Mirkazemi, G Palla, P Angerer`
- :func:`~scanpy.pl.spatial` plot Visium data :pr:`1012` :smaller:`G Palla, P Angerer`
- New spatial data tutorials for basic analysis :tutorial:`spatial/basic-analysis` and integration with single cell data :tutorial:`spatial/integration-scanorama`

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

.. rubric:: Performance

.. warning::

   The new :func:`~scanpy.pp.pca` implementation can result in slightly different results than previous releases when passed a sparse matrix. See the pr (:pr:`1066`) and documentation for more info.

- :func:`~scanpy.pp.pca` now uses efficient implicit centering for sparse matrices. This can lead to signifigantly improved performance for large datasets :pr:`1066` :smaller:`A Tarashansky`
- :func:`~scanpy.tl.score_genes` now has an efficient implementation for sparse matrices with missing values. :pr:`1196` :smaller:`redst4r`.

.. rubric:: Bug fixes

- :func:`~scanpy.pp.combat` now works when `obs_names` aren't unique. :pr:`1215` :smaller:`I Virshup`
- :func:`~scanpy.pp.scale` can now be used on dense arrays without centering :pr:`1160` :smaller:`simonwm`
- :func:`~scanpy.pp.regress_out` now works when some features are constant :pr:`1194` :smaller:`simonwm`
- Fixed bug in :func:`~scanpy.pp.normalize_total`, which would error if the passed object was a view :pr:`1200` :smaller:`I Virshup`
- Fixed bug in :func:`~scanpy.pp.neighbors` which could cause the `n_pcs` argument to not work :pr:`1124` :smaller:`V Bergen`
- Fixed out of date urls in :func:`~scanpy.datasets.ebi_expression_atlas` :pr:`1102` :smaller:`I Virshup`
- Fix :func:`~scanpy.tl.ingest` for UMAP `v0.4+` :pr:`1165` :smaller:`S Rybakov`
- Fix :func:`~scanpy.tl.louvain` for louvain `v0.6+` :pr:`1197` :smaller:`I Virshup`
