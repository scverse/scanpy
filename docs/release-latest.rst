.. role:: small
.. role:: smaller

1.6.0 :small:`2020-07-17`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes several new visualization options and improvements after an
overhaul of the dotplot, matrixplot and stacked_violin functions (see :pr:`1210` :smaller:`F Ramirez`).
Also, this release improves significantly the internals for the
differential expression code (:func:`~scanpy.tl.rank_genes_groups`, :pr:`1156` :smaller:`SR Koncopd`).

.. rubric:: **Plotting improvements for** :func:`~scanpy.pl.dotplot`, :func:`~scanpy.pl.matrixplot` and :func:`~scanpy.pl.stacked_violin`

- Plots are now a wrapper to classes that allow fine-tuning of the images by allowing more options. The classes can be accessed directly (eg. :class:`~scanpy.pl.DotPlot`) or using the new `return_fig` parameter.
- If the plots are called after :func:`scanpy.tl.rank_genes_groups` (eg. :func:`~scanpy.pl.rank_genes_groups_dotplot`) now is also possible to plot log fold change and p-values.
- Added `ax` parameter which allows embedding the plot in other images
- Added option to include a bar plot instead of the dendrogram containing the cell/observation totals per category.
- Return a dictionary of axes for further manipulation. This includes the main plot, legend and dendrogram to totals
- Set a title to the image.
- Legend can be removed
- `groupby` can be a list of categories. E.g. `groupby=[‘tissue’, ‘cell type’]`
- Added padding parameter to dotplot and stacked_violin to address :pr:`1270`
- Updated documentation and tutorial

.. rubric:: :func:`~scanpy.pl.dotplot` **changes**

- Improved the colorbar and size legend for dotplots. Now the colorbar and size have titles, which can be modified using the colorbar_title and size_title arguments. They also align at the bottom of the image and do not shrink if the dotplot image is smaller.
- Allow plotting genes in rows and categories in columns (`swap_axes`).
- Using the :class:`~scanpy.pl.DotPlot` object the `dot_edge_color` and line width can be set up, a grid added as well as several other features
- New style was added in which the dots are replaced by an empty circle and the square behind the circle is colored (like in matrixplots).

.. rubric:: :func:`~scanpy.pl.stacked_violin` **changes**

- violin colors can be colored based on average gene expression as in dotplots
- made the linewidth of the violin plots smaller.
- removed the tics for the y axis as they tend to overlap with each other. Using the style method they can be visualized if needed.

.. rubric:: **Other visualization changes**

- Added title for colorbar and positioned as in dotplot for :func:`~scanpy.pl.matrixplot`
- :func:`~scanpy.pl.heatmap` and :func:`~scanpy.pl.tracksplot` now return a dictionary of axes when `show=False` as for the other plots.
- `interpolation` can be passed as parameter for :func:`~scanpy.pl.heatmap`

.. rubric:: **Additions**

- Added highly variable gene selection strategy from Seurat v3 :pr:`1204` :smaller:`A Gayoso`
- Add `CellRank <https://github.com/theislab/cellrank/>`_ to scanpy ecosystem :pr:`1304` :smaller:`giovp`
- Add backup_url option to :func:`~scanpy.read_10x_h5` :pr:`1296` :smaller:`A Gayoso`
- Allow prefix for :func:`~scanpy.read_10x_mtx` :pr:`1250`  :smaller:`G Sturm`

.. rubric:: **Bug fixes**

- Avoid warning in :func:`~scanpy.tl.rank_genes_groups` if 't-test' is passed :pr:`1303`  :smaller:`A Wolf`
- Restrict sphinx version to < 3.1, > 3.0 :pr:`1297`  :smaller:`I Virshup`
- Clean up _ranks and fix dendrogram for scipy 1.5 :pr:`1290`  :smaller:`SR Koncopd`
- Use raw to translate gene symbols if applicable :pr:`1278`  :smaller:`E Rice`
- Fix diffmap (:issue:`1262`)  :smaller:`G Eraslan`
- Fix neighbors in spring_project :issue:`1260`  :smaller:`SR Koncopd`
- Fix default size of dot in spatial plots :pr:`1255` :issue:`1253`  :smaller:`giovp`
- Bumped version requirement of `scipy` to `scipy>1.4` to support `rmatmat` argument of `LinearOperator` :issue:`1246` :smaller:`I Virshup`

1.5.1 :small:`2020-05-21`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Bug fixes

- Fixed a bug in :func:`~scanpy.pp.pca`, where `random_state` did not have an effect for sparse input :pr:`1240` :smaller:`I Virshup`
- Fixed docstring in :func:`~scanpy.pp.pca` which included an unused argument :pr:`1240` :smaller:`I Virshup`

1.5.0 :small:`2020-05-15`
~~~~~~~~~~~~~~~~~~~~~~~~~

The `1.5.0` release adds a lot of new functionality, much of which takes advantage of :mod:`anndata` updates `0.7.0 - 0.7.2`. Highlights of this release include support for spatial data, dedicated handling of graphs in AnnData, sparse PCA, an interface with scvi, and others.

.. rubric:: Spatial data support

- Basic analysis :tutorial:`spatial/basic-analysis` and integration with single cell data :tutorial:`spatial/integration-scanorama` :smaller:`G Palla`
- :func:`~scanpy.read_visium` read 10x Visium data :pr:`1034` :smaller:`G Palla, P Angerer, I Virshup`
- :func:`~scanpy.datasets.visium_sge` load Visium data directly from 10x Genomics :pr:`1013` :smaller:`M Mirkazemi, G Palla, P Angerer`
- :func:`~scanpy.pl.spatial` plot spatial data :pr:`1012` :smaller:`G Palla, P Angerer`

.. rubric:: New functionality

- Many functions, like :func:`~scanpy.pp.neighbors` and :func:`~scanpy.tl.umap`, now store cell-by-cell graphs in :attr:`~anndata.AnnData.obsp` :pr:`1118` :smaller:`S Rybakov`
- :func:`~scanpy.pp.scale` and :func:`~scanpy.pp.log1p` can be used on any element in :attr:`~anndata.AnnData.layers` or :attr:`~anndata.AnnData.obsm` :pr:`1173` :smaller:`I Virshup`

.. rubric:: External tools

- :func:`~scanpy.external.pp.scvi` for preprocessing with scVI :pr:`1085` :smaller:`G Xing`
- Guide for using :ref:`Scanpy in R <conversion_to_r>` :pr:`1186` :smaller:`L Zappia`

.. rubric:: Performance

- :func:`~scanpy.pp.pca` now uses efficient implicit centering for sparse matrices. This can lead to signifigantly improved performance for large datasets :pr:`1066` :smaller:`A Tarashansky`
- :func:`~scanpy.tl.score_genes` now has an efficient implementation for sparse matrices with missing values :pr:`1196` :smaller:`redst4r`.

.. warning::

   The new :func:`~scanpy.pp.pca` implementation can result in slightly different results for sparse matrices. See the pr (:pr:`1066`) and documentation for more info.

.. rubric:: Code design

- :func:`~scanpy.pl.stacked_violin` can now be used as a subplot :pr:`1084` :smaller:`P Angerer`
- :func:`~scanpy.tl.score_genes` has improved logging :pr:`1119` :smaller:`G Eraslan`
- :func:`~scanpy.pp.scale` now saves mean and standard deviation in the :attr:`~anndata.AnnData.var` :pr:`1173` :smaller:`A Wolf`
- :func:`~scanpy.external.tl.harmony_timeseries` :pr:`1091` :smaller:`A Mousa`

.. rubric:: Bug fixes

- :func:`~scanpy.pp.combat` now works when `obs_names` aren't unique. :pr:`1215` :smaller:`I Virshup`
- :func:`~scanpy.pp.scale` can now be used on dense arrays without centering :pr:`1160` :smaller:`simonwm`
- :func:`~scanpy.pp.regress_out` now works when some features are constant :pr:`1194` :smaller:`simonwm`
- :func:`~scanpy.pp.normalize_total` errored if the passed object was a view :pr:`1200` :smaller:`I Virshup`
- :func:`~scanpy.pp.neighbors` sometimes ignored the `n_pcs` param :pr:`1124` :smaller:`V Bergen`
- :func:`~scanpy.datasets.ebi_expression_atlas` which contained some out-of-date URLs :pr:`1102` :smaller:`I Virshup`
- :func:`~scanpy.tl.ingest` for UMAP `0.4` :pr:`1165` :smaller:`S Rybakov`
- :func:`~scanpy.tl.louvain` for Louvain `0.6` :pr:`1197` :smaller:`I Virshup`
- :func:`~scanpy.pp.highly_variable_genes` which could lead to incorrect results when the `batch_key` argument was used :pr:`1180` :smaller:`G Eraslan`
- :func:`~scanpy.tl.ingest` where an inconsistent number of neighbors was used :pr:`1111` :smaller:`S Rybakov`
