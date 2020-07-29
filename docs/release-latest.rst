.. role:: small
.. role:: smaller

1.6.0 :small:`2020-07-31`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes several new visualization options and improvements after an
overhaul of the `dotplot`, `matrixplot` and `stacked_violin` functions (see :pr:`1210` :smaller:`F Ramirez`).
In addition, the internals for the differential expression code were overhauled (:func:`~scanpy.tl.rank_genes_groups`, :pr:`1156` :smaller:`S Rybakov`).

.. rubric:: Plotting improvements for :func:`~scanpy.pl.dotplot`, :func:`~scanpy.pl.matrixplot` and :func:`~scanpy.pl.stacked_violin`

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

.. rubric:: :func:`~scanpy.pl.dotplot` changes

- Improved the colorbar and size legend for dotplots. Now the colorbar and size have titles, which can be modified using the colorbar_title and size_title arguments. They also align at the bottom of the image and do not shrink if the dotplot image is smaller.
- Allow plotting genes in rows and categories in columns (`swap_axes`).
- Using the :class:`~scanpy.pl.DotPlot` object the `dot_edge_color` and line width can be set up, a grid added as well as several other features
- New style was added in which the dots are replaced by an empty circle and the square behind the circle is colored (like in matrixplots).

.. rubric:: :func:`~scanpy.pl.stacked_violin` changes

- violin colors can be colored based on average gene expression as in dotplots
- made the linewidth of the violin plots smaller.
- removed the tics for the y axis as they tend to overlap with each other. Using the style method they can be visualized if needed.

.. rubric:: Other visualization changes

- Added title for colorbar and positioned as in dotplot for :func:`~scanpy.pl.matrixplot`
- :func:`~scanpy.pl.heatmap` and :func:`~scanpy.pl.tracksplot` now return a dictionary of axes when `show=False` as for the other plots.
- `interpolation` can be passed as parameter for :func:`~scanpy.pl.heatmap`

.. rubric:: Additions

- Added highly variable gene selection strategy from Seurat v3 :pr:`1204` :smaller:`A Gayoso`
- Add `CellRank <https://github.com/theislab/cellrank/>`_ to scanpy ecosystem :pr:`1304` :smaller:`giovp`
- Add backup_url option to :func:`~scanpy.read_10x_h5` :pr:`1296` :smaller:`A Gayoso`
- Allow prefix for :func:`~scanpy.read_10x_mtx` :pr:`1250`  :smaller:`G Sturm`

.. rubric:: Bug fixes

- Avoid warning in :func:`~scanpy.tl.rank_genes_groups` if 't-test' is passed :pr:`1303`  :smaller:`A Wolf`
- Restrict sphinx version to < 3.1, > 3.0 :pr:`1297`  :smaller:`I Virshup`
- Clean up _ranks and fix dendrogram for scipy 1.5 :pr:`1290`  :smaller:`S Rybakov`
- Use raw to translate gene symbols if applicable :pr:`1278`  :smaller:`E Rice`
- Fix diffmap (:issue:`1262`)  :smaller:`G Eraslan`
- Fix neighbors in spring_project :issue:`1260`  :smaller:`S Rybakov`
- Fix default size of dot in spatial plots :pr:`1255` :issue:`1253`  :smaller:`giovp`
- Bumped version requirement of `scipy` to `scipy>1.4` to support `rmatmat` argument of `LinearOperator` :issue:`1246` :smaller:`I Virshup`
