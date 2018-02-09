See all releases `here <https://github.com/theislab/scanpy/releases>`_. The following lists selected improvements.


**Very Soon**

1. Canonical analyses steps like clustering genes, scoring cell cycle, computing correlations...
2. Exporting to Gephi...


**February 9, 2018**: version 0.4.3

1. :func:`~scanpy.api.pl.clustermap`: heatmap from hierarchical clustering,
   based on `seaborn.clustermap
   <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_ [Waskom16]_
2. only return `matplotlib.Axis` in plotting functions when `show=True`, otherwise `None`
3. bug fixes, consistency updates

And due to relying on `anndata <http://anndata.readthedocs.io>`_: version 0.5

1. inform about duplicates in :class:`~scanpy.api.AnnData.var_names` and resolve them using :func:`~scanpy.api.AnnData.var_names_make_unique`
2. by default, generate unique observation names in :func:`~scanpy.api.AnnData.concatenate`
3. automatically remove unused categories after slicing
4. read/write `.loom` files using loompy 2
5. some IDE-backed improvements


**January 7, 2018**: version 0.4.2

1. amendments in `AGA <https://github.com/theislab/graph_abstraction>`_
   and its plotting functions
2. bug fixes


**December 23, 2017**: version 0.4

1. export to `SPRING <https://github.com/AllonKleinLab/SPRING/>`_ [Weinreb17]_
   for interactive visualization of data: `tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`_,
   `docs <https://scanpy.readthedocs.io/en/latest/api/index.html>`_
2. consistency updates, bug fixes, better logging

And due to relying on `anndata <http://anndata.readthedocs.io>`_: version 0.4

1. towards a common file format for exchanging :class:`~scanpy.api.AnnData` with
   packages such as Seurat and SCDE by reading and writing `.loom
   <http://loompy.org>`_ files
2. :class:`~scanpy.api.AnnData`
   provides scalability beyond dataset sizes that fit into memory: see this
   `blog post
   <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_
3. :class:`~scanpy.api.AnnData` has a :class:`~scanpy.api.AnnData.raw` attribute
   that simplifies storing the data matrix when you consider it "raw": see the
   `clustering tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`_


**November 29, 2017**: version 0.3.2

1. finding marker genes via :func:`~scanpy.api.pl.rank_genes_groups_violin` improved: `example <https://github.com/theislab/scanpy/issues/51>`_
2. consistency updates, better logging, docs and bug fixes throughout


**November 16, 2017**: version 0.3

1. :class:`~scanpy.api.AnnData` can be `concatenated <https://scanpy.readthedocs.io/en/latest/api/scanpy.api.AnnData.html>`_
2. :class:`~scanpy.api.AnnData` is available as a `separate package <https://pypi.python.org/pypi/anndata/>`_
3. results of approximate graph abstraction (AGA) are `simplified <https://github.com/theislab/graph_abstraction>`_
4. consistency updates, stability improvements


**October 25, 2017**: version 0.2.9

Initial release of `approximate graph abstraction (AGA) <https://github.com/theislab/graph_abstraction>`_.


**July 24, 2017**: version 0.2.1

Scanpy now includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The implementation efficiently deals with datasets of more than one million cells.


**May 1, 2017**: version 0.1

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing most of Seurat's guided clustering tutorial.
