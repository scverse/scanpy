See all releases `here <https://github.com/theislab/scanpy/releases>`_. The following lists selected improvements.

**To come:** already on GitHub

1. :func:`~scanpy.api.pl.clustermap`: heatmap from hierarchical clustering,
   based on `seaborn.clustermap
   <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_ [Waskom16]_


**January 7, 2018**: version 0.4.2

1. amendments in `AGA <https://github.com/theislab/graph_abstraction>`_
   and its plotting functions
2. bug fixes
   

**December 23, 2017**: version 0.4

1. :class:`~scanpy.api.AnnData` has a `.raw` attribute that simplifies
   interacting with the "raw" data: see the `clustering tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`_
2. :class:`~scanpy.api.AnnData` 
   provides scalability beyond dataset sizes that fit into memory: see this
   `blog post
   <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_ for
   more details
3. towards a common file format for exchanging :class:`~scanpy.api.AnnData` with
   packages such as Seurat and SCDE by reading and writing `.loom
   <http://loompy.org>`_ files
4. export to `SPRING <https://github.com/AllonKleinLab/SPRING/>`_ [Weinreb17]_
   for interactive visualization of data: `tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`_,
   `docs <https://scanpy.readthedocs.io/en/latest/api/index.html>`_
5. consistency updates, bug fixes, better logging  


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
