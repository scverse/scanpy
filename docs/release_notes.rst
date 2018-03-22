See all releases `here <https://github.com/theislab/scanpy/releases>`_. The following lists selected improvements.


**Soon**

- more canonical analyses steps like clustering genes, computing correlations...


**March..., 2018**: version 1.0

.. warning::

   Upgrading to 1.0 isn't fully backwards compat. First, you need to run ``conda
   install numba`` as new dependency (others have disappeared).

   - replace occurances of `group_by` with `groupby` (consistency with
     `pandas`), of `flavor` with `method`

   - cleaner nested structure of storing unstructered annotations as dicts
     

- basic graph now available as ``.uns['neighbors']['connectivities']`` after call of `pp.neighbors` or Neighbors.connectivities
      
- logistic regrssion for finding marker genes :func:`~scanpy.api.rank_genes_groups` with parameter `metfod='logreg'`
      
- scanpy no longer modifies rcParams upon import, simple call `sc.settings.set_figure_params(scanpy=True)` to set the `scanoy` style
      
- new cache directory is now ``./cache/``

- show edges in scatter plots

- neighbors class
  * distances as .distances
  * connectivities as .connectivities
  * transition matrix as .transitions

- changed connectivity measure (dpt will look different in default) set `method='gauss'` to use connectivities from gauss kernel

- downsample_counts function

- default 'louvain_groups' are now called 'louvain'

- X_diffmap now contains the zero component, plotting unchanged
  
- paul15_raw is no longer available and now is called paul15
   
- seurat preprocessing recipe

- removed cython stuff

  
3. exporting to Gephi...

Graph tools now require explicitly constructing the neighborhood graphs as a preprocessing step::

    sc.pp.neighbors(adata, n_neighbors=5, knn=False)
    sc.tl.draw_graph(adata)
    sc.tl.dpt(adata, n_branchings=1)
    sc.tl.louvain(adata, resolution=1.5)

instead of previously::

    sc.tl.draw_graph(adata, n_neighbors=5)
    sc.tl.dpt(adata, n_branchings=1, n_neighbors=5, knn=False)  # n_neighbors not necessary
    sc.tl.louvain(adata, resolution=1.5, n_neighbors=5)  # n_neighbors not necessary


**February 26, 2018**: version 0.4.4

1. embed cells using :func:`~scanpy.api.tl.umap` [McInnes18]_: `examples <https://github.com/theislab/scanpy/pull/92>`_
2. score sets of genes, e.g. for cell cycle, using :func:`~scanpy.api.tl.score_genes` [Satija15]_: `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb>`_


**February 9, 2018**: version 0.4.3

1. :func:`~scanpy.api.pl.clustermap`: heatmap from hierarchical clustering,
   based on `seaborn.clustermap
   <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_ [Waskom16]_
2. only return `matplotlib.Axis` in plotting functions of ``sc.pl`` when `show=False`, otherwise `None`

... and through `anndata v0.5 <http://anndata.readthedocs.io>`_

1. inform about duplicates in :class:`~scanpy.api.AnnData.var_names` and resolve them using :func:`~scanpy.api.AnnData.var_names_make_unique`
2. by default, generate unique observation names in :func:`~scanpy.api.AnnData.concatenate`
3. automatically remove unused categories after slicing
4. read/write `.loom` files using loompy 2


**January 7, 2018**: version 0.4.2

1. amendments in `AGA <https://github.com/theislab/graph_abstraction>`_
   and its plotting functions


**December 23, 2017**: version 0.4

1. export to `SPRING <https://github.com/AllonKleinLab/SPRING/>`_ [Weinreb17]_
   for interactive visualization of data: `tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`_,
   `docs <https://scanpy.readthedocs.io/en/latest/api/index.html>`_

... and through `anndata v0.4 <http://anndata.readthedocs.io>`_

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


**November 16, 2017**: version 0.3

1. :class:`~scanpy.api.AnnData` can be `concatenated <https://scanpy.readthedocs.io/en/latest/api/scanpy.api.AnnData.html>`_
2. :class:`~scanpy.api.AnnData` is available as a `separate package <https://pypi.python.org/pypi/anndata/>`_
3. results of approximate graph abstraction (AGA) are `simplified <https://github.com/theislab/graph_abstraction>`_


**October 25, 2017**: version 0.2.9

Initial release of `approximate graph abstraction (AGA) <https://github.com/theislab/graph_abstraction>`_.


**July 24, 2017**: version 0.2.1

Scanpy now includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The implementation efficiently deals with datasets of more than one million cells.


**May 1, 2017**: version 0.1

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing most of Seurat's guided clustering tutorial.
