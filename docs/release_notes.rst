See the documentation of version 0.4.4 `here <http://scanpy.readthedocs.io/en/0.4.4/>`_. See a list of all releases `here <https://github.com/theislab/scanpy/releases>`_.


**Soon**

- more canonical analyses steps like clustering genes, computing correlations...

- exporting to Gephi from :class:`~scanpy.api.Neighbors`
  

**March 28, 2018**: version 1.0

Scanpy is much faster and more memory efficient. Preprocess, cluster and visualize
1.3M cells in `6 h
<https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/>`_,
130K cells in `14 min
<https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/logfile_130K.txt>`_
and 68K cells in `3 min
<https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170503_zheng17/zheng17.ipynb>`_.

The API gained a preprocessing function :func:`~scanpy.api.pp.neighbors` and a
class :func:`~scanpy.api.Neighbors` to which all basic graph computations are
delegated.

.. warning::

   Upgrading to 1.0 isn't fully backwards compatible - future upgrades will be.

   - the graph-based tools :func:`~scanpy.api.tl.louvain`
     :func:`~scanpy.api.tl.dpt` :func:`~scanpy.api.tl.draw_graph`
     :func:`~scanpy.api.tl.umap` :func:`~scanpy.api.tl.diffmap`
     :func:`~scanpy.api.tl.paga` now require prior computation of the graph:
     
     .. code:: python
     
         sc.pp.neighbors(adata, n_neighbors=5)
         sc.tl.louvain(adata)
     
     instead of previously:
     
     .. code:: python
     
         sc.tl.louvain(adata, n_neighbors=5)
         
   - install `numba` via ``conda install numba``, which replaces cython
      
   - the default connectivity measure (dpt will look different using default
     settings) changed. setting `method='gauss'` in `sc.pp.neighbors` uses
     gauss kernel connectivities and reproduces the previous behavior,
     see, for instance this `example
     <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`_

   - namings of returned annotation have changed for less bloated AnnData
     objects, which means that some of the unstructured annotation of old
     AnnData files is not recognized anymore

   - replace occurances of `group_by` with `groupby` (consistency with
     `pandas`)

   - it is worth checking out the notebook examples to see changes, e.g., `here
     <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb>`_

   - upgrading scikit-learn from 0.18 to 0.19 changed the implementation of PCA,
     some results might therefore look slightly different

Further changes are
   
- UMAP [McInnes18]_ can serve as a first visualization of the data just as tSNE,
  in contrast to tSNE, UMAP directly embeds the single-cell graph and is faster;
  UMAP is now also used for measuring connectivities and computing neighbors,
  see :func:`~scanpy.api.pp.neighbors`

- graph abstraction: AGA is renamed to PAGA: :func:`~scanpy.api.tl.paga`; now,
  it only measures connectivities between partitions of the single-cell graph,
  pseudotime and clustering need to be computed separately via
  :func:`~scanpy.api.tl.louvain` and :func:`~scanpy.api.tl.dpt`, the
  connectivity measure has been improved
           
- logistic regression for finding marker genes
  :func:`~scanpy.api.tl.rank_genes_groups` with parameter `method='logreg'`

- :func:`~scanpy.api.tl.louvain` now provides a better implementation for
  reclustering via `restrict_to`        
        
- scanpy no longer modifies rcParams upon import, call
  `settings.set_figure_params` to set the 'scanpy style'
      
- default cache directory is ``./cache/``, set `settings.cachedir` to change
  this; nested directories in this are now avoided

- show edges in scatter plots based on graph visualization
  :func:`~scanpy.api.tl.draw_graph` and :func:`~scanpy.api.umap` by passing
  `edges=True`

- :func:`~scanpy.api.pp.downsample_counts` function

- default 'louvain_groups' are now called 'louvain'

- 'X_diffmap' now contains the zero component, plotting remains unchanged
     
  

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
