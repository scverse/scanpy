.. note::

    Also see the `release notes <https://anndata.readthedocs.io>`__ of :mod:`anndata`.

.. role:: small

.. role:: smaller

On master :small:`July 31, 2018`
--------------------------------

Plotting of marker genes and quality control:

- :func:`~scanpy.api.pl.dotplot` for visualizing genes across conditions and clusters, see `here <https://gist.github.com/fidelram/2289b7a8d6da055fb058ac9a79ed485c>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.heatmap` for pretty heatmaps, see `here <https://github.com/theislab/scanpy/pull/175>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.violin` now produces very compact overview figures with many panels, see `here <https://github.com/theislab/scanpy/pull/175>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.highest_expr_genes` for quality control, see `here <https://github.com/theislab/scanpy/pull/169>`__; plot genes with highest mean fraction of cells, similar to `plotQC` of *Scater* [McCarthy17]_ :smaller:`thanks to F Ramirez`

There now is a `section <https://scanpy.readthedocs.io/en/latest/api/#imputation>`__ on imputation:

- :func:`~scanpy.api.pp.magic` for imputation using data diffusion [vanDijk18]_ :smaller:`thanks to S Gigante`
- :func:`~scanpy.api.pp.dca` for imputation and latent space construction using an autoencoder [Eraslan18]_

Further changes:

- `frameon=False` enables easy removal of frames in scatter plots and in :func:`~scanpy.api.set_figure_params`


Version 1.2 :small:`June 8, 2018`
---------------------------------

- :func:`~scanpy.api.tl.paga` improved, see `theislab/paga <https://github.com/theislab/paga>`__; the default model changed, restore the previous default model by passing `model='v1.0'`


Version 1.1 :small:`May 31, 2018`
---------------------------------

- :func:`~scanpy.api.set_figure_params` by default passes `vector_friendly=True` and allows you to produce reasonablly sized pdfs by rasterizing large scatter plots
- :func:`~scanpy.api.tl.draw_graph` now defaults to the ForceAtlas2 layout [Jacomy14]_ [Chippada18]_, which is often more visually appealing and whose computation is much faster :smaller:`thanks to S Wollock`
- :func:`~scanpy.api.pl.scatter` also plots along variables axis :smaller:`thanks to MD Luecken`
- :func:`~scanpy.api.pp.pca` and :func:`~scanpy.api.pp.log1p` support chunk processing :smaller:`thanks to S Rybakov`
- :func:`~scanpy.api.pp.regress_out` is back to multiprocessing :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.read` reads compressed text files :smaller:`thanks to G Eraslan`
- :func:`~scanpy.api.queries.mitochondrial_genes` for querying mito genes :smaller:`thanks to FG Brundu`
- :func:`~scanpy.api.pp.mnn_correct` for batch correction [Haghverdi18]_ [Kang18]_
- :func:`~scanpy.api.tl.phate` for low-dimensional embedding [Moon17]_ :smaller:`thanks to S Gigante`
- :func:`~scanpy.api.tl.sandbag`, :func:`~scanpy.api.tl.cyclone` for scoring genes [Scialdone15]_ [Fechtner18]_


Version 1.0 :small:`March 28, 2018`
-----------------------------------

Scanpy is much faster and more memory efficient. Preprocess, cluster and visualize
1.3M cells in `6 h
<https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/>`__,
130K cells in `14 min
<https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/logfile_130K.txt>`__
and 68K cells in `3 min
<https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170503_zheng17/zheng17.ipynb>`__.

The API gained a preprocessing function :func:`~scanpy.api.pp.neighbors` and a
class :func:`~scanpy.api.Neighbors` to which all basic graph computations are
delegated.

Upgrading to 1.0 isn't fully backwards compatible in the following changes:

- the graph-based tools :func:`~scanpy.api.tl.louvain`
  :func:`~scanpy.api.tl.dpt` :func:`~scanpy.api.tl.draw_graph`
  :func:`~scanpy.api.tl.umap` :func:`~scanpy.api.tl.diffmap`
  :func:`~scanpy.api.tl.paga` now require prior computation of the graph:
  ``sc.pp.neighbors(adata, n_neighbors=5); sc.tl.louvain(adata)`` instead of
  previously ``sc.tl.louvain(adata, n_neighbors=5)``
- install `numba` via ``conda install numba``, which replaces cython
- the default connectivity measure (dpt will look different using default
  settings) changed. setting `method='gauss'` in `sc.pp.neighbors` uses
  gauss kernel connectivities and reproduces the previous behavior,
  see, for instance this `example
  <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`__
- namings of returned annotation have changed for less bloated AnnData
  objects, which means that some of the unstructured annotation of old
  AnnData files is not recognized anymore
- replace occurances of `group_by` with `groupby` (consistency with
  `pandas`)
- it is worth checking out the notebook examples to see changes, e.g., `here
  <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb>`__
- upgrading scikit-learn from 0.18 to 0.19 changed the implementation of PCA,
  some results might therefore look slightly different

Further changes are:

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
- :func:`~scanpy.api.pp.downsample_counts` for downsampling counts :smaller:`thanks to MD Luecken`
- default 'louvain_groups' are now called 'louvain'
- 'X_diffmap' now contains the zero component, plotting remains unchanged


Version 0.4.4 :small:`February 26, 2018`
----------------------------------------

- embed cells using :func:`~scanpy.api.tl.umap` [McInnes18]_: `examples <https://github.com/theislab/scanpy/pull/92>`__
- score sets of genes, e.g. for cell cycle, using :func:`~scanpy.api.tl.score_genes` [Satija15]_: `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb>`__


Version 0.4.3 :small:`February 9, 2018`
---------------------------------------

- :func:`~scanpy.api.pl.clustermap`: heatmap from hierarchical clustering,
  based on :func:`seaborn.clustermap` [Waskom16]_
- only return `matplotlib.Axis` in plotting functions of ``sc.pl`` when `show=False`, otherwise `None`


Version 0.4.2 :small:`January 7, 2018`
--------------------------------------

- amendments in `PAGA <https://github.com/theislab/paga>`__ and its plotting
  functions


Version 0.4 :small:`December 23, 2017`
--------------------------------------

- export to `SPRING <https://github.com/AllonKleinLab/SPRING/>`__ [Weinreb17]_
  for interactive visualization of data: `tutorial
  <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`__,
  `docs <https://scanpy.readthedocs.io/en/latest/api/index.html>`__


Version 0.3.2 :small:`November 29, 2017`
----------------------------------------

- finding marker genes via :func:`~scanpy.api.pl.rank_genes_groups_violin` improved: `example <https://github.com/theislab/scanpy/issues/51>`__


Version 0.3 :small:`November 16, 2017`
--------------------------------------

- :class:`~anndata.AnnData` can be :meth:`~anndata.AnnData.concatenate` d.
- :class:`~anndata.AnnData` is available as a `separate package <https://pypi.org/project/anndata/>`__
- results of PAGA are `simplified <https://github.com/theislab/paga>`__


Version 0.2.9 :small:`October 25, 2017`
---------------------------------------

Initial release of `partition-based graph abstraction (PAGA) <https://github.com/theislab/paga>`__.


Version 0.2.1 :small:`July 24, 2017`
---------------------------------------

Scanpy now includes preprocessing, visualization, clustering, pseudotime and
trajectory inference, differential expression testing and simulation of gene
regulatory networks. The implementation efficiently deals with datasets of more
than one million cells.


Version 0.1 :small:`May 1, 2017`
--------------------------------

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing
most of Seurat's guided clustering tutorial.
