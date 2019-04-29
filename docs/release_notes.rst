.. note::

    Also see the `release notes <https://anndata.readthedocs.io>`__ of :mod:`anndata`.

.. role:: small
.. role:: smaller
.. role:: noteversion


On master :small:`April 27, 2019`
---------------------------------

Bug fixes:

- :func:`~scanpy.tl.rank_genes_groups` t-test implementation doesn't return NaN when variance is 0, also changed to scipy's implementation, see `PR <https://github.com/theislab/scanpy/pull/621>`__ :smaller:`thanks to I Virshup`
- :func:`~scanpy.pp.combat` ComBat function now supports additional covariates which may include adjustment variables or biological condition, see `PR <https://github.com/theislab/scanpy/pull/618>`__ :smaller:`thanks to G Eraslan`

Version 1.4.1 :small:`April 27, 2019`
-------------------------------------

New functionality:

- Scanpy has a command line interface again. Invoking it like ``scanpy somecommand [args]`` simply calls ``scanpy-somecommand [args]``, except for builting commands (currently just ``scanpy settings``). Implementation `here <https://github.com/theislab/scanpy/pull/604>`__.
- :func:`~scanpy.datasets.ebi_expression_atlas` allows convenient download of EBI expression atlas :smaller:`thanks to I Virshup`
- :func:`~scanpy.tl.marker_gene_overlap` computes overlaps of marker genes :smaller:`thanks to M Luecken`
- :func:`~scanpy.tl.filter_rank_genes_groups` filters out genes based on fold change and fraction of cells expressing genes :smaller:`thanks to F Ramirez`
- :func:`~scanpy.pp.normalize_total` replaces :func:`~scanpy.pp.normalize_per_cell`, is more efficient and provides a parameter to only normalize using a fraction of expressed genes :smaller:`thanks to S Rybakov`
- :func:`~scanpy.pp.downsample_counts` has been sped up, changed default value of `replace` parameter to `False`, see `here <https://github.com/theislab/scanpy/pull/474>`__ :smaller:`thanks to I Virshup`
- :func:`~scanpy.pl.embedding_density` allows plots of cell densities on embeddings, see `here <https://github.com/theislab/scanpy/pull/543>`__ :smaller:`thanks to M Luecken`
- :func:`~scanpy.external.palantir` interfaces Palantir [Setty18]_, see `here <https://github.com/theislab/scanpy/pull/493>`__ :smaller:`thanks to A Mousa`

Updates:

- `.layers` support of scatter plots :smaller:`thanks to F Ramirez`
- fix double-logarithmization in compute of log fold change in :func:`~scanpy.tl.rank_genes_groups` :smaller:`thanks to A Muñoz-Rojas`
- fix return sections of docs :smaller:`thanks to P Angerer`


Version 1.4 :small:`February 5, 2019`
-------------------------------------

Major updates:

- one can now `import scanpy as sc` instead of `import scanpy.api as sc`, see `here <https://scanpy.readthedocs.io/en/latest/api/>`__ :noteversion:`1.3.7`
- a new plotting gallery for visualizing marker genes, see `here <https://scanpy-tutorials.readthedocs.io/en/latest/visualizing-marker-genes.html>`__ :noteversion:`1.3.6` :smaller:`thanks to F Ramirez`
- tutorials are integrated on ReadTheDocs, see simple `clustering <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>`__ and simple `trajectory inference <https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html>`__ :noteversion:`1.3.6`
- a fully distributed preprocessing backend :noteversion:`1.3.3` :smaller:`thanks to T White and the Laserson Lab`
- changed default compression to `None` in :func:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical :noteversion:`anndata 0.6.16`
- performance gains in :func:`~anndata.AnnData.write_h5ad` due to better handling of strings and categories :noteversion:`anndata 0.6.19` :smaller:`thanks to S Rybakov`

Two new possibilities for interactive exploration of analysis results:

- CZI's `cellxgene <https://github.com/chanzuckerberg/cellxgene>`__ directly reads `.h5ad` files :smaller:`thanks to the cellxgene developers`
- the `UCSC Single Cell Browser <https://github.com/maximilianh/cellBrowser>`__ requires exporting via :func:`~scanpy.external.exporting.cellbrowser` :noteversion:`1.3.6` :smaller:`thanks to M Haeussler`

Further updates:

- :func:`~scanpy.pp.highly_variable_genes` supersedes :func:`~scanpy.pp.filter_genes_dispersion`, it gives the same results but, by default, expects logarithmized data and doesn't subset :noteversion:`1.3.6` :smaller:`thanks to S Rybakov`
- :func:`~scanpy.pp.combat` reimplements Combat for batch effect correction [Johnson07]_ [Leek12]_, heavily based on the Python implementation of [Pedersen12]_, but with performance improvements, see `here <https://github.com/theislab/scanpy/pull/398>`__ :noteversion:`1.3.7` :smaller:`thanks to M Lange`
- :func:`~scanpy.tl.leiden` wraps the recent graph clustering package by [Traag18]_ :noteversion:`1.3.4` :smaller:`thanks to K Polanski`
- :func:`~scanpy.external.pp.bbknn` wraps the recent batch correction package [Park18]_ :noteversion:`1.3.4` :smaller:`thanks to K Polanski`
- :func:`~scanpy.external.tl.phenograph` wraps the graph clustering package Phenograph [Levine15]_  :noteversion:`1.3.7` :smaller:`thanks to A Mousa`
- :func:`~scanpy.pp.calculate_qc_metrics` caculates a number of quality control metrics, similar to `calculateQCMetrics` from *Scater* [McCarthy17]_ :noteversion:`1.3.4` :smaller:`thanks to I Virshup`
- :func:`~scanpy.read_10x_h5` throws more stringent errors and doesn't require speciying default genomes anymore, see `here <https://github.com/theislab/scanpy/pull/442>`__ and `here <https://github.com/theislab/scanpy/pull/444>`__ :noteversion:`1.3.8`  :smaller:`thanks to I Vishrup`
- :func:`~scanpy.read_10x_h5` and :func:`~scanpy.read_10x_mtx` read Cell Ranger 3.0 outputs, see `here <https://github.com/theislab/scanpy/pull/334>`__ :noteversion:`1.3.3`  :smaller:`thanks to Q Gong`


Version 1.3 :small:`September 3, 2018`
--------------------------------------

RNA velocity in single cells [Manno18]_:

- Scanpy and AnnData support loom's layers so that computations for single-cell RNA velocity [Manno18]_ become feasible :smaller:`thanks to S Rybakov and V Bergen`
- the package `scvelo <https://github.com/theislab/scvelo>`__ perfectly harmonizes with Scanpy and is able to process loom files with splicing information produced by Velocyto [Manno18]_, it runs a lot faster than the count matrix analysis of Velocyto and provides several conceptual developments (preprint to come)

Plotting of marker genes and quality control, see this `section <https://scanpy.readthedocs.io/en/latest/api/plotting.html#generic>`__ and scroll down, a few examples are

- :func:`~scanpy.api.pl.dotplot` for visualizing genes across conditions and clusters, see `here <https://gist.github.com/fidelram/2289b7a8d6da055fb058ac9a79ed485c>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.heatmap` for pretty heatmaps, see `here <https://github.com/theislab/scanpy/pull/175>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.violin` now produces very compact overview figures with many panels, see `here <https://github.com/theislab/scanpy/pull/175>`__ :smaller:`thanks to F Ramirez`
- :func:`~scanpy.api.pl.highest_expr_genes` for quality control, see `here <https://github.com/theislab/scanpy/pull/169>`__; plot genes with highest mean fraction of cells, similar to `plotQC` of *Scater* [McCarthy17]_ :smaller:`thanks to F Ramirez`

There is a `section <https://scanpy.readthedocs.io/en/latest/api/#imputation>`__ on imputation:

- :func:`~scanpy.api.pp.magic` for imputation using data diffusion [vanDijk18]_ :smaller:`thanks to S Gigante`
- :func:`~scanpy.api.pp.dca` for imputation and latent space construction using an autoencoder [Eraslan18]_


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
- :func:`~scanpy.external.pp.mnn_correct` for batch correction [Haghverdi18]_ [Kang18]_
- :func:`~scanpy.external.tl.phate` for low-dimensional embedding [Moon17]_ :smaller:`thanks to S Gigante`
- :func:`~scanpy.external.tl.sandbag`, :func:`~scanpy.api.tl.cyclone` for scoring genes [Scialdone15]_ [Fechtner18]_


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
