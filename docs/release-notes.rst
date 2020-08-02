Release Notes
=============

.. role:: small
.. role:: smaller

.. note::

   Also see the release notes of :mod:`anndata`.

.. include:: _links.rst
.. include:: _key_contributors.rst


Version 1.6
-----------

.. include:: release-latest.rst

Version 1.5
-----------

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


Version 1.4
-----------

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

1.4.4 :small:`2019-07-20`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: New functionality

- :mod:`scanpy.get` adds helper functions for extracting data in convenient formats :pr:`619` :smaller:`I Virshup`

.. rubric:: Bug fixes

- Stopped deprecations warnings from AnnData `0.6.22` :smaller:`I Virshup`

.. rubric:: Code design

- :func:`~scanpy.pp.normalize_total` gains param `exclude_highly_expressed`, and `fraction` is renamed to `max_fraction` with better docs :smaller:`A Wolf`

1.4.3 :small:`2019-05-14`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: Bug fixes

- :func:`~scanpy.pp.neighbors` correctly infers `n_neighbors` again from `params`, which was temporarily broken in `v1.4.2` :smaller:`I Virshup`

.. rubric:: Code design

- :func:`~scanpy.pp.calculate_qc_metrics` is single threaded by default for datasets under 300,000 cells -- allowing cached compilation :pr:`615` :smaller:`I Virshup`

1.4.2 :small:`2019-05-06`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: New functionality

- :func:`~scanpy.pp.combat` supports additional covariates which may include adjustment variables or biological condition :pr:`618` :smaller:`G Eraslan`
- :func:`~scanpy.pp.highly_variable_genes` has a `batch_key` option which performs HVG selection in each batch separately to avoid selecting genes that vary strongly across batches :pr:`622` :smaller:`G Eraslan`

.. rubric:: Bug fixes

- :func:`~scanpy.tl.rank_genes_groups` t-test implementation doesn't return NaN when variance is 0, also changed to scipy's implementation :pr:`621` :smaller:`I Virshup`
- :func:`~scanpy.tl.umap` with `init_pos='paga'` detects correct `dtype` :smaller:`A Wolf`
- :func:`~scanpy.tl.louvain` and :func:`~scanpy.tl.leiden` auto-generate `key_added=louvain_R` upon passing `restrict_to`, which was temporarily changed in `1.4.1` :smaller:`A Wolf`

.. rubric:: Code design

- :func:`~scanpy.pp.neighbors` and :func:`~scanpy.tl.umap` got rid of UMAP legacy code and introduced UMAP as a dependency :pr:`576` :smaller:`S Rybakov`

1.4.1 :small:`2019-04-26`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: New functionality

- Scanpy has a command line interface again. Invoking it with `scanpy somecommand [args]` calls `scanpy-somecommand [args]`, except for builtin commands (currently `scanpy settings`)  :pr:`604` :smaller:`P Angerer`
- :func:`~scanpy.datasets.ebi_expression_atlas` allows convenient download of EBI expression atlas :smaller:`I Virshup`
- :func:`~scanpy.tl.marker_gene_overlap` computes overlaps of marker genes :smaller:`M Luecken`
- :func:`~scanpy.tl.filter_rank_genes_groups` filters out genes based on fold change and fraction of cells expressing genes :smaller:`F Ramirez`
- :func:`~scanpy.pp.normalize_total` replaces :func:`~scanpy.pp.normalize_per_cell`, is more efficient and provides a parameter to only normalize using a fraction of expressed genes :smaller:`S Rybakov`
- :func:`~scanpy.pp.downsample_counts` has been sped up, changed default value of `replace` parameter to `False`  :pr:`474` :smaller:`I Virshup`
- :func:`~scanpy.tl.embedding_density` computes densities on embeddings  :pr:`543` :smaller:`M Luecken`
- :func:`~scanpy.external.tl.palantir` interfaces Palantir [Setty18]_  :pr:`493` :smaller:`A Mousa`

.. rubric:: Code design

- `.layers` support of scatter plots :smaller:`F Ramirez`
- fix double-logarithmization in compute of log fold change in :func:`~scanpy.tl.rank_genes_groups` :smaller:`A Muñoz-Rojas`
- fix return sections of docs :smaller:`P Angerer`


Version 1.3
-----------

1.3.8 :small:`2019-02-05`
~~~~~~~~~~~~~~~~~~~~~~~~~
- :func:`~scanpy.read_10x_h5` throws more stringent errors and doesn’t require speciying default genomes anymore. :pr:`442` and :pr:`444` :smaller:`I Vishrup`

1.3.7 :small:`2019-01-02`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: Major updates

- one can `import scanpy as sc` instead of `import scanpy.api as sc`, see :mod:`scanpy`

.. rubric:: New functionality

- :func:`~scanpy.pp.combat` reimplements Combat for batch effect correction [Johnson07]_ [Leek12]_, heavily based on the Python implementation of [Pedersen12]_, but with performance improvements :pr:`398` :smaller:`M Lange`
- :func:`~scanpy.external.tl.phenograph` wraps the graph clustering package Phenograph [Levine15]_ :smaller:`A Mousa`

1.3.6 :small:`2018-12-11`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: Major updates

- a new plotting gallery for :doc:`visualizing-marker-genes` :smaller:`F Ramirez`
- tutorials are integrated on ReadTheDocs, :doc:`pbmc3k` and :doc:`paga-paul15` :smaller:`A Wolf`

.. rubric:: Interactive exploration of analysis results through *manifold viewers*

- CZI’s cellxgene_ directly reads `.h5ad` files :smaller:`the cellxgene developers`
- the `UCSC Single Cell Browser`_ requires exporting via :func:`~scanpy.external.exporting.cellbrowser` :smaller:`M Haeussler`

.. _cellxgene: https://github.com/chanzuckerberg/cellxgene
.. _UCSC Single Cell Browser: https://github.com/maximilianh/cellBrowser

.. rubric:: Code design

- :func:`~scanpy.pp.highly_variable_genes` supersedes :func:`~scanpy.pp.filter_genes_dispersion`, it gives the same results but, by default, expects logarithmized data and doesn’t subset :smaller:`A Wolf`

1.3.5 :small:`2018-12-09`
~~~~~~~~~~~~~~~~~~~~~~~~~

- uncountable figure improvements :pr:`369` :smaller:`F Ramirez`

1.3.4 :small:`2018-11-24`
~~~~~~~~~~~~~~~~~~~~~~~~~

- :func:`~scanpy.tl.leiden` wraps the recent graph clustering package by [Traag18]_ :smaller:`K Polanski`
- :func:`~scanpy.external.pp.bbknn` wraps the recent batch correction package [Polanski19]_ :smaller:`K Polanski`
- :func:`~scanpy.pp.calculate_qc_metrics` caculates a number of quality control metrics, similar to `calculateQCMetrics` from *Scater* [McCarthy17]_ :smaller:`I Virshup`

1.3.3 :small:`2018-11-05`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. rubric:: Major updates

- a fully distributed preprocessing backend :smaller:`T White and the Laserson Lab`

.. rubric:: Code design

- :func:`~scanpy.read_10x_h5` and :func:`~scanpy.read_10x_mtx` read Cell Ranger 3.0 outputs :pr:`334` :smaller:`Q Gong`

.. note::

   .. rubric:: Also see changes in anndata 0.6.

   - changed default compression to `None` in :meth:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical
   - performance gains in :meth:`~anndata.AnnData.write_h5ad` due to better handling of strings and categories :smaller:`S Rybakov`

1.3.1 :small:`2018-09-03`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: RNA velocity in single cells [Manno18]_

- Scanpy and AnnData support loom’s layers so that computations for single-cell RNA velocity [Manno18]_ become feasible :smaller:`S Rybakov and V Bergen`
- scvelo_ harmonizes with Scanpy and is able to process loom files with splicing information produced by Velocyto [Manno18]_, it runs a lot faster than the count matrix analysis of Velocyto and provides several conceptual developments

.. _scvelo: https://github.com/theislab/scvelo

.. rubric:: Plotting (:ref:`pl-generic`)

- :func:`~scanpy.pl.dotplot` for visualizing genes across conditions and clusters, see `here`__ :pr:`199` :smaller:`F Ramirez`
- :func:`~scanpy.pl.heatmap` for pretty heatmaps :pr:`175` :smaller:`F Ramirez`
- :func:`~scanpy.pl.violin` produces very compact overview figures with many panels :pr:`175` :smaller:`F Ramirez`

.. __: https://gist.github.com/fidelram/2289b7a8d6da055fb058ac9a79ed485c

.. rubric:: There now is a section on imputation in :doc:`external <external/index>`:

- :func:`~scanpy.external.pp.magic` for imputation using data diffusion [vanDijk18]_ :pr:`187` :smaller:`S Gigante`
- :func:`~scanpy.external.pp.dca` for imputation and latent space construction using an autoencoder [Eraslan18]_ :pr:`186` :smaller:`G Eraslan`


Version 1.2
-----------

1.2.1 :small:`2018-06-08`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Plotting of :ref:`pl-generic` marker genes and quality control.

- :func:`~scanpy.pl.highest_expr_genes` for quality control; plot genes with highest mean fraction of cells, similar to `plotQC` of *Scater* [McCarthy17]_ :pr:`169` :smaller:`F Ramirez`

1.2.0 :small:`2018-06-08`
~~~~~~~~~~~~~~~~~~~~~~~~~

- :func:`~scanpy.tl.paga` improved, see PAGA_; the default model changed, restore the previous default model by passing `model='v1.0'`


Version 1.1
-----------

1.1.0 :small:`2018-06-01`
~~~~~~~~~~~~~~~~~~~~~~~~~

- :func:`~scanpy.set_figure_params` by default passes `vector_friendly=True` and allows you to produce reasonablly sized pdfs by rasterizing large scatter plots :smaller:`A Wolf`
- :func:`~scanpy.tl.draw_graph` defaults to the ForceAtlas2 layout [Jacomy14]_ [Chippada18]_, which is often more visually appealing and whose computation is much faster :smaller:`S Wollock`
- :func:`~scanpy.pl.scatter` also plots along variables axis :smaller:`MD Luecken`
- :func:`~scanpy.pp.pca` and :func:`~scanpy.pp.log1p` support chunk processing :smaller:`S Rybakov`
- :func:`~scanpy.pp.regress_out` is back to multiprocessing :smaller:`F Ramirez`
- :func:`~scanpy.read` reads compressed text files :smaller:`G Eraslan`
- :func:`~scanpy.queries.mitochondrial_genes` for querying mito genes :smaller:`FG Brundu`
- :func:`~scanpy.external.pp.mnn_correct` for batch correction [Haghverdi18]_ [Kang18]_
- :func:`~scanpy.external.tl.phate` for low-dimensional embedding [Moon17]_ :smaller:`S Gigante`
- :func:`~scanpy.external.tl.sandbag`, :func:`~scanpy.external.tl.cyclone` for scoring genes [Scialdone15]_ [Fechtner18]_


Version 1.0
-----------

1.0.0 :small:`2018-03-30`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Major updates

- Scanpy is much faster and more memory efficient: preprocess, cluster and
  visualize 1.3M cells in 6h_, 130K cells in 14min_, and 68K cells in 3min_ :smaller:`A Wolf`
- the API gained a preprocessing function :func:`~scanpy.pp.neighbors` and a
  class :func:`~scanpy.Neighbors` to which all basic graph computations are
  delegated :smaller:`A Wolf`

.. _6h: https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/
.. _14min: https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/logfile_130K.txt
.. _3min: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170503_zheng17/zheng17.ipynb


.. warning::

   .. rubric:: Upgrading to 1.0 isn’t fully backwards compatible in the following changes

   - the graph-based tools :func:`~scanpy.tl.louvain`
     :func:`~scanpy.tl.dpt` :func:`~scanpy.tl.draw_graph`
     :func:`~scanpy.tl.umap` :func:`~scanpy.tl.diffmap`
     :func:`~scanpy.tl.paga` require prior computation of the graph:
     ``sc.pp.neighbors(adata, n_neighbors=5); sc.tl.louvain(adata)`` instead of
     previously ``sc.tl.louvain(adata, n_neighbors=5)``
   - install `numba` via ``conda install numba``, which replaces cython
   - the default connectivity measure (dpt will look different using default
     settings) changed. setting `method='gauss'` in `sc.pp.neighbors` uses
     gauss kernel connectivities and reproduces the previous behavior,
     see, for instance in the example paul15_.
   - namings of returned annotation have changed for less bloated AnnData
     objects, which means that some of the unstructured annotation of old
     AnnData files is not recognized anymore
   - replace occurances of `group_by` with `groupby` (consistency with
     `pandas`)
   - it is worth checking out the notebook examples to see changes, e.g.
     the seurat_ example.
   - upgrading scikit-learn from 0.18 to 0.19 changed the implementation of PCA,
     some results might therefore look slightly different

.. _paul15: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb
.. _seurat: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb

.. rubric:: Further updates

- UMAP [McInnes18]_ can serve as a first visualization of the data just as tSNE,
  in contrast to tSNE, UMAP directly embeds the single-cell graph and is faster;
  UMAP is also used for measuring connectivities and computing neighbors,
  see :func:`~scanpy.pp.neighbors` :smaller:`A Wolf`
- graph abstraction: AGA is renamed to PAGA_: :func:`~scanpy.tl.paga`; now,
  it only measures connectivities between partitions of the single-cell graph,
  pseudotime and clustering need to be computed separately via
  :func:`~scanpy.tl.louvain` and :func:`~scanpy.tl.dpt`, the
  connectivity measure has been improved :smaller:`A Wolf`
- logistic regression for finding marker genes
  :func:`~scanpy.tl.rank_genes_groups` with parameter `method='logreg'` :smaller:`A Wolf`
- :func:`~scanpy.tl.louvain` provides a better implementation for
  reclustering via `restrict_to` :smaller:`A Wolf`
- scanpy no longer modifies rcParams upon import, call
  `settings.set_figure_params` to set the 'scanpy style' :smaller:`A Wolf`
- default cache directory is ``./cache/``, set `settings.cachedir` to change
  this; nested directories in this are avoided :smaller:`A Wolf`
- show edges in scatter plots based on graph visualization
  :func:`~scanpy.tl.draw_graph` and :func:`~scanpy.tl.umap` by passing `edges=True` :smaller:`A Wolf`
- :func:`~scanpy.pp.downsample_counts` for downsampling counts :smaller:`MD Luecken`
- default `'louvain_groups'` are called `'louvain'` :smaller:`A Wolf`
- `'X_diffmap'` contains the zero component, plotting remains unchanged :smaller:`A Wolf`


Version 0.4
-----------

0.4.4 :small:`2018-02-26`
~~~~~~~~~~~~~~~~~~~~~~~~~

- embed cells using :func:`~scanpy.tl.umap` [McInnes18]_ :pr:`92` :smaller:`G Eraslan`
- score sets of genes, e.g. for `cell cycle`_, using :func:`~scanpy.tl.score_genes` [Satija15]_ :smaller:`D Cittaro`

.. _cell cycle: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb

0.4.3 :small:`2018-02-09`
~~~~~~~~~~~~~~~~~~~~~~~~~

- :func:`~scanpy.pl.clustermap`: heatmap from hierarchical clustering,
  based on :func:`seaborn.clustermap` [Waskom16]_ :smaller:`A Wolf`
- only return :class:`matplotlib.axes.Axes` in plotting functions of `sc.pl`
  when `show=False`, otherwise `None` :smaller:`A Wolf`

0.4.2 :small:`2018-01-07`
~~~~~~~~~~~~~~~~~~~~~~~~~

- amendments in PAGA_ and its plotting functions :smaller:`A Wolf`

0.4.0 :small:`2017-12-23`
~~~~~~~~~~~~~~~~~~~~~~~~~

- export to SPRING_ [Weinreb17]_ for interactive visualization of data:
  `spring tutorial`_ :smaller:`S Wollock`

.. _SPRING: https://github.com/AllonKleinLab/SPRING/
.. _spring tutorial: https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export


Version 0.3
-----------

0.3.2 :small:`2017-11-29`
~~~~~~~~~~~~~~~~~~~~~~~~~

- finding marker genes via :func:`~scanpy.pl.rank_genes_groups_violin` improved,
  see :issue:`51` :smaller:`F Ramirez`

0.3.0 :small:`2017-11-16`
~~~~~~~~~~~~~~~~~~~~~~~~~

- :class:`~anndata.AnnData` gains method :meth:`~anndata.AnnData.concatenate` :smaller:`A Wolf`
- :class:`~anndata.AnnData` is available as the separate anndata_ package :smaller:`P Angerer, A Wolf`
- results of PAGA_ simplified :smaller:`A Wolf`

.. _anndata: https://pypi.org/project/anndata/


Version 0.2
-----------

0.2.9 :small:`2017-10-25`
~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Initial release of the new trajectory inference method PAGA_

- :func:`~scanpy.tl.paga` computes an abstracted, coarse-grained (PAGA) graph of the neighborhood graph :smaller:`A Wolf`
- :func:`~scanpy.pl.paga_compare` plot this graph next an embedding :smaller:`A Wolf`
- :func:`~scanpy.pl.paga_path` plots a heatmap through a node sequence in the PAGA graph :smaller:`A Wolf`

0.2.1 :small:`2017-07-24`
~~~~~~~~~~~~~~~~~~~~~~~~~

Scanpy includes preprocessing, visualization, clustering, pseudotime and
trajectory inference, differential expression testing and simulation of gene
regulatory networks. The implementation efficiently deals with `datasets of more
than one million cells
<https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells>`__. :smaller:`A Wolf, P Angerer`


Version 0.1
-----------

0.1.0 :small:`2017-05-17`
~~~~~~~~~~~~~~~~~~~~~~~~~

Scanpy computationally outperforms and allows reproducing both the `Cell Ranger
R kit's <https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17>`__
and most of `Seurat’s
<https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`__
clustering workflows. :smaller:`A Wolf, P Angerer`
