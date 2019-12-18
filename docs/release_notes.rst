.. note::
   Also see the `release notes`__ of :mod:`anndata`.

.. __: https://anndata.readthedocs.io

.. role:: small
.. role:: smaller

.. sidebar:: Key Contributors

   `anndata graph`_ & `scanpy graph`_;
   the ☀ signifies current maintainers

   * `Isaac Virshup`_: anndata overhaul, diverse contributions ☀
   * Gökcen Eraslan: diverse contributions ☀
   * Sergei Rybakov: diverse contributions ☀
   * Fidel Ramirez: plotting ☀
   * `Tom White`_: distributed computing
   * Philipp Angerer: initial anndata conception/development, software quality ☀
   * `Alex Wolf`_: initial anndata & scanpy conception/development ☀
   * `Fabian Theis`_ & lab: enabling guidance, support and environment

.. _anndata graph: https://github.com/theislab/anndata/graphs/contributors
.. _scanpy graph: https://github.com/theislab/scanpy/graphs/contributors
.. _Isaac Virshup: https://twitter.com/ivirshup
.. _Tom White: https://twitter.com/tom_e_white
.. _Alex Wolf: https://twitter.com/falexwolf
.. _Fabian Theis: https://twitter.com/fabian_theis


On master
---------

- :mod:`scanpy.pp.downsample_counts` now always preserves the dtype of it's input, instead of converting to floats to int :noteversion:`1.4.5` :pr:`865` :smaller:`thanks to I Virshup`
- :mod:`scanpy.queries` recieved many updates. This includes enrichment through gprofiler_ and more advanced biomart queries :pr:`467` :smaller:`thanks to I Virshup`
- Allow specifying a base for :func:`~scanpy.pp.log1p` :pr:`931` :smaller:`thanks to G Eraslan`
- :func:`~scanpy.tl.ingest` allows to map labels and embeddings from one adata to another :pr:`651` :smaller:`thanks to S Rybakov`

.. _gprofiler: https://biit.cs.ut.ee/gprofiler/


Version 1.4.*
-------------

1.4.4
~~~~~
New functionality:

- :mod:`scanpy.get` adds helper functions for extracting data in convenient formats :noteversion:`1.4.4` :pr:`619` :smaller:`thanks to I Virshup`

Bug fixes:

- Stopped deprecations warnings from AnnData `0.6.22` :noteversion:`1.4.4` :smaller:`thanks to I Virshup`

Code design:

- :func:`~scanpy.pp.normalize_total` gains param `exclude_highly_expressed`, and `fraction` is renamed to `max_fraction` with better docs :noteversion:`1.4.4` :smaller:`thanks to A Wolf`

1.4.3
~~~~~
Bug fixes:

- :func:`~scanpy.pp.neighbors` correctly infers `n_neighbors` again from `params`, which was temporarily broken in `v1.4.2` :noteversion:`1.4.3` :smaller:`thanks to I Virshup`

Code design:

- :func:`~scanpy.pp.calculate_qc_metrics` is single threaded by default for datasets under 300,000 cells -- allowing cached compilation :noteversion:`1.4.3` :pr:`615 :smaller:`thanks to I Virshup`

1.4.2
~~~~~
New functionality:

- :func:`~scanpy.pp.combat` supports additional covariates which may include adjustment variables or biological condition :noteversion:`1.4.2` :pr:`618 :smaller:`thanks to G Eraslan`
- :func:`~scanpy.pp.highly_variable_genes` has a `batch_key` option which performs HVG selection in each batch separately to avoid selecting genes that vary strongly across batches :noteversion:`1.4.2` :pr:`622 :smaller:`thanks to G Eraslan`

Bug fixes:

- :func:`~scanpy.tl.rank_genes_groups` t-test implementation doesn't return NaN when variance is 0, also changed to scipy's implementation :noteversion:`1.4.2` :pr:`621 :smaller:`thanks to I Virshup`
- :func:`~scanpy.tl.umap` with `init_pos='paga'` detects correct `dtype` :noteversion:`1.4.2` :smaller:`thanks to A Wolf`
- :func:`~scanpy.tl.louvain` and :func:`~scanpy.tl.leiden` auto-generate `key_added=louvain_R` upon passing `restrict_to`, which was temporarily changed in `v1.4.1` :noteversion:`1.4.2` :smaller:`thanks to A Wolf`

Code design:

- :func:`~scanpy.pp.neighbors` and :func:`~scanpy.tl.umap` got rid of UMAP legacy code and introduced UMAP as a dependency :noteversion:`1.4.2` :pr:`576 :smaller:`thanks to S Rybakov`

1.4.1
~~~~~
New functionality:

- Scanpy has a command line interface again. Invoking it with `scanpy somecommand [args]` calls `scanpy-somecommand [args]`, except for builtin commands (currently `scanpy settings`) :noteversion:`1.4.1` :pr:`604` :smaller:`thanks to P Angerer`
- :func:`~scanpy.datasets.ebi_expression_atlas` allows convenient download of EBI expression atlas :noteversion:`1.4.1` :smaller:`thanks to I Virshup`
- :func:`~scanpy.tl.marker_gene_overlap` computes overlaps of marker genes :noteversion:`1.4.1` :smaller:`thanks to M Luecken`
- :func:`~scanpy.tl.filter_rank_genes_groups` filters out genes based on fold change and fraction of cells expressing genes :noteversion:`1.4.1` :smaller:`thanks to F Ramirez`
- :func:`~scanpy.pp.normalize_total` replaces :func:`~scanpy.pp.normalize_per_cell`, is more efficient and provides a parameter to only normalize using a fraction of expressed genes :noteversion:`1.4.1` :smaller:`thanks to S Rybakov`
- :func:`~scanpy.pp.downsample_counts` has been sped up, changed default value of `replace` parameter to `False` :noteversion:`1.4.1` :pr:`474 :smaller:`thanks to I Virshup`
- :func:`~scanpy.pl.embedding_density` allows plots of cell densities on embeddings :noteversion:`1.4.1` :pr:`543 :smaller:`thanks to M Luecken`
- :func:`~scanpy.external.tl.palantir` interfaces Palantir [Setty18]_ :noteversion:`1.4.1` :pr:`493 :smaller:`thanks to A Mousa`

Code design:

- `.layers` support of scatter plots :noteversion:`1.4.1` :smaller:`thanks to F Ramirez`
- fix double-logarithmization in compute of log fold change in :func:`~scanpy.tl.rank_genes_groups` :noteversion:`1.4.1` :smaller:`thanks to A Muñoz-Rojas`
- fix return sections of docs :noteversion:`1.4.1` :smaller:`thanks to P Angerer`


Version 1.3.*
-------------

1.3.8
~~~~~
- :func:`~scanpy.read_10x_h5` throws more stringent errors and doesn’t require speciying default genomes anymore. :noteversion:`1.3.8` :pr:`442` and :pr:`444  :smaller:`thanks to I Vishrup`

1.3.7
~~~~~
Major updates:

- one can `import scanpy as sc` instead of `import scanpy.api as sc`, see :mod:`scanpy` :noteversion:`1.3.7`

Further updates:

- :func:`~scanpy.pp.combat` reimplements Combat for batch effect correction [Johnson07]_ [Leek12]_, heavily based on the Python implementation of [Pedersen12]_, but with performance improvements, see :noteversion:`1.3.7` :pr:`398 :smaller:`thanks to M Lange`
- :func:`~scanpy.external.tl.phenograph` wraps the graph clustering package Phenograph [Levine15]_  :noteversion:`1.3.7` :smaller:`thanks to A Mousa`

1.3.6
~~~~~
Major updates:

- a new plotting gallery for :doc:`visualizing-marker-genes` :noteversion:`1.3.6` :smaller:`thanks to F Ramirez`
- tutorials are integrated on ReadTheDocs, :doc:`pbmc3k` and :doc:`paga-paul15` :noteversion:`1.3.6`

Two new possibilities for interactive exploration of analysis results:

- CZI’s cellxgene_ directly reads `.h5ad` files :smaller:`thanks to the cellxgene developers`
- the `UCSC Single Cell Browser`_ requires exporting via :func:`~scanpy.external.exporting.cellbrowser` :noteversion:`1.3.6` :smaller:`thanks to M Haeussler`

.. _cellxgene: https://github.com/chanzuckerberg/cellxgene
.. _UCSC Single Cell Browser: https://github.com/maximilianh/cellBrowser


Further updates:

- :func:`~scanpy.pp.highly_variable_genes` supersedes :func:`~scanpy.pp.filter_genes_dispersion`, it gives the same results but, by default, expects logarithmized data and doesn’t subset :noteversion:`1.3.6`

1.3.5
~~~~~

- Uncountable figure improvements :noteversion:`1.3.5` :pr:`369` :smaller:`thanks to F Ramirez`

1.3.4
~~~~~

- :func:`~scanpy.tl.leiden` wraps the recent graph clustering package by [Traag18]_ :noteversion:`1.3.4` :smaller:`thanks to K Polanski`
- :func:`~scanpy.external.pp.bbknn` wraps the recent batch correction package [Polanski19]_ :noteversion:`1.3.4` :smaller:`thanks to K Polanski`
- :func:`~scanpy.pp.calculate_qc_metrics` caculates a number of quality control metrics, similar to `calculateQCMetrics` from *Scater* [McCarthy17]_ :noteversion:`1.3.4` :smaller:`thanks to I Virshup`

1.3.3
~~~~~

Major updates:

- a fully distributed preprocessing backend :noteversion:`1.3.3` :smaller:`thanks to T White and the Laserson Lab`

Further updates:

- :func:`~scanpy.read_10x_h5` and :func:`~scanpy.read_10x_mtx` read Cell Ranger 3.0 outputs, see :noteversion:`1.3.3` :pr:`334` :smaller:`thanks to Q Gong`

AnnData 0.6.*
~~~~~~~~~~~~~

- changed default compression to `None` in :meth:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical :noteversion:`anndata 0.6.16`
- performance gains in :meth:`~anndata.AnnData.write_h5ad` due to better handling of strings and categories :noteversion:`anndata 0.6.19` :smaller:`thanks to S Rybakov`

1.3
~~~

RNA velocity in single cells [Manno18]_:

- Scanpy and AnnData support loom’s layers so that computations for single-cell RNA velocity [Manno18]_ become feasible :smaller:`thanks to S Rybakov and V Bergen`
- the package scvelo_ perfectly harmonizes with Scanpy and is able to process loom files with splicing information produced by Velocyto [Manno18]_, it runs a lot faster than the count matrix analysis of Velocyto and provides several conceptual developments (preprint to come)

.. _scvelo: https://github.com/theislab/scvelo

Plotting of :ref:`pl-generic` marker genes and quality control.

- :func:`~scanpy.pl.dotplot` for visualizing genes across conditions and clusters, see `here`__. :noteversion:`1.3` :pr:`199` :smaller:`thanks to F Ramirez`
- :func:`~scanpy.pl.heatmap` for pretty heatmaps. :noteversion:`1.3` :pr:`175` :smaller:`thanks to F Ramirez`
- :func:`~scanpy.pl.violin` produces very compact overview figures with many panels. :noteversion:`1.3` :pr:`175` :smaller:`thanks to F Ramirez`

.. __: https://gist.github.com/fidelram/2289b7a8d6da055fb058ac9a79ed485c

There is now a section on :ref:`pp-imputation`:

- :func:`~scanpy.external.pp.magic` for imputation using data diffusion [vanDijk18]_. :noteversion:`1.3` :pr:`187` :smaller:`thanks to S Gigante`
- :func:`~scanpy.external.pp.dca` for imputation and latent space construction using an autoencoder [Eraslan18]_. :noteversion:`1.3` :pr:`186` :smaller:`thanks to G Eraslan`


Version 1.2 :small:`June 8, 2018`
---------------------------------

1.2.1
~~~~~

Plotting of :ref:`pl-generic` marker genes and quality control.

- :func:`~scanpy.pl.highest_expr_genes` for quality control; plot genes with highest mean fraction of cells, similar to `plotQC` of *Scater* [McCarthy17]_. :noteversion:`1.2.1` :pr:`169` :smaller:`thanks to F Ramirez`

1.2
~~~

- :func:`~scanpy.tl.paga` improved, see `theislab/paga`_; the default model changed, restore the previous default model by passing `model='v1.0'`


Version 1.1 :small:`May 31, 2018`
---------------------------------

- :func:`~scanpy.set_figure_params` by default passes `vector_friendly=True` and allows you to produce reasonablly sized pdfs by rasterizing large scatter plots
- :func:`~scanpy.tl.draw_graph` defaults to the ForceAtlas2 layout [Jacomy14]_ [Chippada18]_, which is often more visually appealing and whose computation is much faster :smaller:`thanks to S Wollock`
- :func:`~scanpy.pl.scatter` also plots along variables axis :smaller:`thanks to MD Luecken`
- :func:`~scanpy.pp.pca` and :func:`~scanpy.pp.log1p` support chunk processing :smaller:`thanks to S Rybakov`
- :func:`~scanpy.pp.regress_out` is back to multiprocessing :smaller:`thanks to F Ramirez`
- :func:`~scanpy.read` reads compressed text files :smaller:`thanks to G Eraslan`
- :func:`~scanpy.queries.mitochondrial_genes` for querying mito genes :smaller:`thanks to FG Brundu`
- :func:`~scanpy.external.pp.mnn_correct` for batch correction [Haghverdi18]_ [Kang18]_
- :func:`~scanpy.external.tl.phate` for low-dimensional embedding [Moon17]_ :smaller:`thanks to S Gigante`
- :func:`~scanpy.external.tl.sandbag`, :func:`~scanpy.external.tl.cyclone` for scoring genes [Scialdone15]_ [Fechtner18]_


Version 1.0 :small:`March 28, 2018`
-----------------------------------

Scanpy is much faster and more memory efficient. Preprocess, cluster and visualize
1.3M cells in 6h_, 130K cells in 14min_, and 68K cells in 3min_.

.. _6h: https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/
.. _14min: https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/logfile_130K.txt
.. _3min: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170503_zheng17/zheng17.ipynb

The API gained a preprocessing function :func:`~scanpy.pp.neighbors` and a
class :func:`~scanpy.Neighbors` to which all basic graph computations are
delegated.

Upgrading to 1.0 isn’t fully backwards compatible in the following changes:

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

Further changes are:

- UMAP [McInnes18]_ can serve as a first visualization of the data just as tSNE,
  in contrast to tSNE, UMAP directly embeds the single-cell graph and is faster;
  UMAP is also used for measuring connectivities and computing neighbors,
  see :func:`~scanpy.pp.neighbors`
- graph abstraction: AGA is renamed to PAGA: :func:`~scanpy.tl.paga`; now,
  it only measures connectivities between partitions of the single-cell graph,
  pseudotime and clustering need to be computed separately via
  :func:`~scanpy.tl.louvain` and :func:`~scanpy.tl.dpt`, the
  connectivity measure has been improved
- logistic regression for finding marker genes
  :func:`~scanpy.tl.rank_genes_groups` with parameter `method='logreg'`
- :func:`~scanpy.tl.louvain` provides a better implementation for
  reclustering via `restrict_to`
- scanpy no longer modifies rcParams upon import, call
  `settings.set_figure_params` to set the 'scanpy style'
- default cache directory is ``./cache/``, set `settings.cachedir` to change
  this; nested directories in this are avoided
- show edges in scatter plots based on graph visualization
  :func:`~scanpy.tl.draw_graph` and :func:`~scanpy.tl.umap` by passing `edges=True`
- :func:`~scanpy.pp.downsample_counts` for downsampling counts :smaller:`thanks to MD Luecken`
- default `'louvain_groups'` are called `'louvain'`
- `'X_diffmap'` contains the zero component, plotting remains unchanged


Version 0.4.4 :small:`February 26, 2018`
----------------------------------------

- embed cells using :func:`~scanpy.tl.umap` [McInnes18]_: :pr:`92`
- score sets of genes, e.g. for `cell cycle`_, using :func:`~scanpy.tl.score_genes` [Satija15]_.

.. _cell cycle: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb

Version 0.4.3 :small:`February 9, 2018`
---------------------------------------

- :func:`~scanpy.pl.clustermap`: heatmap from hierarchical clustering,
  based on :func:`seaborn.clustermap` [Waskom16]_
- only return :class:`matplotlib.axes.Axes` in plotting functions of `sc.pl`
  when `show=False`, otherwise `None`


Version 0.4.2 :small:`January 7, 2018`
--------------------------------------

- amendments in `theislab/paga`_ and its plotting functions


Version 0.4 :small:`December 23, 2017`
--------------------------------------

- export to SPRING_ [Weinreb17]_ for interactive visualization of data:
  `spring tutorial`_, docs :mod:`scanpy.api`.

.. _SPRING: https://github.com/AllonKleinLab/SPRING/
.. _spring tutorial: https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export

Version 0.3.2 :small:`November 29, 2017`
----------------------------------------

- finding marker genes via :func:`~scanpy.pl.rank_genes_groups_violin` improved:
  For an example, see :issue:`51`.


Version 0.3 :small:`November 16, 2017`
--------------------------------------

- :class:`~anndata.AnnData` can be :meth:`~anndata.AnnData.concatenate` d.
- :class:`~anndata.AnnData` is available as the anndata_ package.
- results of PAGA are simplified: `theislab/paga`_

.. _anndata: https://pypi.org/project/anndata/


Version 0.2.9 :small:`October 25, 2017`
---------------------------------------

Initial release of *partition-based graph abstraction (PAGA)*: `theislab/paga`_

.. _theislab/paga: https://github.com/theislab/paga


Version 0.2.1 :small:`July 24, 2017`
---------------------------------------

Scanpy includes preprocessing, visualization, clustering, pseudotime and
trajectory inference, differential expression testing and simulation of gene
regulatory networks. The implementation efficiently deals with datasets of more
than one million cells.


Version 0.1 :small:`May 1, 2017`
--------------------------------

Scanpy computationally outperforms the Cell Ranger R kit and allows reproducing
most of Seurat’s guided clustering tutorial.
