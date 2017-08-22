Examples
--------

Good starting points are the following examples, which build on establish results from the literature. All use cases are collected on the `GitHub <scanpy_usage_>`__.

.. _scanpy_usage: https://github.com/theislab/scanpy_usage

------------

`Example 1: <17-05-05_>`__ `Seurat's <Seurat_>`__ [Satija15]_ guided clustering `tutorial <http://satijalab.org/seurat/pbmc3k_tutorial.html>`_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/louvain.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/NKG7.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/violin.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/cell_types.png" style="width: 200px">

------------

`Example 2: <17-05-02_>`__ The Diffusion Pseudotime (DPT) analyses as of [Haghverdi16]_. Note that DPT has recently been very `favorably discussed`_ by the authors of Monocle_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/scatter.png" style="width: 350px"><img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/heatmap.png" style="width: 80px">


------------

`Example 3: <17-05-03_>`__ Analyzing 68 000 cells from [Zheng17]_, we find that Scanpy is about a factor 5 to 16 faster and more memory efficient than the `Cell Ranger`_ R kit for secondary analysis.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170503_zheng17/speedup.png" style="width: 300px"><img src="http://falexwolf.de/img/scanpy_usage/170503_zheng17/scatter.png" style="width: 100px">
   
------------

.. _17-05-05: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
.. _17-05-03: https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17
.. _17-05-02: https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16
.. _17-04-30: https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11

.. _Guided Clustering tutorial: http://satijalab.org/seurat/pbmc3k_tutorial.html
.. _Seurat: http://satijalab.org/seurat
.. _Cell Ranger: https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
.. _favorably discussed: https://doi.org/10.1101/110668
.. _Monocle: http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/
