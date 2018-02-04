Examples
--------

Good starting points are the following examples, which build on established results from the literature. All examples are versioned on `GitHub <scanpy_usage_>`__.

.. _scanpy_usage: https://github.com/theislab/scanpy_usage


------------

Clustering
~~~~~~~~~~

`Seurat's <Seurat_>`__ [Satija15]_ guided clustering `tutorial <http://satijalab.org/seurat/pbmc3k_tutorial.html>`_ for 2700 `PBMCs <https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell>`_ from 10x Genomics. This consists in preprocessing, clustering and the identification of cell types via known marker genes.

- `overview <17-05-05_>`_ with benchmarks
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb>`_  

.. raw:: html

   <a href="https://github.com/theislab/scanpy_usage/tree/master/170505_seurat"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"></a><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/louvain.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/NKG7.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/violin.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/cell_types.png" style="width: 200px">

Compare this to preprocessing and clustering of 68K PBMCs as in Cell Ranger [Zheng17]_: `overview <17-05-03_>`__.
   

------------

Simple Pseudotime
~~~~~~~~~~~~~~~~~

The Diffusion Pseudotime (DPT) analyses of [Haghverdi16]_ for hematopoietic data of [Paul15]_ and [Moignard15]_.

- `overview <17-05-02_>`_
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`_ for [Paul15]_
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170501_moignard15/moignard15.ipynb>`_ for [Moignard15]_

.. raw:: html

   <a href="https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16"><img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/scatter.png" style="width: 350px"></a><img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/heatmap.png" style="width: 80px">

   
------------

Scaling Computations
~~~~~~~~~~~~~~~~~~~~

Visualizing and clustering 1.3M `neurons <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1M_neurons>`_ from 10x Genomics.

- `overview <17-05-22_>`_
- run the computation in a `script <https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/cluster.py>`_
- plot the results in a `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/plot.ipynb>`_  

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170522_visualizing_one_million_cells/tsne_1.3M.png" style="width: 120px">
   
------------

Simulating
~~~~~~~~~~

Simulating single cells using literature-curated gene regulatory networks [Wittmann09]_; here, myeloid differentiation [Krumsiek11]_.

- `overview <17-04-30_>`_

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/timeseries.png" style="width: 200px"><img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/tsne.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/draw_graph.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/diffmap.png" style="width: 100px">
   
------------

Images
~~~~~~

Pseudotime-based vs. deep-learning based reconstruction of cell cycle from image data [Eulenberg17]_.

- `overview <17-04-30_>`_

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170529_images/dpt_DNA_content.png" style="width: 200px">
   
------------

.. _17-04-30: https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11
.. _17-05-03: https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17
.. _17-05-02: https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16
.. _17-05-05: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
.. _17-05-22: https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells

.. _Guided Clustering tutorial: http://satijalab.org/seurat/pbmc3k_tutorial.html
.. _Seurat: http://satijalab.org/seurat
.. _favorably discussed: https://doi.org/10.1101/110668
.. _Monocle: http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/



User Examples
~~~~~~~~~~~~~

January 12, 2018: `Exploring the mouse cell atlas <https://github.com/dpcook/fun_analysis/blob/master/tabula_muris/mouse_atlas_scanpy.ipynb>`_ by `David P. Cook <https://twitter.com/DavidPCook>`_. Data by `Tabula Muris Consortium <https://www.biorxiv.org/content/early/2017/12/20/237446>`_.
