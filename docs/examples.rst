Examples
--------

To explore different versions of the example notebooks, use their history on GitHub.

------------

Clustering
~~~~~~~~~~

`Seurat's <http://satijalab.org/seurat>`_ [Satija15]_ guided clustering `tutorial <http://satijalab.org/seurat/pbmc3k_tutorial.html>`_ for 2700 `PBMCs <https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell>`_ from 10x Genomics.

- `summary <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`_ with benchmarks
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb>`_ containing preprocessing, clustering and the identification of cell types via known marker genes

.. raw:: html

   <a href="https://github.com/theislab/scanpy_usage/tree/master/170505_seurat"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"></a><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/louvain.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/NKG7.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/violin.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/cell_types.png" style="width: 200px">

Compare this to `preprocessing and clustering of 68K PBMCs <https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17>`_ as in Cell Ranger [Zheng17]_.


------------

Simple Pseudotime
~~~~~~~~~~~~~~~~~

The Diffusion Pseudotime (DPT) analyses of [Haghverdi16]_ for hematopoietic data of [Paul15]_ and [Moignard15]_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/scatter.png" style="width: 350px; margin: -30px 0px 0px 0px" align="right"><img src="http://falexwolf.de/img/scanpy_usage/170501_moignard/heatmap.png" style="width: 80px; margin: -30px 0px 0px -150px" align="right">

- `summary <https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16>`_
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`_ for [Paul15]_
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170501_moignard15/moignard15.ipynb>`_ for [Moignard15]_


------------

Scaling Computations
~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170522_visualizing_one_million_cells/tsne_1.3M.png" style="width: 120px; margin: -15px 100px 0px 0px" align="right">

Visualizing and clustering 1.3M `neurons <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1M_neurons>`_ from 10x Genomics.

- `summary <https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells>`_
- `script <https://github.com/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/cluster.py>`_ for running the computation
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170522_visualizing_one_million_cells/plot.ipynb>`_ for plotting the results


------------

Simulations
~~~~~~~~~~~

Simulating single cells using literature-curated gene regulatory networks [Wittmann09]_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/timeseries.png" style="width: 200px; margin: -15px 0px 0px 0px" align="right"><img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/draw_graph.png" style="width: 100px; margin: -15px 0px 0px -100px" align="right">
  
- `summary <https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11>`_ for myeloid differentiation [Krumsiek11]_
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170430_krumsiek11/krumsiek11.ipynb>`_ for myeloid differentiation
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170430_krumsiek11/toggleswitch.ipynb>`_ for simple toggleswitch


------------

Images
~~~~~~

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170529_images/dpt_DNA_content.png" style="width: 200px" align="right">

Pseudotime-based vs. deep-learning based reconstruction of cell cycle from image data [Eulenberg17]_.

- `summary <https://github.com/theislab/scanpy_usage/tree/master/170529_images>`_


------------

User Examples
~~~~~~~~~~~~~

January 12, 2018: `Exploring the mouse cell atlas <https://github.com/dpcook/fun_analysis/blob/master/tabula_muris/mouse_atlas_scanpy.ipynb>`_ by `David P. Cook <https://twitter.com/DavidPCook>`_. Data by `Tabula Muris Consortium <https://www.biorxiv.org/content/early/2017/12/20/237446>`_.
