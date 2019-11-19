Tutorials
=========


----------
Clustering
----------

For getting started, we recommend `Scanpy's reimplementation <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>`__ of Seurat's [Satija15]_ clustering tutorial for 3K PBMCs from 10x Genomics, containing preprocessing, clustering and the identification of cell types via known marker genes.

For more possibilities on visualizing marker genes, see this `plotting gallery <https://scanpy-tutorials.readthedocs.io/en/latest/visualizing-marker-genes.html>`__.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/louvain.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/NKG7.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/violin.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/cell_types.png" style="width: 200px">


--------------------
Trajectory Inference
--------------------

For trajectory inference on complex datasets, we offer several examples `here <https://github.com/theislab/paga>`__. Get started `here <https://nbviewer.jupyter.org/github/theislab/paga/blob/master/blood/paul15/paul15.ipynb>`__ for the following result on hematopoiesis.

.. raw:: html

   <img src="http://www.falexwolf.de/img/paga_paul15.png" style="width: 450px">

You can extend this to multi-resolution analyses of whole animals, such as `here <https://nbviewer.jupyter.org/github/theislab/paga/blob/master/planaria/planaria.ipynb>`__.

.. raw:: html

   <img src="http://www.falexwolf.de/img/paga_planaria.png" style="width: 350px">

The PAGA method behind this is described in [Wolf19]_. As a reference for simple pseudotime analyses, we provide the diffusion pseudotime analyses of [Haghverdi16]_ for two hematopoiesis datasets: `here <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`__ for [Paul15]_ and `here <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170501_moignard15/moignard15.ipynb>`__ for [Moignard15]_.


-----------------
Further Tutorials
-----------------

Conversion: AnnData, SingleCellExperiment, and Seurat objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See this `notebook <https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb>`__
for a tutorial on `anndata2ri`.

Regressing out cell cycle
~~~~~~~~~~~~~~~~~~~~~~~~~

See this `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb>`__.


Scaling Computations
~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170522_visualizing_one_million_cells/tsne_1.3M.png" style="width: 120px; margin: -100px 50px 0px 0px" align="right">

Visualize and cluster 1.3M neurons from 10x Genomics `here <https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells>`__.


Simulations
~~~~~~~~~~~

Simulating single cells using literature-curated gene regulatory networks [Wittmann09]_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/timeseries.png" style="width: 200px; margin: -15px 0px 0px 0px" align="right"><img src="http://falexwolf.de/img/scanpy_usage/170430_krumsiek11/draw_graph.png" style="width: 100px; margin: -15px 0px 0px -100px" align="right">

- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170430_krumsiek11/krumsiek11.ipynb>`__ for myeloid differentiation
- `notebook <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170430_krumsiek11/toggleswitch.ipynb>`__ for simple toggleswitch


Images
~~~~~~

See a pseudotime-based vs. deep-learning based reconstruction of cell cycle from image data `here <https://github.com/theislab/scanpy_usage/tree/master/170529_images>`__ [Eulenberg17]_.


..
    User Examples
    ~~~~~~~~~~~~~

    January 12, 2018: `Exploring the mouse cell atlas <https://github.com/dpcook/fun_analysis/blob/master/tabula_muris/mouse_atlas_scanpy.ipynb>`__ by `David P. Cook <https://twitter.com/DavidPCook>`__. Data by `Tabula Muris Consortium <https://www.biorxiv.org/content/early/2017/12/20/237446>`__.

