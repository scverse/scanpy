`Getting started`_ \| Features_ \| Installation_ \| References_

|Build Status|

.. |Build Status| image:: https://travis-ci.org/theislab/scanpy.svg?branch=master
   :target: https://travis-ci.org/theislab/scanpy

Scanpy -- Single-Cell Analysis in Python
========================================

Efficient tools for analyzing and simulating large-scale single-cell data that aim at an understanding of dynamic biological processes from snapshots of transcriptome or proteome. The draft `Wolf, Angerer & Theis (2017) <http://falexwolf.de/docs/scanpy.pdf>`__ explains conceptual ideas of the package. Any comments are appreciated!

Getting started
---------------

Download or clone the repository -- green button on top of the page -- and ``cd`` into its root directory. With Python 3.5 or 3.6 (preferably Miniconda_) installed, type::

    pip install -e .

Aside from enabling ``import scanpy as sc`` anywhere on your system, you can also work with the top-level command ``scanpy`` on the command-line (more info `here <Installation_>`__).

Then go through the use cases compiled in scanpy_usage_, in particular, the recent additions

.. _scanpy_usage: https://github.com/theislab/scanpy_usage

17-05-05_
  We reproduce some of the recent `Guided Clustering tutorial`_ of Seurat_ [macosco15]_.
17-05-03_
  Analyzing 64 000 cells from [zheng]_, we find that Scanpy is about a factor 5 to 10 faster and more memory efficient than the optimized `Cell Ranger`_ pipeline. For large-scale data, this becomes crucial for interactive analysis.
17-05-01_
  Diffusion Pseudotime analysis resolves developmental processes in data of [moignard15]_, reproducing results of [haghverdi16]_. Also, note that DPT has recently been very `favorably discussed`_ by the authors of Monocle_.

.. _17-05-05: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
.. _17-05-03: https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17
.. _17-05-01: https://github.com/theislab/scanpy_usage/tree/master/170501_moignard15/notebook.ipynb

.. _Guided Clustering tutorial: http://satijalab.org/seurat/pbmc-tutorial.html
.. _Seurat: http://satijalab.org/seurat
.. _Cell Ranger: https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
.. _favorably discussed: https://doi.org/10.1101/110668
.. _Monocle: http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/


Features 
---------

Let us give an Overview_ of the toplevel user functions, followed by a few words on Scanpy's `Basic Features`_ and more `details <Visualizations_>`__.

Overview
~~~~~~~~

Scanpy user functions are grouped into the following modules

sc.tools_
  Machine Learning and statistics tools. Abbreviation ``sc.tl``.
sc.preprocessing_
  Preprocessing. Abbreviation ``sc.pp``.
sc.plotting_
  Plotting. Abbreviation ``sc.pl``.
sc.settings_
  Settings.

.. _sc.tools: scanpy/tools
.. _sc.preprocessing: scanpy/preprocessing
.. _sc.plotting: scanpy/plotting
.. _sc.settings: scanpy/settings.py

Preprocessing
^^^^^^^^^^^^^

`pp.* <sc.preprocessing_>`__
  Filtering of highly-variable genes, batch-effect correction, per-cell (UMI normalization).

Visualizations
^^^^^^^^^^^^^^

`tl.pca <pca_>`__
  PCA ([pedregosa11]_).
`tl.diffmap <diffmap_>`__
  Diffusion Maps ([coifman05]_; [haghverdi15]_; [wolf17]_).
`tl.tsne <tsne_>`__
  t-SNE ([maaten08]_; [amir13]_; [pedregosa11]_).
`tl.spring <spring_>`__
  `Force-directed graph drawing`_ ([weinreb16]_).

.. _Force-directed graph drawing: https://en.wikipedia.org/wiki/Force-directed_graph_drawing

Branching trajectories and pseudotime, clustering, differential expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`tl.dpt <dpt_>`__
  Infer progression of cells, identify *branching* subgroups ([haghverdi16]_; [wolf17]_).
`tl.dbscan <dbscan_>`__
  Cluster cells into subgroups ([ester96]_; [pedregosa11]_).
`tl.diffrank <diffrank_>`__
  Rank genes according to differential expression ([wolf17]_).

Simulation
^^^^^^^^^^

`tl.sim <sim_>`__
  Simulate dynamic gene expression data ([wittmann09]_; [wolf17]_).

Basic Features
~~~~~~~~~~~~~~

The typical workflow consists of subsequent calls of data analysis tools
of the form::

    sc.tl.tool(adata, **params)

where ``adata`` is an ``AnnData`` object and ``params`` is a dictionary that stores optional parameters. Each of these calls adds annotation to an expression matrix *X*, which stores *n* *d*-dimensional gene expression measurements. By default, Scanpy tools operate *inplace* and return ``None``. If you want to copy the ``AnnData`` object, pass the ``copy`` argument::

    adata_copy = sc.tl.tool(adata, copy=True, **params)

Reading and writing data files and AnnData objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One usually calls::

    adata = sc.read(filename)

to initialize an AnnData object, possibly adds further annotation, e.g. by::

    annotation = np.genfromtxt(filename_annotation)
    adata.smp['cell_groups'] = annotation[:, 2]  # categorical annotation of type str
    adata.smp['time'] = annotation[:, 3]         # numerical annotation of type float

and uses::

    sc.write(filename, adata)

to save the ``adata`` to a file. Reading foresees filenames with extensions *h5*, *xlsx*, *mtx*, *txt*, *csv* and others. Writing foresees writing *h5*, *csv* and *txt*. Instead of providing a filename, you can provide a *filekey*, i.e., any string that does *not* end on a valid file extension. By default, Scanpy writes to ``./write/filekey.h5``, an *hdf5* file which is configurable by setting ``sc.settings.writedir`` and ``sc.settings.file_format_data``.

AnnData objects
^^^^^^^^^^^^^^^

An ``AnnData`` instance stores an array-like data matrix as ``adata.X``, dict-like sample annotation as ``adata.smp``, dict-like variable annotation as ``adata.var`` and additional unstructured dict-like annotation as ``adata.add``. While ``adata.add`` is a conventional dictionary, ``adata.smp`` and ``adata.var`` are instances of a low-level Pandas dataframe-like class.

Values can be retrieved and appended via ``adata.smp[key]`` and ``adata.var[key]``. Sample and variable names can be accessed via ``adata.smp_names`` and ``adata.var_names``, respectively. AnnData objects can be sliced like Pandas dataframes, for example, ``adata = adata[:, list_of_gene_names]``. The AnnData class is similar to R's ExpressionSet ([huber15]_); the latter though is not implemented for sparse data.

Plotting
^^^^^^^^

For each tool, there is an associated plotting function::

    sc.pl.tool(adata)

that retrieves and plots the elements of ``adata`` that were previously written by ``sc.tl.tool(adata)``. To not display figures interactively but save all plots to default locations, you can set ``sc.sett.savefigs = True``. By default, figures are saved as *png* to ``./figs/``. Reset ``sc.sett.file_format_figs`` and ``sc.sett.figdir`` if you want to change this. Scanpy's plotting module can be viewed similar to Seaborn_: an extension of matplotlib_ that allows visualizing certain frequent tasks with one-line commands. Detailed configuration has to be done via matplotlib functions, which is easy as Scanpy's plotting functions usually return a ``Matplotlib.Axes`` object.

.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/

Builtin examples
^^^^^^^^^^^^^^^^

Show all builtin example data using ``sc.show_exdata()`` and all builtin example use cases via ``sc.show_examples()``. Load annotated and preprocessed data using an *example key*, here 'paul15', via::

    adata = sc.get_example('paul15')

The key 'paul15' can also be used within ``sc.read('paul15')`` and ``sc.write('paul15')`` to write the current state of the AnnData object to disk.

Visualization
~~~~~~~~~~~~~

pca
^^^

`[source] <scanpy/tools/pca.py>`__ Computes the PCA representation ``X_pca`` of data, principal components and variance decomposition. Uses the implementation of the ``scikit-learn`` package ([pedregosa11]_).

tsne
^^^^

`[source] <scanpy/tools/tsne.py>`__ Computes the tSNE representation ``X_tsne`` of data.

The algorithm has been introduced by [maaten08]_ and proposed for single-cell data by [amir13]_. By default, Scanpy uses the implementation of the ``scikit-learn`` package ([pedregosa11]_). You can achieve a huge speedup if you install the Multicore-tSNE package by [ulyanov16]_, which will be automatically detected by Scanpy.

diffmap
^^^^^^^

`[source] <scanpy/tools/diffmap.py>`__ Computes the diffusion maps representation ``X_diffmap`` of data.

Diffusion maps ([coifman05]_) has been proposed for visualizing single-cell data by [haghverdi15]_. The tool uses the adapted Gaussian kernel suggested by [haghverdi16]_. The Scanpy implementation is due to [wolf17]_.

spring
^^^^^^

Beta version.

`[source] <scanpy/tools/spring.py>`__ Force-directed graph drawing is a long-established algorithm for visualizing graphs, see `Force-directed graph drawing`_. It has been suggested for visualizing single-cell data by [weinreb16]_.

Here, the Fruchterman & Reingold ([fruchterman91]_) algorithm is used. The implementation uses elements of the NetworkX implementation ([hagberg08]_).

Discrete clustering of subgroups and continuous progression through subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dpt
^^^

`[source] <scanpy/tools/dpt.py>`__ Reconstruct the progression of a biological process from snapshot data and detect branching subgroups. Diffusion Pseudotime analysis has been introduced by [haghverdi16]_ and implemented for Scanpy by [wolf17]_.

The functionality of diffmap and dpt compare to the R package destiny_ of [angerer16]_, but run faster and scale to much higher cell numbers.

*Examples:* See one of the early examples [notebook_, `command line`_] dealing with data of [moignard15]_.


.. _destiny: http://bioconductor.org/packages/destiny
.. _notebook: https://github.com/theislab/scanpy_usage/tree/master/170503_moignard15.ipynb
.. _command line: https://github.com/theislab/scanpy_usage/tree/master/EXAMPLES.md#moignard15

dbscan
^^^^^^

`[source] <scanpy/tools/dbscan.py>`__ Cluster cells using DBSCAN_ ([ester96]_), in the implementation of ``scikit-learn`` ([pedregosa11]_).

This is a very simple clustering method. A better one -- in the same framework as DPT and Diffusion Maps -- will come soon.

.. _DBSCAN: https://en.wikipedia.org/wiki/DBSCAN

Differential expression
~~~~~~~~~~~~~~~~~~~~~~~

diffrank
^^^^^^^^

`[source] <scanpy/tools/diffrank.py>`__ Rank genes by differential expression.

Simulation
~~~~~~~~~~

sim
^^^

`[source] <scanpy/tools/sim.py>`__ Sample from a stochastic differential equation model built from literature-curated boolean gene regulatory networks, as suggested by [wittmann09]_. The Scanpy implementation is due to [wolf17]_.

The tool compares to the Matlab tool *Odefy* of [krumsiek10]_.

Installation 
-------------

If you use Windows or Mac OS X and do not have a current Python distribution (Python 3.5 or 3.6), download and install Miniconda_ (see below). If you use Linux, use your package manager to obtain a current python distribution.

Then, download or clone the repository -- green button on top of the page -- and ``cd`` into its root directory. To install with symbolic links (stay up to date with your cloned version after you update with ``git pull``) call::

    pip install -e .

and work with the top-level command ``scanpy`` or::

    import scanpy as sc

in any directory.

Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions. Either reopen a new terminal or ``source ~/.bashrc`` on Linux/ ``source ~/.bash_profile`` on Mac. The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html

PyPi
~~~~

The package is registered_ in the `Python Packaging Index`_, but
versioning has not started yet. In the future, installation will also be
possible without reference to GitHub via ``pip install scanpy``.

.. _registered: https://pypi.python.org/pypi/scanpy
.. _Python Packaging Index: https://pypi.python.org/pypi

References
----------

.. [amir13] Amir *et al.* (2013),
   *viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia*
   `Nature Biotechnology 31, 545 <http://dx.doi.org/10.1038/nbt.2594>`__.

.. [angerer16] Angerer *et al.* (2016),
   *destiny -- diffusion maps for large-scale single-cell data in R*,
   `Bioinformatics 32, 1241-1243 <https://doi.org/10.1093/bioinformatics/btv715>`__.

.. [coifman05] Coifman *et al.* (2005),
   *Geometric diffusions as a tool for harmonic analysis and structure definition of data: Diffusion maps*,
   `PNAS 102, 7426 <http://dx.doi.org/10.1038/nmeth.3971>`__.

.. [ester96] Ester *et al.* (1996),
   *A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise*
   `Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining,
   Portland, OR, pp. 226-231 <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.121.9220>`__.

.. [fruchterman91] Fruchterman & Reingold (1991),
   *Graph drawing by force-directed placement*
   `Software: Practice & Experience <http://doi.org:10.1002/spe.4380211102>`__

.. [hagberg08] Hagberg *et al.* (2008),
   *Exploring Network Structure, Dynamics, and Function using NetworkX*
   `Scipy Conference <http://conference.scipy.org/proceedings/SciPy2008/paper_2/>`__

.. [haghverdi15] Haghverdi *et al.* (2015),
   *Diffusion maps for high-dimensional single-cell analysis of differentiation data*,
   `Bioinformatics 31, 2989 <http://dx.doi.org/10.1093/bioinformatics/btv325>`__.

.. [haghverdi16] Haghverdi *et al.* (2016),
   *Diffusion pseudotime robustly reconstructs branching cellular lineages*,
   `Nature Methods 13, 845 <http://dx.doi.org/10.1038/nmeth.3971>`__.

.. [huber15] Huber *et al.* (2015),
   *Orchestrating high-throughput genomic analysis with Bioconductor*
   `Nature Methods <https://doi.org/10.1038/nmeth.3252>`__

.. [krumsiek10] Krumsiek *et al.* (2010),
   *Odefy -- From discrete to continuous models*,
   `BMC Bioinformatics 11, 233 <http://dx.doi.org/10.1186/1471-2105-11-233>`__.

.. [krumsiek11] Krumsiek *et al.* (2011),
   *Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*,
   `PLoS ONE 6, e22649 <http://dx.doi.org/10.1371/journal.pone.0022649>`__.

.. [maaten08] Maaten & Hinton (2008),
   *Visualizing data using t-SNE*,
   `JMLR 9, 2579 <http://www.jmlr.org/papers/v9/vandermaaten08a.html>`__.

.. [macosco15] Macosko *et al.* (2015)
   *Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets*
   `Cell <https://doi.org/10.1016/j.cell.2015.05.002>`__

.. [moignard15] Moignard *et al.* (2015),
   *Decoding the regulatory network of early blood development from single-cell gene expression measurements*,
   `Nature Biotechnology 33, 269 <http://dx.doi.org/10.1038/nbt.3154>`__.

.. [pedregosa11] Pedregosa *et al.* (2011),
   *Scikit-learn: Machine Learning in Python*,
   `JMLR 12, 2825 <http://www.jmlr.org/papers/v12/pedregosa11a.html>`__.

.. [paul15] Paul *et al.* (2015),
   *Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*,
   `Cell 163, 1663 <http://dx.doi.org/10.1016/j.cell.2015.11.013>`__.

.. [ulyanov16] Ulyanov (2016),
   *Multicore t-SNE*
   `GitHub <https://github.com/DmitryUlyanov/Multicore-TSNE>`__

.. [weinreb16] Weinreb *et al.* (2016),
   *SPRING: a kinetic interface for visualizing high dimensional single-cell expression data*
   `bioRXiv <https://doi.org/10.1101/090332>`__

.. [wittmann09] Wittmann *et al.* (2009),
   *Transforming Boolean models to continuous models: methodology and application to T-cell receptor signaling*,
   `BMC Systems Biology 3, 98 <http://dx.doi.org/10.1186/1752-0509-3-98>`__.

.. [wolf17] Wolf *et al* (2017),
   TBD

.. [zheng] Zheng *et al.* (2017),
   *Massively parallel digital transcriptional profiling of single cells*
   `Nature Communications <https://doi.org/10.1038/ncomms14049>`__
