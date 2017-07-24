`Getting started`_ \| Features_ \| Installation_ \| References_

|Build Status|

.. |Build Status| image:: https://travis-ci.org/theislab/scanpy.svg?branch=master
   :target: https://travis-ci.org/theislab/scanpy

Scanpy – Single-Cell Analysis in Python
=======================================

Efficient tools for analyzing and simulating large-scale single-cell data that aim at an understanding of dynamic biological processes from snapshots of transcriptome or proteome. The draft `Wolf, Angerer & Theis (2017) <http://falexwolf.de/docs/scanpy.pdf>`__ explains conceptual ideas of the package. Any comments are appreciated!

Getting started
---------------

Download or clone the repository – green button on top of the page – and ``cd`` into its root directory. With Python 3.5 or 3.6 (preferably Miniconda_) installed, type::

    pip install -e .

Aside from enabling ``import scanpy as sc`` anywhere on your system, you can also work with the top-level command ``scanpy`` on the command-line (more info `here <Installation_>`__).

Then go through the use cases compiled in scanpy_usage_, in particular, the recent additions

.. _scanpy_usage: https://github.com/theislab/scanpy_usage

17-05-05_
  We reproduce most of the `Guided Clustering tutorial`_ of Seurat_ [Macosco15]_.
17-05-03_
  Analyzing 68 000 cells from [Zheng17]_, we find that Scanpy is about a factor 5 to 16 faster and more memory efficient than the `Cell Ranger`_ R kit for secondary analysis.
17-05-02_
  We reproduce the results of the Diffusion Pseudotime (DPT) paper of [Haghverdi16]_. Note that DPT has recently been very `favorably discussed`_ by the authors of Monocle_.

.. _17-05-05: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
.. _17-05-03: https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17
.. _17-05-02: https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16

.. _Guided Clustering tutorial: http://satijalab.org/seurat/pbmc-tutorial.html
.. _Seurat: http://satijalab.org/seurat
.. _Cell Ranger: https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
.. _favorably discussed: https://doi.org/10.1101/110668
.. _Monocle: http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/


Features 
---------

Let us give an Overview_ of the toplevel user functions, followed by a few words on Scanpy's `Basic Features`_ and more `details <Visualization_>`__.

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
  Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization.

Visualizations
^^^^^^^^^^^^^^

`tl.pca <pca_>`__
  PCA [Pedregosa11]_.
`tl.diffmap <diffmap_>`__
  Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf17]_.
`tl.tsne <tsne_>`__
  t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.
`tl.draw_graph <draw_graph_>`__
  Force-directed graph drawing [Csardi06]_ [Weinreb17]_.

Branching trajectories and pseudotime, clustering, differential expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`tl.dpt <dpt_>`__
  Infer progression of cells, identify *branching* subgroups [Haghverdi16]_ [Wolf17]_.
`tl.louvain <louvain_>`__
  Cluster cells into subgroups [Blondel08]_ [Traag17]_.
`tl.rank_genes_groups <rank_genes_groups_>`__
  Rank genes according to differential expression [Wolf17]_.

Simulation
^^^^^^^^^^

`tl.sim <sim_>`__
  Simulate dynamic gene expression data [Wittmann09]_ [Wolf17]_.

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

to save the ``adata`` to a file. Reading foresees filenames with extensions *h5*, *xlsx*, *mtx*, *txt*, *csv* and others. Writing foresees writing *h5*, *csv* and *txt*. Instead of providing a filename, you can provide a *filekey*, i.e., any string that does *not* end on a valid file extension.

AnnData objects
^^^^^^^^^^^^^^^

An ``AnnData`` instance stores an array-like data matrix as ``adata.X``, dict-like sample annotation as ``adata.smp``, dict-like variable annotation as ``adata.var`` and additional unstructured dict-like annotation as ``adata.add``. While ``adata.add`` is a conventional dictionary, ``adata.smp`` and ``adata.var`` are instances of a low-level Pandas dataframe-like class.

Values can be retrieved and appended via ``adata.smp[key]`` and ``adata.var[key]``. Sample and variable names can be accessed via ``adata.smp_names`` and ``adata.var_names``, respectively. AnnData objects can be sliced like Pandas dataframes, for example, ``adata = adata[:, list_of_gene_names]``. The AnnData class is similar to R's ExpressionSet [Huber15]_ the latter though is not implemented for sparse data.

Plotting
^^^^^^^^

For each tool, there is an associated plotting function::

    sc.pl.tool(adata)

that retrieves and plots the elements of ``adata`` that were previously written by ``sc.tl.tool(adata)``. Scanpy's plotting module can be viewed similar to Seaborn_: an extension of matplotlib_ that allows visualizing operations on AnnData objects with one-line commands. Detailed configuration has to be done via matplotlib functions, which is easy as Scanpy's plotting functions accept and return a ``Matplotlib.Axes`` object.

.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/


Visualization
~~~~~~~~~~~~~

pca
^^^

`[source] <scanpy/tools/pca.py>`__ Computes the PCA representation ``X_pca`` of data, principal components and variance decomposition. Uses the implementation of the ``scikit-learn`` package [Pedregosa11]_.

tsne
^^^^

`[source] <scanpy/tools/tsne.py>`__ Computes the tSNE representation ``X_tsne`` of data.

The algorithm has been introduced by [Maaten08]_ and proposed for single-cell data by [Amir13]_. By default, Scanpy uses the implementation of the ``scikit-learn`` package [Pedregosa11]_. You can achieve a huge speedup if you install the Multicore-tSNE package by [Ulyanov16]_, which will be automatically detected by Scanpy.

diffmap
^^^^^^^

`[source] <scanpy/tools/diffmap.py>`__ Computes the diffusion maps representation ``X_diffmap`` of data.

Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested by [Haghverdi16]_. The Scanpy implementation is due to [Wolf17]_.

draw_graph
^^^^^^^^^^

`[source] <scanpy/tools/draw_graph.py>`__ Force-directed graph drawing is a long-established algorithm for visualizing graphs, see `Force-directed graph drawing`_. It has been suggested for visualizing single-cell data by [Weinreb17]_.

Here, the Fruchterman & Reingold [Fruchterman91]_ algorithm is used by default, but many other layouts are available. We use the igraph implementation [Csardi06]_.

Discrete clustering of subgroups, continuous progression through subgroups, differential expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dpt
^^^

`[source] <scanpy/tools/dpt.py>`__ Reconstruct the progression of a biological process from snapshot data and detect branching subgroups. Diffusion Pseudotime analysis has been introduced by [Haghverdi16]_ and implemented for Scanpy by [Wolf17]_.

The functionality of diffmap and dpt compare to the R package destiny_ of [Angerer16]_, but run faster and scale to much higher cell numbers.

*Examples:* See this `use case <https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16>`__.

.. _destiny: http://bioconductor.org/packages/destiny

louvain
^^^^^^

`[source] <scanpy/tools/louvain.py>`__ Cluster cells using the Louvain algorithm [Blondel08]_ in the implementation of [Traag17]_.

The Louvain algorithm has been proposed for single-cell analysis by [Levine15]_.

*Examples:* See this `use case <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`__.


rank_genes_groups
^^^^^^^^^^^^^^^^^

`[source] <scanpy/tools/rank_genes_groups.py>`__ Rank genes by differential expression.

*Examples:* See this `use case <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`__.


Simulation
~~~~~~~~~~

sim
^^^

`[source] <scanpy/tools/sim.py>`__ Sample from a stochastic differential equation model built from literature-curated boolean gene regulatory networks, as suggested by [Wittmann09]_. The Scanpy implementation is due to [Wolf17]_.

The tool compares to the Matlab tool *Odefy* of [Krumsiek10]_.

*Examples:* See this `use case <https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11>`__.


Installation 
-------------

If you use Windows or Mac OS X and do not have a current Python distribution (Python 3.5 or 3.6), download and install Miniconda_ (see below). If you use Linux, use your package manager to obtain a current python distribution.

Then, download or clone the repository – green button on top of the page – and ``cd`` into its root directory. To install with symbolic links (stay up to date with your cloned version after you update with ``git pull``) call::

    pip install -e .

and work with the top-level command ``scanpy`` or::

    import scanpy.api as sc

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

.. [Amir13] Amir *et al.* (2013),
   *viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia*,
   `Nature Biotechnology <https://doi.org/10.1038/nbt.2594>`__.

.. [Angerer16] Angerer *et al.* (2016),
   *destiny – diffusion maps for large-scale single-cell data in R*,
   `Bioinformatics <https://doi.org/10.1093/bioinformatics/btv715>`__.

.. [Blondel08] Blondel *et al.* (2008),
   *Fast unfolding of communities in large networks*,
   `J. Stat. Mech. <https://doi.org/10.1088/1742-5468/2008/10/P10008>`__.   

.. [Coifman05] Coifman *et al.* (2005),
   *Geometric diffusions as a tool for harmonic analysis and structure definition of data: Diffusion maps*,
   `PNAS <https://doi.org/10.1038/nmeth.3971>`__.

.. [Csardi06] Csardi *et al.* (2006),
   *The igraph software package for complex network researc*,
   `InterJournal Complex Systems <http://igraph.org>`__.

   
.. [Ester96] Ester *et al.* (1996),
   *A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise*,
   `Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining,
   Portland, OR <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.121.9220>`__.

.. [Fruchterman91] Fruchterman & Reingold (1991),
   *Graph drawing by force-directed placement*,
   `Software: Practice & Experience <http://doi.org:10.1002/spe.4380211102>`__.

.. [Hagberg08] Hagberg *et al.* (2008),
   *Exploring Network Structure, Dynamics, and Function using NetworkX*,
   `Scipy Conference <http://conference.scipy.org/proceedings/SciPy2008/paper_2/>`__.

.. [Haghverdi15] Haghverdi *et al.* (2015),
   *Diffusion maps for high-dimensional single-cell analysis of differentiation data*,
   `Bioinformatics <https://doi.org/10.1093/bioinformatics/btv325>`__.

.. [Haghverdi16] Haghverdi *et al.* (2016),
   *Diffusion pseudotime robustly reconstructs branching cellular lineages*,
   `Nature Methods <https://doi.org/10.1038/nmeth.3971>`__.

.. [Huber15] Huber *et al.* (2015),
   *Orchestrating high-throughput genomic analysis with Bioconductor*,
   `Nature Methods <https://doi.org/10.1038/nmeth.3252>`__.

.. [Krumsiek10] Krumsiek *et al.* (2010),
   *Odefy – From discrete to continuous models*,
   `BMC Bioinformatics <https://doi.org/10.1186/1471-2105-11-233>`__.

.. [Krumsiek11] Krumsiek *et al.* (2011),
   *Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*,
   `PLoS ONE <https://doi.org/10.1371/journal.pone.0022649>`__.

.. [Levine15] Levine *et al.* (2015),
   *Data-Driven Phenotypic Dissection of AML Reveals Progenitor--like Cells that Correlate with Prognosis*,
   `Cell <https://doi.org/10.1016/j.cell.2015.05.047>`__.
   
.. [Maaten08] Maaten & Hinton (2008),
   *Visualizing data using t-SNE*,
   `JMLR <http://www.jmlr.org/papers/v9/vandermaaten08a.html>`__.

.. [Macosco15] Macosko *et al.* (2015),
   *Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets*,
   `Cell <https://doi.org/10.1016/j.cell.2015.05.002>`__.

.. [Moignard15] Moignard *et al.* (2015),
   *Decoding the regulatory network of early blood development from single-cell gene expression measurements*,
   `Nature Biotechnology <https://doi.org/10.1038/nbt.3154>`__.

.. [Pedregosa11] Pedregosa *et al.* (2011),
   *Scikit-learn: Machine Learning in Python*,
   `JMLR <http://www.jmlr.org/papers/v12/pedregosa11a.html>`__.

.. [Paul15] Paul *et al.* (2015),
   *Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*,
   `Cell <https://doi.org/10.1016/j.cell.2015.11.013>`__.

.. [Traag17] Traag (2017),
   *Louvain*,
   `GitHub <https://doi.org/10.5281/zenodo.35117>`__.
   
.. [Ulyanov16] Ulyanov (2016),
   *Multicore t-SNE*,
   `GitHub <https://github.com/DmitryUlyanov/Multicore-TSNE>`__.

.. [Weinreb17] Weinreb *et al.* (2016),
   *SPRING: a kinetic interface for visualizing high dimensional single-cell expression data*,
   `bioRXiv <https://doi.org/10.1101/090332>`__.

.. [Wittmann09] Wittmann *et al.* (2009),
   *Transforming Boolean models to continuous models: methodology and application to T-cell receptor signaling*,
   `BMC Systems Biology <https://doi.org/10.1186/1752-0509-3-98>`__.

.. [Wolf17] Wolf *et al* (2017),
   TBD.

.. [Zheng17] Zheng *et al.* (2017),
   *Massively parallel digital transcriptional profiling of single cells*,
   `Nature Communications <https://doi.org/10.1038/ncomms14049>`__.
