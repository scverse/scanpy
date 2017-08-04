`Getting started`_ \| Examples_ \| Docs_ \| Installation_ \| References_

|Build Status|

.. |Build Status| image:: https://travis-ci.org/theislab/scanpy.svg?branch=master
   :target: https://travis-ci.org/theislab/scanpy

Scanpy – Single-Cell Analysis in Python
=======================================

Scanpy is a scalable toolkit for analyzing single-cell gene expression data. It includes preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing and simulation of gene regulatory networks. The Python-based implementation efficiently deals with data sets of more than one million cells and enables easy integration of advanced machine learning algorithms.

For conceptual ideas and context, see our `draft <http://falexwolf.de/docs/scanpy.pdf>`__; comments are highly appreciated.

Getting started
---------------

With Python 3.5 or 3.6 installed, get `releases on PyPI <https://pypi.python.org/pypi/scanpy>`__ via (more information on installation `here <Installation_>`__)::

  pip3 install scanpy

To work with the latest version on `GitHub <https://github.com/theislab/scanpy>`__: clone the repository – green button on top of the page – and ``cd`` into its root directory and type::

    pip3 install -e .

You can now ``import scanpy.api as sc`` anywhere on your system and work with the command ``scanpy`` on the command-line.

Examples
--------

Examples are collected in the repo scanpy_usage_. Good starting points are the following use cases:

.. _scanpy_usage: https://github.com/theislab/scanpy_usage

17-05-05_
  We reproduce most of the `Guided Clustering tutorial`_ of Seurat_ [Satija15]_.
17-05-03_
  Analyzing 68 000 cells from [Zheng17]_, we find that Scanpy is about a factor 5 to 16 faster and more memory efficient than the `Cell Ranger`_ R kit for secondary analysis.
17-05-02_
  We reproduce the results of the Diffusion Pseudotime (DPT) paper of [Haghverdi16]_. Note that DPT has recently been very `favorably discussed`_ by the authors of Monocle_.

.. _17-05-05: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
.. _17-05-03: https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17
.. _17-05-02: https://github.com/theislab/scanpy_usage/tree/master/170502_haghverdi16
.. _17-04-30: https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11

.. _Guided Clustering tutorial: http://satijalab.org/seurat/pbmc3k_tutorial.html
.. _Seurat: http://satijalab.org/seurat
.. _Cell Ranger: https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
.. _favorably discussed: https://doi.org/10.1101/110668
.. _Monocle: http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/


Docs
----

Here, we give an Overview_ of the toplevel user functions, describe `Basic Features`_ and the context of the `Tools <Visualization_>`__. For detailed help on the functions, use Python's ``help``. A separate docs page will soon be established.

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

.. _sc.tools:         https://github.com/theislab/scanpy/tree/master/scanpy/tools
.. _sc.preprocessing: https://github.com/theislab/scanpy/tree/master/scanpy/preprocessing
.. _sc.plotting:      https://github.com/theislab/scanpy/tree/master/scanpy/plotting
.. _sc.settings:      https://github.com/theislab/scanpy/tree/master/scanpy/settings.py

Preprocessing
^^^^^^^^^^^^^

`pp.* <sc.preprocessing_>`__
  Filtering of highly-variable genes, batch-effect correction, per-cell (UMI) normalization, preprocessing recipes.

Visualizations
^^^^^^^^^^^^^^

`tl.pca <pca_>`__
  PCA [Pedregosa11]_.
`tl.diffmap <diffmap_>`__
  Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf17]_.
`tl.tsne <tsne_>`__
  t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.
`tl.draw_graph <draw_graph_>`__
  Force-directed graph drawing [Fruchterman91]_ [Weinreb17]_ [Csardi06]_.

Branching trajectories and pseudotime, clustering, differential expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`tl.dpt <dpt_>`__
  Infer progression of cells, identify *branching* subgroups [Haghverdi16]_ [Wolf17]_.
`tl.louvain <louvain_>`__
  Cluster cells into subgroups [Blondel08]_ [Levine15]_ [Traag17]_.
`tl.rank_genes_groups <rank_genes_groups_>`__
  Rank genes according to differential expression [Wolf17]_.

Simulations
^^^^^^^^^^^

`tl.sim <sim_>`__
  Simulate dynamic gene expression data [Wittmann09]_ [Wolf17]_.

Basic Features
~~~~~~~~~~~~~~

The typical workflow consists of subsequent calls of data analysis tools
of the form::

    sc.tl.tool(adata, **params)

where ``adata`` is an ``AnnData`` object and ``params`` are optional parameters. Each of these calls adds annotation to an expression matrix *X*, which stores *n* *d*-dimensional gene expression measurements. To facilitate writing memory-efficient pipelines, by default, Scanpy tools operate *inplace* on ``adata`` and return ``None``. If you want to copy the ``AnnData`` object, pass the ``copy`` argument::

    adata_copy = sc.tl.tool(adata, copy=True, **params)

Reading and writing data files and AnnData objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One usually calls::

    adata = sc.read(filename)

to initialize an AnnData object, possibly adds further annotation using, e.g., ``np.genfromtxt`` or ``pd.read_csv``::

    annotation = pd.read_csv(filename_annotation)
    adata.smp['cell_groups'] = annotation['cell_groups']  # categorical annotation of type str or int
    adata.smp['time'] = annotation['time']                # numerical annotation of type float

and uses::

    sc.write(filename, adata)

to save the ``adata`` as a collection of data arrays to a file in a platform and language-independent way. Reading foresees filenames with extensions *h5*, *xlsx*, *mtx*, *txt*, *csv* and others. Writing foresees writing *h5*, *csv* and *txt*. Instead of providing a filename, you can provide a *filekey*, i.e., any string that does *not* end on a valid file extension.

AnnData objects
^^^^^^^^^^^^^^^

An ``AnnData`` instance stores an array-like data matrix as ``adata.X``, dict-like sample annotation as ``adata.smp``, dict-like variable annotation as ``adata.var`` and additional unstructured dict-like annotation as ``adata.add``. While ``adata.add`` is a conventional dictionary, ``adata.smp`` and ``adata.var`` are instances of a low-level Pandas dataframe-like class.

Values can be retrieved and appended via ``adata.smp[key]`` and ``adata.var[key]``. Sample and variable names can be accessed via ``adata.smp_names`` and ``adata.var_names``, respectively. AnnData objects can be sliced like Pandas dataframes, for example, ``adata = adata[:, list_of_gene_names]``. The AnnData class is similar to R's ExpressionSet [Huber15]_ the latter though is not implemented for sparse data.

Plotting
^^^^^^^^

For each tool, there is an associated plotting function::

    sc.pl.tool(adata)

that retrieves and plots annotation in ``adata`` that has been added by ``sc.tl.tool(adata)``. Scanpy's plotting module can be viewed similar to Seaborn_: an extension of matplotlib_ that allows visualizing operations on AnnData objects with one-line commands. Detailed configuration has to be done via matplotlib functions, which is easy as Scanpy's plotting functions accept and return a ``Matplotlib.Axes`` object.

.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/


Visualization
~~~~~~~~~~~~~

pca
^^^

`[source] <tl.pca_>`__ Computes PCA coordinates, loadings and variance decomposition. Uses the implementation of *scikit-learn* [Pedregosa11]_.

tsne
^^^^

`[source] <tl.tsne_>`__ t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been proposed for single-cell data by [Amir13]_. By default, Scanpy uses the implementation of *scikit-learn* [Pedregosa11]_. You can achieve a huge speedup if you install *Multicore-tSNE* by [Ulyanov16]_, which will be automatically detected by Scanpy.

diffmap
^^^^^^^

`[source] <tl.diffmap_>`__ Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested by [Haghverdi16]_. Uses the implementation of [Wolf17]_.

draw_graph
^^^^^^^^^^

`[source] <tl.draw_graph_>`__ `Force-directed graph drawing`_ describes a class of long-established algorithms for visualizing graphs. It has been suggested for visualizing single-cell data by [Weinreb17]_. Here, by default, the Fruchterman & Reingold [Fruchterman91]_ algorithm is used; many other layouts are available. Uses the igraph implementation [Csardi06]_.

.. _Force-directed graph drawing: https://en.wikipedia.org/wiki/Force-directed_graph_drawing

Discrete clustering of subgroups, continuous progression through subgroups, differential expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dpt
^^^

`[source] <tl.dpt_>`__ Reconstruct the progression of a biological process from snapshot data and detect branching subgroups. Diffusion Pseudotime analysis has been introduced by [Haghverdi16]_. Here, we use a further developed version, which is able to detect multiple branching events [Wolf17]_.

The possibilities of *diffmap* and *dpt* are similar to those of the R package destiny_ of [Angerer16]_. The Scanpy tools though run faster and scale to much higher cell numbers.

*Examples:* See this `use case <17-05-02_>`__.

.. _destiny: http://bioconductor.org/packages/destiny

louvain
^^^^^^^

`[source] <tl.louvain_>`__ Cluster cells using the Louvain algorithm [Blondel08]_ in the implementation of [Traag17]_. The Louvain algorithm has been proposed for single-cell analysis by [Levine15]_.

*Examples:* See this `use case <17-05-05_>`__.

rank_genes_groups
^^^^^^^^^^^^^^^^^

`[source] <tl.rank_genes_groups_>`__ Rank genes by differential expression.

*Examples:* See this `use case <17-05-05_>`__.


Simulation
~~~~~~~~~~

sim
^^^

`[source] <scanpy/tools/sim.py>`__ Sample from a stochastic differential equation model built from literature-curated boolean gene regulatory networks, as suggested by [Wittmann09]_. The Scanpy implementation is due to [Wolf17]_.

The tool is similar to the Matlab tool *Odefy* of [Krumsiek10]_.

*Examples:* See this `use case <17-04-30_>`__.

.. _tl.pca:               https://github.com/theislab/scanpy/tree/master/scanpy/tools/pca.py
.. _tl.tsne:              https://github.com/theislab/scanpy/tree/master/scanpy/tools/tsne.py
.. _tl.diffmap:           https://github.com/theislab/scanpy/tree/master/scanpy/tools/diffmap.py
.. _tl.draw_graph:        https://github.com/theislab/scanpy/tree/master/scanpy/tools/draw_graph.py
.. _tl.dpt:               https://github.com/theislab/scanpy/tree/master/scanpy/tools/dpt.py
.. _tl.louvain:           https://github.com/theislab/scanpy/tree/master/scanpy/tools/louvain.py
.. _tl.rank_genes_groups: https://github.com/theislab/scanpy/tree/master/scanpy/tools/rank_genes_groups.py


Installation
------------

If you use Windows or Mac OS X and do not have Python 3.5 or 3.6, download and install Miniconda_ (see below). If you use Linux, use your package manager to obtain a current Python distribution.

Get `releases on PyPI <https://pypi.python.org/pypi/scanpy>`__ via::

  pip3 install scanpy

To work with the latest version on `GitHub <https://github.com/theislab/scanpy>`__: clone the repository – green button on top of the page – and ``cd`` into its root directory. To install with symbolic links (stay up to date with your cloned version after you update with ``git pull``) call::

    pip3 install -e .

You can now ``import scanpy.api as sc`` anywhere on your system and work with the command ``scanpy`` on the command-line.


Installing Miniconda
~~~~~~~~~~~~~~~~~~~~

After downloading Miniconda_, in a unix shell (Linux, Mac), run

.. code:: shell

    cd DOWNLOAD_DIR
    chmod +x Miniconda3-latest-VERSION.sh
    ./Miniconda3-latest-VERSION.sh

and accept all suggestions. Either reopen a new terminal or ``source ~/.bashrc`` on Linux/ ``source ~/.bash_profile`` on Mac. The whole process takes just a couple of minutes.

.. _Miniconda: http://conda.pydata.org/miniconda.html

Trouble shooting
~~~~~~~~~~~~~~~~

If you do not have sudo rights (you get a ``Permission denied`` error)::

    pip install --user scanpy

On MacOS, you probably need to install the C core of igraph via homebrew first

- ``brew install igraph``
- If python-igraph still fails to install, see `here <https://stackoverflow.com/questions/29589696/problems-compiling-c-core-of-igraph-with-python-2-7-9-anaconda-2-2-0-on-mac-osx>`__ or consider installing gcc via ``brew install gcc --without-multilib`` and exporting ``export CC="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"; export CXX="/usr/local/Cellar/gcc/X.x.x/bin/gcc-X"``, where ``X`` and ``x`` refers to the version of ``gcc``; in my case, the path reads ``/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6``.


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

.. [Satija15] Satija *et al.* (2015),
   *Spatial reconstruction of single-cell gene expression data*,
   `Nature Biotechnology <https://doi.org/10.1038/nbt.3192>`__.

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
