[Getting Started](#getting_started) |
[Features](#features) |
[Installation](#install) |
[References](#references)

[![Build Status](https://travis-ci.org/theislab/scanpy.svg?branch=master)](https://travis-ci.org/theislab/scanpy)

# Scanpy - Single-Cell Analysis in Python

Efficient tools for analyzing and simulating large-scale single-cell data that aim at an understanding
of dynamic biological processes from snapshots of transcriptome or
proteome. The draft [Wolf, Angerer & Theis (2017)](http://falexwolf.de/docs/scanpy.pdf) explains conceptual ideas of the package. Any comments are appreciated!


## Getting started <a id="getting_started"></a>

Download or clone the repository - green button on top of the page - and `cd`
into its root directory. With Python 3.5 or 3.6 (preferably [Miniconda](http://conda.pydata.org/miniconda.html)) installed, type
```
pip install -e .
```
Aside from enabling `import scanpy.api as sc` anywhere on your system, you can also work
with the top-level command `scanpy` on the command-line (more info on installation [here](#install)).

Then go through the use cases compiled in
[scanpy_usage](https://github.com/theislab/scanpy_usage), in particular, the recent additions

* 17-05-22 | [link](https://github.com/theislab/scanpy_usage/tree/master/170522_visualizing_one_million_cells) | Visualizing the [publicly available](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1M_neurons) 10X 1.3 mio brain cell data set using the command line.

* 17-05-05 | [link](https://github.com/theislab/scanpy_usage/tree/master/170505_seurat) | We reproduce parts of the recent [Guided Clustering tutorial](http://satijalab.org/seurat/pbmc-tutorial.html) of [*Seurat*](http://satijalab.org/seurat/) [(Macosko *et al.*, Cell 2015)](https://doi.org/10.1016/j.cell.2015.05.002).

* 17-05-03 | [link](https://github.com/theislab/scanpy_usage/tree/master/170503_zheng17) | Analyzing 68 000 cells from [(Zheng *et al.*, Nat. Comms. 2017)](https://doi.org/10.1038/ncomms14049), we find that Scanpy is about a factor 5 to 10 faster and more memory efficient than the [*Cell Ranger*](https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis) R kit for secondary analysis. For large-scale data, this becomes crucial for an interactive analysis.

* 17-05-01 | [link](https://github.com/theislab/scanpy_usage/tree/master/170501_moignard15) | Diffusion Pseudotime analysis resolves developmental processes in data of [Moignard *et al*, Nat. Biotechn. (2015)](https://doi.org/10.1038/nbt.3154), reproducing results of [Haghverdi *et al.*, Nat. Meth. (2016)](https://doi.org/10.1038/nmeth.3971). Also, note that DPT has recently been very [favorably discussed](https://doi.org/10.1101/110668) by the authors of [Monocle](http://cole-trapnell-lab.github.io/monocle-release/articles/v2.0.0/).


## Features <a id="features"></a>

We first give an [overview](#overview) of the toplevel user functions of `scanpy.api`, followed by a few words on Scanpy's [basic features](#basic-features) and more [details](#visualization-1). For usage of the command-line interface, which is parallels usage of the API, see this introductory [example](https://github.com/theislab/scanpy_usage/blob/master/170501_moignard15).

### Overview

Scanpy user functions are grouped into the following modules

* [`sc.tools`](scanpy/tools) - Machine Learning and statistics tools. Abbreviation `sc.tl`.
* [`sc.preprocessing`](scanpy/preprocessing) - Preprocessing. Abbreviation `sc.pp`.
* [`sc.plotting`](scanpy/plotting) - Plotting. Abbreviation `sc.pl`.
* [`sc.settings`](scanpy/settings.py) - Settings. Abbreviation `sc.sett`.

#### Preprocessing

* [`pp.*`](scanpy/preprocessing) - Filtering of highly-variable genes,
batch-effect correction, per-cell (UMI) normalization.

#### Visualization

* [tl.pca](#pca) - PCA ([Pedregosa *et al.*, 2011](#ref_pedregosa11)).

* [tl.diffmap](#diffmap) - Diffusion Maps
([Coifman *et al.*, 2005](#ref_coifman05); [Haghverdi *et al.*,
2015](#ref_haghverdi15); [Wolf *et al.*, 2017](#ref_wolf17)).

* [tl.tsne](#tsne) - t-SNE ([Maaten & Hinton, 2008](#ref_maaten08); [Amir *et al.*, 2013](#ref_amir13);
  [Pedregosa *et al.*, 2011](#ref_pedregosa11)).

* [tl.spring](#spring) - Force-directed graph drawing
  ([Wikipedia;](https://en.wikipedia.org/wiki/Force-directed_graph_drawing)
  [Weinreb *et al.*, 2016](https://doi.org/10.1101/090332)).

#### Branching trajectories and pseudotime, clustering, differential expression

* [tl.dpt](#dpt) - Infer progression of cells, identify *branching*
subgroups ([Haghverdi *et al.*, 2016](#ref_haghverdi16); [Wolf *et al.*, 2017](#ref_wolf17)).

* [tl.dbscan](#dbscan) - Cluster cells into subgroups ([Ester *et al.*,
1996](#ref_ester96); [Pedregosa *et al.*, 2011](#ref_pedregosa11)).

* [tl.diffrank](#diffrank) - Rank genes according to differential
  expression ([Wolf *et al.*, 2017](#ref_wolf17)).

#### Simulation

* [tl.sim](#sim) - Simulate dynamic gene expression data ([Wittmann
*et al.*, 2009](#ref_wittmann09); [Wolf *et al.*, 2017](#ref_wolf17)).

### Basic features

The typical workflow consists of subsequent calls of data analysis tools of the form
```
sc.tl.tool(adata, **params)
```
where `adata` is an `AnnData` object and  `params` is a dictionary that stores optional parameters. Each of these calls adds annotation to an expression matrix *X*, which stores *n* *d*-dimensional gene expression measurements. By default, Scanpy tools operate *inplace* and return `None`. If you want to copy the `AnnData` object, pass the `copy` argument
```
adata_copy = sc.tl.tool(adata, copy=True, **params)
```

#### AnnData objects

Instantiate `AnnData` via
```
adata = sc.AnnData(X[, smp][, var][, add])
```
The instance `adata` stores *X* as `adata.X` and sample annotation as `adata.smp`, variable annotation as `adata.var` and additional unstructured annotation as `adata.add`. While `adata.X` is array-like and `adata.add` is a conventional dictionary, `adata.smp` and  `adata.var` are instances of a Pandas dataframe-like class, which is based on a Numpy array and, respectively.  Values can be retrieved and appended via `adata.smp['foo_key']` and `adata.var['bar_key']`. Sample and variable names can be accessed via `adata.smp_names` and `adata.var_names`, respectively. AnnData objects can be sliced like Pandas dataframes, for example, `adata = adata[:, list_of_gene_names]`. The AnnData class is similar to R's ExpressionSet ([Huber *et al.*, 2015](https://doi.org/10.1038/nmeth.3252)); the latter though is not implemented for sparse data.

#### Reading and writing data files and AnnData objects

Instead of invoking the explicit constructor, one usually calls
```
adata = sc.read(filename)
```
to initialize an AnnData object, and `sc.write(filename, adata)` to write it back to a file. Reading foresees filenames with extensions *h5*,  *xlsx*,  *mtx*,  *txt*, *csv* and others. Writing foresees writing *h5*,  *csv* and *txt*. Scanpy is smart about file storage and extensions. Instead of providing a full filename, you can provide *filekeys*. By default, Scanpy writes to `./write/filekey.h5`, an *hdf5* file, which is configurable by setting `sc.sett.writedir` and `sc.sett.file_format_data`.

#### Plotting

For each tool, there is an associated plotting function
```
sc.pl.tool(adata)
```
that retrieves and plots the elements of `adata` that were previously written by `sc.tl.tool(adata)`. To not display figures interactively but save all plots to default locations, you can set `sc.sett.savefigs = True`.
By default, figures are saved as *png* to `./figs/`. Reset `sc.sett.file_format_figs` and `sc.sett.figdir` if you want to change this. Scanpy's plotting module can be seen similar to [Seaborn](http://seaborn.pydata.org/): an extension of [matplotlib](http://matplotlib.org/) that allows visualizing certain frequent tasks with one-line commands. Detailed configuration has to be done via  matplotlib functions, which is easy as Scanpy's plotting functions usually return a `Matplotlib.Axes` object.

#### Builtin examples

Show all builtin example data using `sc.show_exdata()`
and all builtin example use cases via `sc.show_examples()`. Load annotated and preprocessed data using an *example key*, here 'paul15', via
```
adata = sc.get_example('paul15')
```
The key 'paul15' can also be used within `sc.read('paul15')` and `sc.write('paul15', adata)` to write the current state of the AnnData object to disk.

### Visualization

#### pca <a id="pca"></a>

[[source]](scanpy/tools/pca.py) Computes the PCA representation `X_pca` of data, principal components
and variance decomposition. Uses the implementation of the `scikit-learn`
package ([Pedregosa *et al.*, 2011](#ref_pedregosa11)).

#### tsne <a id="tsne"></a>

[[source]](scanpy/tools/tsne.py) Computes the tSNE representation `X_tsne` of data.

The algorithm has been introduced by [Maaten & Hinton (2008)](#ref_maaten08) and
  proposed for single-cell data by [Amir *et al.* (2013)](#ref_amir13). By
  default, Scanpy uses the implementation of the `scikit-learn` package
  ([Pedregosa *et al.*, 2011](#ref_pedregosa11)). You can achieve a huge speedup
  if you install the Multicore-TSNE package by [Ulyanov
  (2016)](https://github.com/DmitryUlyanov/Multicore-TSNE), which will be
  automatically detected by Scanpy.

#### diffmap <a id="diffmap"></a>

[[source]](scanpy/tools/diffmap.py) Computes the diffusion maps representation `X_diffmap` of data.

Diffusion maps ([Coifman *et al.*, 2005](#ref_coifman05)) has been proposed for
visualizing single-cell data by [Haghverdi *et al.*
(2015)](#ref_haghverdi15). The tool uses the adapted Gaussian kernel suggested by [Haghverdi *et
al.* (2016)](#ref_haghverdi16). The Scanpy implementation is due to [Wolf *et
al.* (2017)](#ref_wolf17).

#### spring <a id="spring"></a>

Beta version.

[[source]](scanpy/tools/spring.py) Force-directed graph drawing is a
long-established algorithm for visualizing graphs, see [Wikipedia](https://en.wikipedia.org/wiki/Force-directed_graph_drawing).
It has been suggested for visualizing single-cell data by [Weinreb *et al.* (2016)](https://doi.org/10.1101/090332).

Here, the [Fruchterman & Reingold (1991)](http://doi.org:10.1002/spe.4380211102)
algorithm is used. The implementation uses elements of the NetworkX
implementation [(Hagberg
*et al.*, 2008)](http://conference.scipy.org/proceedings/SciPy2008/paper_2/).

### Discrete clustering of subgroups and continuous progression through subgroups

#### dpt <a id="dpt"></a>

[[source]](scanpy/tools/dpt.py) Reconstruct the progression of a biological process
from snapshot data and detect branching subgroups. Diffusion Pseudotime analysis
has been introduced by [Haghverdi *et al.* (2016)](#ref_haghverdi16) and
implemented for Scanpy by [Wolf *et al.* (2017)](#ref_wolf17).

The functionality of diffmap and dpt compare to the R package
[destiny](http://bioconductor.org/packages/release/bioc/html/destiny.html) of
[Angerer *et al.* (2015)](#ref_angerer16).

*Examples:* See one of the early examples [[notebook](https://github.com/theislab/scanpy_usage/tree/master/170503_moignard15.ipynb)/[command line](https://github.com/theislab/scanpy_usage/tree/master/EXAMPLES.md#moignard15)] dealing with data of [Moignard *et al.*, Nat. Biotechn. (2015)](http://doi.org/10.1038/nbt.3154).

#### dbscan <a id="dbscan"></a>

[[source]](scanpy/tools/dbscan.py) Cluster cells using [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) ([Ester *et al.*, 1996](#ref_ester96)), in the implementation of
`scikit-learn` ([Pedregosa *et al.*, 2011](#ref_pedregosa11)).

This is a very simple clustering method. A better one - in the same framework as DPT and Diffusion Maps - will come soon.

### Differential expression

#### diffrank <a id="diffrank"></a>

[[source]](scanpy/tools/diffrank.py) Rank genes by differential expression.

### Simulation

#### sim <a id="sim"></a>

[[source]](scanpy/tools/sim.py) Sample from a stochastic differential equation
model built from literature-curated boolean gene regulatory networks, as
suggested by [Wittmann *et al.* (2009)](#ref_wittmann09). The Scanpy implementation is
due to [Wolf *et al.* (2017)](#ref_wolf17).

The tool compares to the Matlab tool *Odefy* of [Krumsiek *et al.*
(2010)](#ref_krumsiek10).

## Installation <a id="install"></a>

If you do not have a current Python distribution (Python 3.5 or 3.6), download
and install [Miniconda](http://conda.pydata.org/miniconda.html) (see below).

Then, download or clone the repository - green button on top of the page - and `cd`
into its root directory. To install with symbolic links (stay up to date with
your cloned version after you update with `git pull`) call
```
pip install -e .
```
and work with the top-level command `scanpy` or
```python
import scanpy as sc
```
in any directory.

#### Installing Miniconda

After downloading [Miniconda](http://conda.pydata.org/miniconda.html), in a unix shell (Linux, Mac), run
```shell
cd DOWNLOAD_DIR
chmod +x Miniconda3-latest-VERSION.sh
./Miniconda3-latest-VERSION.sh
```
and accept all suggestions. Either reopen a new terminal or `source
~/.bashrc` on Linux/ `source ~/.bash_profile` on Mac. The whole process takes just a couple of minutes.

#### PyPi

The package is [registered](https://pypi.python.org/pypi/scanpy/0.1) in the
[Python Packaging Index](https://pypi.python.org/pypi), but versioning has not
started yet. In the future, installation will also be possible without reference
to GitHub via `pip install scanpy`.

## References <a id="references"></a>

<a id="ref_amir13"></a>
Amir *et al.* (2013),
*viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia*
[Nature Biotechnology 31, 545](
http://dx.doi.org/10.1038/nbt.2594).

<a id="ref_angerer15"></a>
Angerer *et al.* (2015),
*destiny - diffusion maps for large-scale single-cell data in R*,
[Bioinformatics 32, 1241](
 https://doi.org/10.1093/bioinformatics/btv715).

<a id="ref_coifman05"></a>
Coifman *et al.* (2005),
*Geometric diffusions as a tool for harmonic analysis and structure definition of data: Diffusion maps*,
[PNAS 102, 7426](
http://www.pnas.org/content/102/21/7432.full).

<a id="ref_ester96"></a>
Ester *et al.* (1996),
*A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise*
[Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining, Portland, OR, pp. 226-231](
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.121.9220).

<a id="ref_haghverdi15"></a>
Haghverdi *et al.* (2015),
*Diffusion maps for high-dimensional single-cell analysis of differentiation data*,
[Bioinformatics 31, 2989](
http://dx.doi.org/10.1093/bioinformatics/btv325).

<a id="ref_haghverdi16"></a>
Haghverdi *et al.* (2016),
*Diffusion pseudotime robustly reconstructs branching cellular lineages*,
[Nature Methods 13, 845](
http://dx.doi.org/10.1038/nmeth.3971).

<a id="ref_krumsiek10"></a>
Krumsiek *et al.* (2010),
*Odefy - From discrete to continuous models*,
[BMC Bioinformatics 11, 233](http://dx.doi.org/10.1186/1471-2105-11-233).

<a id="ref_krumsiek11"></a>
Krumsiek *et al.* (2011),
*Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*,
[PLoS ONE 6, e22649](http://dx.doi.org/10.1371/journal.pone.0022649).

<a id="ref_maaten08"></a>
Maaten & Hinton (2008),
*Visualizing data using t-SNE*,
[JMLR 9, 2579](
http://www.jmlr.org/papers/v9/vandermaaten08a.html).

<a id="ref_moignard15"></a>
Moignard *et al.* (2015),
*Decoding the regulatory network of early blood development from single-cell gene expression measurements*,
[Nature Biotechnology 33, 269](
http://dx.doi.org/10.1038/nbt.3154).

<a id="ref_pedregosa11"></a>
Pedregosa *et al.* (2011),
*Scikit-learn: Machine Learning in Python*,
[JMLR 12, 2825](
http://www.jmlr.org/papers/v12/pedregosa11a.html).

<a id="ref_paul15"></a>
Paul *et al.* (2015),
*Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*,
[Cell 163, 1663](
http://dx.doi.org/10.1016/j.cell.2015.11.013).

<a id="ref_wittmann09"></a>
Wittmann *et al.* (2009),
*Transforming Boolean models to continuous models: methodology and application to T-cell receptor signaling*,
[BMC Systems Biology 3, 98](
http://dx.doi.org/10.1186/1752-0509-3-98).
