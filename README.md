[Quick Start](#quick_start) | 
[Example Library](examples) | 
[Tools](#tools) | 
[Installation](#install) | 
[References](#references)

# Scanpy - Single-Cell Analysis in Python

Tools for analyzing and simulating single-cell data that aim at an understanding
of dynamic biological processes from snapshots of transcriptome or proteome.

* [pca](scanpy/tools/pca.py) - Visualize data using PCA.

* [diffmap](#diffmap) - Visualize data using Diffusion Maps
([Coifman *et al.*, 2005](#ref_coifman05); [Haghverdi *et al.*,
2015](#ref_haghverdi15)).

* [tsne](scanpy/tools/tsne.py) - Visualize data using t-SNE ([van
  der Maaten & Hinton, 2008](#ref_vandermaaten08); [Amir *et al.*, 2013](#ref_amir13)).

* [dpt](#dpt) - Infer progression of cells, identify *branching*
subgroups ([Haghverdi *et al.*, 2016](#ref_haghverdi16)).

*  [difftest](scanpy/tools/difftest.py) - Test for differential gene
  expression.

* [sim](#sim) - Simulate dynamic gene expression data ([Wittmann
*et al.*, 2009](#ref_wittmann09)).

Please, cite the original references. The draft [Wolf & Theis
(2017)](http://falexwolf.de/docs/scanpy.pdf) explains conceptual ideas and usage
as a library. Potential coauthors who would like to work on software and
manuscript are welcome! Contributors, who just want to add their own example or
tool are welcome, too! Any comments are appreciated!

## Quick Start <a id="quick_start"></a>

Download or clone the repository - green button on top of the page - and `cd`
into its root directory. Type `pip install -e .` and you can immediately work
with the top-level command `scanpy` in any directory (more info [here](#install)).

Packages you might need (all default in
[Anaconda](https://www.continuum.io/downloads)) can be easily installed using
[Miniconda](http://conda.pydata.org/miniconda.html). Then run `conda install
scipy matplotlib h5py scikit-learn pandas xlrd`. Scanpy is written in Python 3
and compatible with Python 2.

#### Data of [Moignard *et al.* (2015)](#ref_moignard15) <a id="moignard15"></a>

Early mesoderm cells in mouse differentiate through three subsequent stages (PS,
NP, HF) and then branch into erythorytes (4SG) and endothelial cells (4SFG).
```
scanpy moignard15 pca
scanpy moignard15 tsne
scanpy moignard15 diffmap
```
<img src="http://falexwolf.de/scanpy/figs0/moignard15_pca_groups.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_tsne_groups.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_diffmap_groups.png" height="175">

Coloring samples/cells by gene expression works analogously,
```
scanpy moignard15 pca -p smp HbbbH1
scanpy moignard15 tsne -p smp HbbbH1
scanpy moignard15 diffmap -p smp HbbbH1
```
<img src="http://falexwolf.de/scanpy/figs0/moignard15_pca_HbbbH1.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_tsne_HbbbH1.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_diffmap_HbbbH1.png" height="175">

Diffusion Pseudotime (DPT) analysis reveals differentation and branching. It
detects the *trunk* of progenitor cells (segment 0) and the *branches* of endothelial
cells (segment 1/2) and erythrocytes (segment 3). The inferred *pseudotime*
traces the degree of cells' progression in the differentiation process. By default,
this is plotted using Diffusion Maps, but you might just as well plot the subgroups using another
plotting tool.
```
scanpy moignard15 dpt
scanpy moignard15 dpt tsne
```
<img src="http://falexwolf.de/scanpy/figs0/moignard15_dpt_diffmap_pseudotimes_segments_groups.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_dpt_segpt.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/moignard15_dpt_tsne_pseudotimes_segments_groups.png" height="175">

This orders cells by segment, and within each segment, by pseudotime.

<img src="http://falexwolf.de/scanpy/figs0/moignard15_dpt_heatmap.png" height="250">

With this, we reproduced most of Fig. 1 from [Haghverdi *et al.*
(2016)](#ref_haghverdi16). See this [notebook](examples/moignard15.ipynb) for
more information.

If you want to use the results externally, simply read the corresponding hdf5
file, which you can inspect using `h5ls write/moignard15_dpt.h5`. If you prefer
reading and writing csv/txt files, which is much slower, however, use the
option `--fileformat csv` (or `txt`).

#### More examples and help

For more examples, read [this](examples), or display them on the command line
(example data and example use cases, respectively).
```shell
scanpy exdata
scanpy examples
```

Get general and tool-specific help, respectively.
```shell
scanpy --help
scanpy dpt --help
```

#### Work on your own examples <a id="add_example"></a>

To work on your own example, make a copy and edit the following
[notebook](examples/myexample_template.ipynb). If you want to call user examples
from the command-line, create a file `user_exs.py` in your current working
directory, e.g., by downloading and renaming
[user_exs_template.py](scanpy/exs/user_exs_template.py) and changing the function
`myexample()` to your needs. Consider using copy and paste from
[scanpy/exs/builtin.py](scanpy/exs/builtin.py). Call your example using `scanpy
myexample pca`. For the previous example (`moignard15`) you would define the
following
```python
def moignard15():
    filename = 'data/moignard15/nbt.3154-S3.xlsx'
    ddata = sc.read(filename, sheet='dCt_values.txt')
    adata = sc.AnnData(ddata)
    # filter out genes: the 4th column (Eif2b1), the 31nd (Mrpl19), the 36th
    # (Polr2a) and the 45th (last,UBC), as done by Haghverdi et al. (2016)
    genes = np.array([g not in [4, 31, 36, 45] for g in range(adata.X.shape[1])])
    adata = adata[:, genes] # filter adata
    # choose root cell as in Haghverdi et al. (2016)
    adata['iroot'] = iroot = 532 # note that in Matlab/R, counting starts at 1
    adata['xroot'] = adata.X[iroot]
    # annotate with Moignard et al. (2015) experimental cell groups
    groups_names = ['HF', 'NP', 'PS', '4SG', '4SFG']
    # annotate each sample/cell
    adata.smp['groups'] = [
        next(gname for gname in groups_names if sname.startswith(gname))
        for sname in adata.smp_names]
    # fix the order and colors of names in "groups"
    adata['groups_names'] = groups_names
    adata['groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    return adata
```

Also, it'd be awesome if you provide you upload your example to the
[examples](examples) library: copy your example from `user_exs.py` to
[scanpy/exs/builtin.py](scanpy/exs/builtin.py) with a link to the public data
and make a pull request. If you have questions or prefer sending your script by
email, contact Alex.

## Tools <a id="tools"></a>

### [diffmap](scanpy/tools/diffmap.py) <a id="diffmap"></a>

This implements *diffusion maps* ([Coifman
*et al.*, 2005](#ref_coifman05)), which has been proposed for visualizing
single-cell data by [Haghverdi *et al.* (2015)](#ref_haghverdi15). Also, it uses
the kernel suggested by [Haghverdi *et al.* (2016)](#ref_haghverdi16).

### [dpt](scanpy/tools/dpt.py) <a id="dpt"></a>

Reconstruct progression in a biological process from snapshot data and detect
branching subgroups. Diffusion Pseudotime analysis has been introduced by
[Haghverdi *et al.* (2016)](#ref_haghverdi16).

The functionality of diffmap and dpt compare to the R package
[destiny](http://bioconductor.org/packages/release/bioc/html/destiny.html) of
[Angerer *et al.* (2015)](#ref_angerer16).

### [sim](scanpy/tools/sim.py) <a id="sim"></a>

Sample from a stochastic differential equation model built from
literature-curated boolean gene regulatory networks, as suggested by [Wittmann
*et al.* (2009)](#ref_wittmann09).

The tool compares to the Matlab tool *Odefy* of [Krumsiek *et al.*
(2010)](#ref_krumsiek10).

## Installation <a id="install"></a>

Download or clone the repository - green button on top of the page - and `cd`
into its root directory. To install with symlinks (stay up to date with your
cloned version after you update with `git pull`) with ```pip install -e .``` and
work with the top-level command `scanpy` in any directory.

Or, call the wrapper `python scanpy.py` from within the root of the repository,
which works **without** [installation](#install).

Packages you might need (all default in
[Anaconda](https://www.continuum.io/downloads)) can be easily installed via 
[Miniconda](http://conda.pydata.org/miniconda.html). Run `conda install scipy
matplotlib h5py pandas xlrd scikit-learn`. Scanpy is written in Python 3 and compatible with
Python 2.

#### Package Managment via Miniconda

After downloading [Miniconda](http://conda.pydata.org/miniconda.html), in a unix shell (Linux, Mac), run
```shell
cd DOWNLOAD_DIR
chmod +x Miniconda3-latest-VERSION.sh
./Miniconda3-latest-VERSION.sh
```
and accept all suggestions. Either reopen a new terminal or `source
~/.bashrc` on Linux/ `source ~/.bash_profile` on Mac. Then run
`conda install scipy matplotlib h5py pandas xlrd scikit-learn`. The whole process takes about 5 min.

<!--- You can install Scanpy, hence `import scanpy` from anywhere on your system, via
```shell
pip install .
```
Then, `scanpy` becomes a top-level command available from anywhere in the
system, call `scanpy --help`. The package is
[registered](https://pypi.python.org/pypi/scanpy/0.1) in the [Python Packaging
Index](https://pypi.python.org/pypi), but versioning has not started yet. In the
future, installation will be possible without reference to GitHub via 
`pip install scanpy`.
-->

## References <a id="references"></a>

<a id="ref_amir13"></a>
Amir *et al.* (2013), 
*viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia*
[Nature Biotechnology 31, 545](http://dx.doi.org/10.1038/nbt.2594).

<a id="ref_angerer15"></a>
Angerer *et al.* (2015), *destiny - diffusion maps for large-scale single-cell
data in R*, [Bioinformatics 32, 1241](http://dx.doi.org/10.1038/nmeth.3971).

<a id="ref_coifman05"></a>
Coifman *et al.* (2005), *Geometric diffusions as a tool for harmonic analysis and
structure definition of data: Diffusion maps*, [PNAS 102, 7426](
http://dx.doi.org/10.1038/nmeth.3971).

<a id="ref_haghverdi15"></a>
Haghverdi *et al.* (2015), *Diffusion maps for high-dimensional single-cell
analysis of differentiation data*, [Bioinformatics 31, 2989](
http://dx.doi.org/10.1093/bioinformatics/btv325).

<a id="ref_haghverdi16"></a>
Haghverdi *et al.* (2016), *Diffusion pseudotime robustly
reconstructs branching cellular lineages*, [Nature Methods 13, 845](
http://dx.doi.org/10.1038/nmeth.3971).

<a id="ref_krumsiek10"></a>
Krumsiek *et al.* (2010), 
*Odefy - From discrete to continuous models*, 
[BMC Bioinformatics 11, 233](http://dx.doi.org/10.1186/1471-2105-11-233).

<a id="ref_krumsiek11"></a>
Krumsiek *et al.* (2011), 
*Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*, 
[PLoS ONE 6, e22649](http://dx.doi.org/10.1371/journal.pone.0022649).

<a id="ref_moignard15"></a>
Moignard *et al.* (2015), *Decoding the regulatory network of early blood
development from single-cell gene expression measurements*, [Nature Biotechnology 33,
269](
http://dx.doi.org/10.1038/nbt.3154).

<a id="ref_paul15"></a>
Paul *et al.* (2015), *Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*, [Cell 163,
1663](
http://dx.doi.org/10.1016/j.cell.2015.11.013).

<a id="ref_vandermaaten08"></a>
van der Maaten & Hinton (2008), *Visualizing data using t-SNE*, [JMLR 9, 2579](
http://www.jmlr.org/papers/v9/vandermaaten08a.html).

<a id="ref_wittmann09"></a>
Wittmann *et al.* (2009), *Transforming Boolean models to
continuous models: methodology and application to T-cell receptor signaling*,
[BMC Systems Biology 3, 98](
http://dx.doi.org/10.1186/1752-0509-3-98).

