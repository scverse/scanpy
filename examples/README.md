# Examples

<!--- call "./scripts/scanpy.py exdata markup" to get this output -->
Examples using experimental data.
* [moignard15](#moignard15) - [Moignard *et al.*, Nature Biotechnology 33, 269 (2015)](http://dx.doi.org/10.1038/nbt.3154)   
*Decoding the regulatory network of early blood development from single-cell gene expression measurements*
* [paul15](#paul15) - [Paul *et al.*, Cell 163, 1663 (2015)](http://dx.doi.org/10.1016/j.cell.2015.11.013)   
*Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*

Examples using simulated data.
* [krumsiek11](#krumsiek11) - [Krumsiek *et al.*, PLoS ONE 6, e22649 (2011)](http://dx.doi.org/10.1371/journal.pone.0022649)   
*Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*
* [toggleswitch](#toggleswitch)   
*Simple toggle switch model.*   


#### Data of [Moignard *et al.* (2015)](#ref_moignard15) <a id="moignard15"></a>

Early mesoderm cells in mouse differentiate through three subsequent stages (PS,
NP, HF) and then branch into erythorytes (4SG) and endothelial cells (4SFG).
```shell
python scripts/scanpy.py moignard15 pca
python scripts/scanpy.py moignard15 tsne
python scripts/scanpy.py moignard15 diffmap
```
<img src="http://falexwolf.de/scanpy/figs/moignard15_pca.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/moignard15_tsne.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/moignard15_diffmap.png" height="175">

Diffusion Pseudotime (DPT) analysis reveals differentation and branching. It
detects the *trunk* of progenitor cells (segment 0) and the *branches* of endothelial
cells (segment 1/2) and erythrocytes (segment 3). The inferred *pseudotime*
traces the degree of cells' progression in the differentiation process.
```shell
python scripts/scanpy.py moignard15 dpt
```
<img src="http://falexwolf.de/scanpy/figs/moignard15_dpt_diffmap.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/moignard15_dpt_segpt.png" height="175">

This orders cells by segment, and within each segment, by pseudotime and outputs

<img src="http://falexwolf.de/scanpy/figs/moignard15_dpt_heatmap.png" height="250">

With this, we reproduced Fig. 1 from [Haghverdi *et al.*
(2016)](#ref_haghverdi16). See this [notebook](examples/moignard15.ipynb) for
more information.

#### Data of [Paul *et al.* (2015)](#ref_paul15)

Diffusion Pseudotime (DPT) analysis detects the branch of granulocyte/macrophage
progenitors (GMP), and the branch of megakaryocyte/erythrocyte progenitors
(MEP). There are two small further subgroups (*segments* 0 and 2).
```shell
python scripts/scanpy.py paul15 dpt
```
<img src="http://falexwolf.de/scanpy/figs/paul15_dpt_diffmap.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/paul15_dpt_segpt.png" height="175">

We can now test for differential gene expression.
```shell
python scripts/scanpy.py paul15 difftest
```
<img src="http://falexwolf.de/scanpy/figs/paul15_difftest.png" height="175">

See the [notebook](examples/paul15.ipynb) for more information.

#### Simulated myeloid progenitor data ([Krumsiek *et al.*, 2011](#ref_krumsiek11)) <a id="krumsiek11"></a>

Here, we are going to simulate some data using a literature-curated boolean gene
regulatory network, which is believed to describe myeloid differentiation
([Krumsiek *et al.*, 2011](#ref_krumsiek11)). Using [sim.py](scanpy/sim.py), the
[boolean model](models/krumsiek11.txt) is translated into a stochastic 
differential equation ([Wittmann *et al.*, 2009](#ref_wittmann09)). Simulations result
in branching time series of gene expression, where each branch corresponds to a
certain cell fate of common myeloid progenitors (megakaryocytes, erythrocytes,
granulocytes and monocytes).
```shell
python scripts/scanpy.py krumsiek11 sim
```
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_sim.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_sim_shuffled.png" height="175">

If the order is shuffled, as in a snapshot, the same data looks as on the right.
Let us reconstruct the process using DPT and obtain the branching lineages
```shell
python scripts/scanpy.py krumsiek11 dpt --plotparams layout 3d
```
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_vsorder.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_segpt.png" height="175">

The left panel illustrates how the data is organized according to a *pseudotime*
and different *segments*. Pseudotime 'estimates geodesic distance on the
manifold' from a root cell. Segments are discrete partitions of the data. 

<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_diffmap.pngRELOAD" height="175">

See the [notebook](examples/krumsiek11.ipynb) for more.


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

