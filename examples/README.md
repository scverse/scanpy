# Examples

<!--- ----------------------------------------------------------------------- -->
<!--- call "./scanpy.py exdata markup" to get the list below          -->

#### Example Data

* [moignard15](#moignard15) - [Moignard *et al.*, Nature Biotechnology 33, 269 (2015)](http://dx.doi.org/10.1038/nbt.3154)   
*Decoding the regulatory network of early blood development from single-cell gene expression measurements*
* [paul15](#paul15) - [Paul *et al.*, Cell 163, 1663 (2015)](http://dx.doi.org/10.1016/j.cell.2015.11.013)   
*Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors*

Examples using simulated data.
* [krumsiek11](#krumsiek11) - [Krumsiek *et al.*, PLoS ONE 6, e22649 (2011)](http://dx.doi.org/10.1371/journal.pone.0022649)   
*Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network*
* [toggleswitch](#toggleswitch)   
*Simple toggle switch model.*   


<!--- ----------------------------------------------------------------------- -->
<!--- start the technical topics list here                                    -->

#### Technical Topics

* preprocessing with PCA - [paul15pca](#paul15pca)



<!--- ----------------------------------------------------------------------- -->
<!--- list the actual description of examples here                            -->

#### Data of [Moignard *et al.* (2015)](#ref_moignard15) <a id="moignard15"></a>

See Scanpy [README.md](https://github.com/theislab).

#### Data of [Paul *et al.* (2015)](#ref_paul15)

Load the data in `user_exs.py` as follows.
```python
def paul15():
    adata = paul15_raw()
    adata.X = sc.pp.log(adata.X)
    adata['xroot'] = adata.X[adata['iroot']]
    return adata

def paul15_raw():
    filename = 'data/paul15/paul15.h5'
    ddata = sc.read(filename, 'data.debatched')
    adata = sc.AnnData(ddata)
    # the data has to be transposed (in the hdf5 and R files, each row
    # corresponds to one gene, we use the opposite convention)
    adata = adata.transpose()
    # cluster assocations identified by Paul et al.
    # groups = sc.read(filename,'cluster.id')['X']
    infogenes_names = sc.read(filename, 'info.genes_strings')['X']
    # just keep the first of the two equivalent names per gene
    adata.var_names = np.array([gn.split(';')[0] for gn in adata.var_names])
    # index array for the informative genes
    infogenes_idcs = np.array([gn in infogenes_names for gn in adata.var_names])
    # restrict data array to the 3451 informative genes
    adata = adata[:, infogenes_idcs]
    # set root cell as in Haghverdi et al. (2016)
    adata['iroot'] = iroot = 840 # note that in Matlab/R, counting starts at 1
    adata['xroot'] = adata.X[iroot]
    return adata
```

##### DPT analysis

Diffusion Pseudotime (DPT) analysis detects the branch of granulocyte/macrophage
progenitors (GMP), and the branch of megakaryocyte/erythrocyte progenitors
(MEP). There are two small further subgroups (*segments* 0 and 2).
```
scanpy paul15 dpt
```
<img src="http://falexwolf.de/scanpy/figs0/paul15_dpt_diffmap_pseudotimes_segments.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/paul15_dpt_segpt.png" height="175">

We can now test for differential gene expression.
```
scanpy paul15 difftest --prev dpt
```
<img src="http://falexwolf.de/scanpy/figs/paul15_difftest.png" height="175">

See the [notebook](examples/paul15.ipynb) for more information.

##### Preprocessing with PCA <a id="paul15pca"></a>

As for most methods that build a data graph, i.e., compute a distance matrix, preprocessing by PCA speeds up
computations tremendously. Distance matrix computations in high dimensions are a very demanding problems.
```python
def paul15pca():
    adata = paul15_raw()
    adata.X = sc.pp.log(adata.X)
    # reduce to 50 components
    adata['Xpca'] = sc.pca(adata.X, nr_comps=50)
    # adjust expression vector of root cell
    adata['xroot'] = adata['Xpca'][adata['iroot']]
    return adata    
```

Does this preprocessing change the biology?
```
scanpy paul15pca dpt
scanpy paul15pca dpt difftest --prev dpt
```
<img src="http://falexwolf.de/scanpy/figs0/paul15pca_dpt_diffmap_pseudotimes_segments.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/paul15pca_difftest.png" height="175">

Even though the diffmap representation has changed a lot, the differential genes are still the same.
Let us check whether the subgroups look much different in tSNE.
```
scanpy paul15 dpt tsne
scanpy paul15pca dpt tsne
```
<img src="http://falexwolf.de/scanpy/figs0/paul15_dpt_tsne_pseudotimes_segments.png" height="175">
<img src="http://falexwolf.de/scanpy/figs0/paul15pca_dpt_tsne_pseudotimes_segments.png" height="175">

Even though diffusion components 1 and 2 were considerably affected by preprocessing
with PCA, subgroup identification has not!

Note that, as tSNE yields non-deterministic output, we used the file
`write/paul15_tsne.h5` also for plotting `paul15pca` (`cp write/paul15_tsne.h5
write/paul15pca_tsne.h5`).

#### Simulated myeloid progenitor data ([Krumsiek *et al.*, 2011](#ref_krumsiek11)) <a id="krumsiek11"></a>

Here, we are going to simulate some data using a literature-curated boolean gene
regulatory network, which is believed to describe myeloid differentiation
([Krumsiek *et al.*, 2011](#ref_krumsiek11)). Using [sim.py](scanpy/sim.py), the
[boolean model](models/krumsiek11.txt)
```
Gata2 = Gata2 and not (Gata1 and Fog1) and not Pu.1
Gata1 = (Gata1 or Gata2 or Fli1) and not Pu.1
Fog1 = Gata1
EKLF = Gata1 and not Fli1
Fli1 = Gata1 and not EKLF
SCL = Gata1 and not Pu.1
Cebpa = Cebpa and not (Gata1 and Fog1 and SCL)
Pu.1 = (Cebpa or Pu.1) and not (Gata1 or Gata2)
cJun = Pu.1 and not Gfi1
EgrNab = (Pu.1 and cJun) and not Gfi1
Gfi1 = Cebpa and not EgrNab
```
is translated into a stochastic differential equation ([Wittmann *et al.*, 2009](#ref_wittmann09)). Simulations result
in branching time series of gene expression, where each branch corresponds to a
certain cell fate of common myeloid progenitors (megakaryocytes, erythrocytes,
granulocytes and monocytes).
```
scanpy krumsiek11 sim
```
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_sim.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_sim_shuffled.png" height="175">

If the order is shuffled, as in a snapshot, the same data looks as on the right.
Let us reconstruct the process using DPT and obtain the branching lineages
```
scanpy krumsiek11 dpt -p layout 3d
```
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_vsorder.png" height="175">
<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_segpt.png" height="175">

The left panel illustrates how the data is organized according to a *pseudotime*
and different *segments*. Pseudotime 'estimates geodesic distance on the
manifold' from a root cell. Segments are discrete partitions of the data. 

<img src="http://falexwolf.de/scanpy/figs/krumsiek11_dpt_diffmap_new.png" height="175">

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

