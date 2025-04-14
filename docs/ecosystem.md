# Ecosystem

```{warning}
We are no longer accepting new tools on this page.
Instead, please submit your tool to the [scverse ecosystem package listing](https://scverse.org/packages/#ecosystem).
```

## Viewers

Interactive manifold viewers.

- [cellxgene](https://github.com/chanzuckerberg/cellxgene) via direct reading of `.h5ad` {small}`CZI`
- [cirrocumulus](https://cirrocumulus.readthedocs.io/) via direct reading of `.h5ad` {small}`Broad Inst.`
- [cell browser](https://cells.ucsc.edu/) via exporing through {func}`~scanpy.external.exporting.cellbrowser` {small}`UCSC`
- [SPRING](https://github.com/AllonKleinLab/SPRING) via exporting through {func}`~scanpy.external.exporting.spring_project` {small}`Harvard Med`
- [vitessce](https://github.com/vitessce/vitessce#readme) for purely browser based viewing of zarr formatted AnnData files {smaller}`Harvard Med`

## Portals

- the [Gene Expression Analysis Resource](https://umgear.org/) {small}`U Maryland`
- the [Galaxy Project](https://humancellatlas.usegalaxy.eu) for the Human Cell Atlas [\[tweet\]](https://twitter.com/ExpressionAtlas/status/1151797848469626881) {small}`U Freiburg`
- the [Expression Atlas](https://www.ebi.ac.uk/gxa/sc/help.html) {small}`EMBL-EBI`

## Modalities

### RNA velocity

- [scVelo](https://scvelo.org) {small}`Helmholtz Munich`

### Spatial Transcriptomics Tools

- [squidpy](https://squidpy.readthedocs.io/en/stable/) {small}`Helmholtz Munich`

  > Squidpy is a comprehensive toolkit for working with spatial single cell omics data.

- [PASTE](https://github.com/raphael-group/paste) {small}`Princeton`

  > PASTE is a computational method to align and integrate spatial transcriptomics data across adjacent tissue slices by leveraging both gene expression similarity and spatial distances between spots.

- [bento](https://bento-tools.readthedocs.io/en/latest/) 🍱 {small}`UC San Diego`

  > Bento is an accessible Python toolkit for performing subcellular analysis of spatial transcriptomics data.

### Multimodal integration

- [MUON](https://muon.readthedocs.io/en/latest/) and [MuData](https://mudata.readthedocs.io/en/latest/) {small}`EMBL/ DKFZ`

  > MUON, and it's associated data structure MuData are designed to organise, analyse, visualise, and exchange multimodal data.
  > MUON enables a range of analyses for ATAC and CITE-seq, from data preprocessing to flexible multi-omics alignment.

### Adaptive immune receptor repertoire (AIRR)

- [scirpy](https://github.com/icbi-lab/scirpy) {small}`Medical University of Innsbruck`

  > scirpy is a scanpy extension to expore single-cell T-cell receptor (TCR) and B-cell receptor (BCR) repertoires.

- [dandelion](https://github.com/zktuong/dandelion) {small}`University of Cambridge`

  > dandelion is a single-cell BCR-seq network analysis package that integrates with transcriptomic data analyzed via scanpy.

### Long reads

- [Swan](https://freese.gitbook.io/swan/tutorials/data_processing) {small}`UC Irvine`

  > Swan is a Python library designed for the analysis and visualization of transcriptomes, especially with long-read transcriptomes in mind.
  > Users can add transcriptomes from different datasets and explore distinct splicing and expression patterns across datasets.

## Analysis methods

### scvi-tools

- [scvi-tools](https://github.com/YosefLab/scvi-tools) {small}`Berkeley`

  > scvi-tools hosts deep generative models (DGM) for end-to-end analysis of single-cell
  > omics data (e.g., scVI, scANVI, totalVI). It also contains several primitives to build novel DGMs.

### Fate mapping

- [CellRank](https://cellrank.org) {small}`Helmholtz Munich`

  > CellRank is a framework to uncover cellular dynamics based on single-cell data.
  > It incorporates modalities such as RNA velocity, pseudotime, developmental potential, real-time information, etc.

### Differential expression

- [diffxpy](https://github.com/theislab/diffxpy) {small}`Helmholtz Munich`

(eco-data-integration)=

### Data integration

- [scanaroma](https://github.com/brianhie/scanorama) {small}`MIT`

### Modeling perturbations

- [scGen](https://github.com/theislab/scgen) / [trVAE](https://github.com/theislab/trvae) {small}`Helmholtz Munich`

### Feature selection

- [triku 🦔](https://gitlab.com/alexmascension/triku) {small}`Biodonostia Health Research Institute`
- [CIARA](https://github.com/ScialdoneLab/CIARA_python) {small}`Helmholtz Munich`

  > CIARA is an algorithm for feature selection, that aims for the identification of rare cell types via scRNA-Seq data in scanpy.

### Annotation/ Enrichment Analysis

Analyses using curated prior knowledge

- [decoupler](https://github.com/saezlab/decoupler-py) is a collection of footprint enrichment methods that allows to infer transcription factor or pathway activities. {small}`Institute for Computational Biomedicine, Heidelberg University`
- [Cubé](https://github.com/connerlambden/Cube) {small}`Harvard University`

  > Intuitive Nonparametric Gene Network Search Algorithm that learns from existing biological pathways & multiplicative gene interference patterns.
