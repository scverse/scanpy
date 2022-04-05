# Tutorials

## Clustering

For getting started, we recommend Scanpy’s reimplementation {tutorial}`pbmc3k`
of Seurat’s [^cite_satija15] clustering tutorial for 3k PBMCs from 10x Genomics,
containing preprocessing, clustering and the identification of cell types via
known marker genes.

```{image} _static/img/tutorials/170505_seurat/filter_genes_dispersion.png
:width: 100px
```

```{image} _static/img/tutorials/170505_seurat/louvain.png
:width: 100px
```

```{image} _static/img/tutorials/170505_seurat/NKG7.png
:width: 100px
```

```{image} _static/img/tutorials/170505_seurat/violin.png
:width: 100px
```

```{image} _static/img/tutorials/170505_seurat/cell_types.png
:width: 200px
```

## Visualization

This tutorial shows how to visually explore genes using scanpy. {tutorial}`plotting/core`

```{image} _static/img/stacked_violin_dotplot_matrixplot.png
:width: 550px
```

## Trajectory inference

Get started with the following example for hematopoiesis for data of [^cite_paul15]: {tutorial}`paga-paul15`

```{image} _static/img/tutorials/paga_paul15.png
:width: 450px
```

More examples for trajectory inference on complex datasets can be found in the
[PAGA](https://github.com/theislab/paga) repository [^cite_wolf19], for instance, multi-resolution analyses of whole
animals, such as for [planaria] for data of [^cite_plass18].

```{image} _static/img/tutorials/paga_planaria.png
:width: 350px
```

As a reference for simple pseudotime analyses, we provide the diffusion pseudotime (DPT) analyses of [^cite_haghverdi16]
for two hematopoiesis datasets: [DPT example 1] [^cite_paul15] and [DPT example 2] [^cite_moignard15].

## Integrating datasets

Map labels and embeddings of reference data to new data: {tutorial}`integrating-data-using-ingest`

```{image} https://scanpy-tutorials.readthedocs.io/en/latest/_images/integrating-data-using-ingest_21_0.png
:width: 350px
```

## Spatial data

- Basic analysis of spatial data: {tutorial}`spatial/basic-analysis`
- Integrating spatial data with scRNA-seq using scanorama: {tutorial}`spatial/integration-scanorama`

```{image} _static/img/spatial-basic-analysis.png
:width: 250px
```

## Further Tutorials

(conversion-to-r)=

### Conversion: AnnData, SingleCellExperiment, and Seurat objects

```{image} https://github.com/scverse/scanpy-in-R/raw/master/logo.png
:align: right
:width: 200px
```

- See [Seurat to AnnData] for a tutorial on `anndata2ri`.
- See the [Scanpy in R] guide for a tutorial on interacting with Scanpy from R.

### Regressing out cell cycle

See the [cell cycle] notebook.

```{image} _static/img/tutorials/170522_visualizing_one_million_cells/tsne_1.3M.png
:align: right
:width: 120px
```

### Normalization with Pearson Residuals

Normalization of scRNA-seq data with Pearson Residuals, from [^cite_lause21]: {tutorial}`tutorial_pearson_residuals`

### Scaling Computations

- Visualize and cluster [1.3M neurons] from 10x Genomics.

### Simulations

Simulating single cells using literature-curated gene regulatory networks [^cite_wittmann09].

```{image} _static/img/tutorials/170430_krumsiek11/timeseries.png
:align: right
:width: 200px
```

- Notebook for [myeloid differentiation]
- Notebook for simple [toggleswitch]

### Images

See pseudotime-time inference on deep-learning based features for [cell cycle reconstruction] from image data [^cite_eulenberg17].

% User Examples
% ~~~~~~~~~~~~~
%
% January 12, 2018: `Exploring the mouse cell atlas`_ by `David P. Cook`_.
% Data by `Tabula Muris Consortium`_.
%
% .. _Exploring the mouse cell atlas: https://github.com/dpcook/fun_analysis/blob/master/tabula_muris/mouse_atlas_scanpy.ipynb
% .. _David P. Cook: https://twitter.com/DavidPCook
% .. _Tabula Muris Consortium: https://www.biorxiv.org/content/early/2017/12/20/237446

[1.3m neurons]: https://github.com/scverse/scanpy_usage/tree/master/170522_visualizing_one_million_cells
[cell cycle]: https://nbviewer.jupyter.org/github/scverse/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
[cell cycle reconstruction]: https://github.com/scverse/scanpy_usage/tree/master/170529_images
[dpt example 1]: https://nbviewer.jupyter.org/github/scverse/scanpy_usage/blob/master/170502_paul15/paul15.ipynb
[dpt example 2]: https://nbviewer.jupyter.org/github/scverse/scanpy_usage/blob/master/170501_moignard15/moignard15.ipynb
[myeloid differentiation]: https://nbviewer.jupyter.org/github/scverse/scanpy_usage/blob/master/170430_krumsiek11/krumsiek11.ipynb
[planaria]: https://nbviewer.jupyter.org/github/theislab/paga/blob/master/planaria/planaria.ipynb
[scanpy in r]: https://theislab.github.io/scanpy-in-R/
[seurat to anndata]: https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb
[toggleswitch]: https://nbviewer.jupyter.org/github/scverse/scanpy_usage/blob/master/170430_krumsiek11/toggleswitch.ipynb
