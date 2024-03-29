# This script generates the reference test data "seurat_hvg_v3.csv" (this file is compressed to .gz afterwards) and "seurat_hvg_v3_batch.csv"
# Requires to load the data from this link: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

library(dplyr)
library(Seurat)
library(patchwork)

################################################################################
### FindVariableFeatures (no batch covariate)

# Load the PBMC dataset - load the data from the link above!
# pbmc.data <- Read10X(data.dir = "<INSERT_PATH_TO_DATA_HERE>/filtered_gene_bc_matrices/hg19/")
pbmc.data <- Read10X(data.dir = "/Users/eljas.roellin/Documents/R_stuff/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc <- FindVariableFeatures(pbmc,  mean.function=ExpMean, selection.method = 'vst', nfeatures = 2000)

hvf_info <- HVFInfo(pbmc)

write.csv(hvf_info, "seurat_hvg_v3.csv")

################################################################################
### SelectIntegrationFeatures (with batch covariate)

# introduce dummy "technical covariates"
pbmc.data <- Read10X(data.dir = "/Users/eljas.roellin/Documents/R_stuff/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k_2", min.cells = 3, min.features = 200)
pbmc

pbmc@meta.data["dummy_tech"] = "source_1"
pbmc@meta.data[501:1001, "dummy_tech"] = "source_2"
pbmc@meta.data[1001:1500, "dummy_tech"] = "source_3"
pbmc@meta.data[1501:2000, "dummy_tech"] = "source_4"
pbmc@meta.data[2001:nrow(pbmc@meta.data), "dummy_tech"] = "source_5"

pbmc_list = SplitObject(pbmc, split.by='dummy_tech')

features = SelectIntegrationFeatures(pbmc_list)

write.csv(features, "seurat_hvg_v3_batch.csv")
