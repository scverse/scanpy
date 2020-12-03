library(reticulate)
library(Seurat)

sc <- import("scanpy.api", convert=FALSE)

ad <- sc$datasets$pbmc68k_reduced()
ad$X <- ad$raw$X$toarray()

ad$var_names_make_unique()

sr <- Convert(ad, to="seurat")

sc$pp$normalize_per_cell(ad, counts_per_cell_after=1e4)
sc$pp$log1p(ad)

sr <- NormalizeData(object=sr, normalization.method="LogNormalize", scale.factor=1e4)

sc$pp$highly_variable_genes(ad, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5)

ad_hvg <- py_to_r(ad$var["highly_variable"])
ad_hvg <- names(ad_hvg[ad_hvg])

sr <- FindVariableGenes(object=sr, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)

sr_hvg <- sr@var.genes

print(all(sr_hvg==ad_hvg))

hvg_info <- as.data.frame(sr@hvg.info)
hvg_info <- cbind(hvg_info, rownames(hvg_info) %in% sr_hvg)
colnames(hvg_info) <- c("means", "dispersions", "dispersions_norm", "highly_variable")
hvg_info <- hvg_info[rownames(sr@data),]

write.table(hvg_info, "seurat_hvg.csv")
