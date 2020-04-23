library(reticulate)
use_python("/usr/local/bin/python3")
library(Seurat)

sc <- import("scanpy", convert=FALSE)

ad <- sc$datasets$pbmc3k()
ad$X <- ad$X$toarray()

ad$var_names_make_unique()

sc$pp$highly_variable_genes_seurat_v3(ad, n_top_genes=50, lowess_frac=0.15)


X <- py_to_r(ad$X$T)
bc <- py_to_r(ad$obs_names$tolist())
colnames(X) <- bc
genes <- py_to_r(ad$var_names$tolist())
rownames(X) <- genes

sr <- CreateSeuratObject(counts = X)
ad_hvg <- py_to_r(ad$var["highly_variable"])
ad_hvg <- names(ad_hvg[ad_hvg])

sr <- FindVariableFeatures(object = sr, selection.method = "vst", nfeatures = 50)
sr_hvg <- sr@assays$RNA@var.features

print(all(sr_hvg==ad_hvg))

hvg_info <- as.data.frame(HVFInfo(sr))
hvg_info <- cbind(hvg_info, rownames(hvg_info) %in% sr_hvg)
colnames(hvg_info) <- c("means", "variances", "variances_norm", "highly_variable")

write.table(hvg_info, "seurat_hvg_v3.csv")
