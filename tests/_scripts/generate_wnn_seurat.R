# Called by generate_wnn_seurat.py; all embeddings are synthetic float32 inputs.
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
stopifnot(as.character(packageVersion("Seurat")) == "5.5.1")
directory <- commandArgs(trailingOnly = TRUE)[[1]]

# Change only neighbor search: exact RANN makes the reference independent of
# Annoy approximation and comparable to the port's default brute-force search.
original_nn_helper <- Seurat:::NNHelper
exact_nn_helper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  original_nn_helper(
    data = data, query = query, k = k, method = "rann", cache.index = cache.index
  )
}
assignInNamespace("NNHelper", exact_nn_helper, ns = "Seurat")

pca <- as.matrix(read.csv(file.path(directory, "pca.csv"), header = FALSE))
apca <- as.matrix(read.csv(file.path(directory, "apca.csv"), header = FALSE))
n <- nrow(pca)
cells <- paste0("cell", seq_len(n))
rownames(pca) <- rownames(apca) <- cells
colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
colnames(apca) <- paste0("apca_", seq_len(ncol(apca)))
counts <- sparseMatrix(
  i = rep(1, n), j = seq_len(n), x = 1, dims = c(2, n),
  dimnames = list(c("gene1", "gene2"), cells)
)
object <- CreateSeuratObject(counts = counts, assay = "RNA")
object[["ADT"]] <- CreateAssayObject(counts = counts)
object[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
object[["apca"]] <- CreateDimReducObject(embeddings = apca, key = "apca_", assay = "ADT")
object <- FindMultiModalNeighbors(
  object, reduction.list = list("pca", "apca"),
  dims.list = list(seq_len(ncol(pca)), seq_len(ncol(apca))),
  modality.weight.name = c("RNA.weight", "ADT.weight"),
  k.nn = 20, knn.range = 40, l2.norm = TRUE, sd.scale = 1,
  smooth = FALSE, prune.SNN = 1 / 15, verbose = FALSE
)

write_result <- function(value, name) {
  write.table(
    value, file.path(directory, paste0(name, ".csv")),
    sep = ",", row.names = FALSE, col.names = FALSE
  )
}
write_result(cbind(object$RNA.weight, object$ADT.weight), "weights")
write_result(Indices(object[["weighted.nn"]]) - 1L, "indices")
write_result(Distances(object[["weighted.nn"]]), "distances")
write_result(as.matrix(object[["wsnn"]]), "snn")
writeLines(c(
  "Synthetic inputs: NumPy default_rng(0), 200 cells, 12 RNA PCs / 8 ADT PCs.",
  paste("Seurat", packageVersion("Seurat"), "; RANN", packageVersion("RANN")),
  R.version.string,
  "FindMultiModalNeighbors with NNHelper redirected to exact RANN search.",
  "k.nn=20; knn.range=40; l2.norm=TRUE; sd.scale=1; smooth=FALSE; prune.SNN=1/15.",
  "Cross-modality constants use the Seurat default (1e-4).",
  "Indices are zero-based; weights, distances and SNN values retain R precision.",
  "Regenerate with tests/_scripts/generate_wnn_seurat.py."
), file.path(directory, "provenance.txt"))
