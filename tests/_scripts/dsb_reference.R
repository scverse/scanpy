# Reference values for tests/test_dsb.py from the dsb R package (2.0.1,
# mclust 6.1.3, limma >= 3.66, which centers the covariates of
# removeBatchEffect). Uses a subset of the example data shipped with dsb
# (CC0 / public domain).
#   Rscript tests/_scripts/dsb_reference.R tests/_data/dsb
suppressPackageStartupMessages(library(dsb))
stopifnot(packageVersion("limma") >= "3.66.0")
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- dsb::cells_citeseq_mtx[, 1:100]
empty <- dsb::empty_drop_citeseq_mtx[, 1:500]
iso <- rownames(cells)[67:70]
w <- function(m, name) {
  con <- gzfile(file.path(out, paste0(name, ".csv.gz")), "w")
  write.csv(m, con, row.names = FALSE)
  close(con)
}
w(t(as.matrix(cells)), "cells")
w(t(as.matrix(empty)), "empty")
quiet <- function(expr) { invisible(capture.output(r <- expr)); r }
r <- quiet(DSBNormalizeProtein(cells, empty, isotype.control.name.vec = iso,
                               return.stats = TRUE))
w(t(r$dsb_normalized_matrix), "default")
w(r$technical_stats, "default_technical_stats")
r <- quiet(DSBNormalizeProtein(cells, empty, isotype.control.name.vec = iso,
                               define.pseudocount = TRUE, pseudocount.use = 5,
                               scale.factor = "mean.subtract", quantile.clipping = TRUE))
w(t(r), "meansub_clip")
r <- quiet(ModelNegativeADTnorm(cells, isotype.control.name.vec = iso))
w(t(r), "negative")
# without isotype controls, step 2 depends on the covariate centering
r <- quiet(DSBNormalizeProtein(cells, empty, use.isotype.control = FALSE))
w(t(r), "no_isotypes")
r <- quiet(ModelNegativeADTnorm(cells, use.isotype.control = FALSE))
w(t(r), "negative_no_isotypes")
