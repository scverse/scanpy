# Scientific audit of the standard Scanpy tutorial

This report examines calculations on the full public tutorial dataset.
The source baseline is `a6f1a2d23496e97a01860029f7bd57b0bcce6e27`.
It is an audit of this development checkout, not every Scanpy release.

The tutorial starts with 17,125 cells from two samples.
Its basic filters leave 17,041 cells and 23,427 genes.
Data import and file handling are preparation, outside the scientific review.

## Results

The audit is in progress. Completed reviews are listed here.

| Calculation | Independent result | Scientific interpretation |
| --- | --- | --- |
| Highly variable genes | Exact agreement for all 2,000 selected genes. Maximum normalized-dispersion error: `4.316e-7`. | Only 628 qualify in both samples. The list includes 184 genes detected in ten cells or fewer. |
| Differential expression | All 398,259 gene/group comparisons pass. Default p-values match exactly; maximum BH error is `2.221e-16`. | Tie correction changes the number of pairs with adjusted p < 0.05 from 80,557 to 165,700. |
| Neighbors and UMAP | Full fuzzy-graph support matches; maximum weight error is `3.55e-6`. Exact-neighbor recall is 98.40% on 1,024 query cells. | The two-dimensional embedding retains 14.95% of exact neighbors for those queries. |

The HVG result supports the implementation of the stated method on this dataset.
It does not establish that each selected gene is a reliable biological marker.
The pooled-sample method selects a different list, with 1,445 genes in common.
See the [HVG review](05-hvg.md) for the formulas, source paths, tolerances, and limits.

The original tutorial DGE call completes with all genes and groups, including both plots.
Cluster 7 has NK-associated markers consistent with the notebook interpretation.
Changing the fold-change definition reverses 25,882 signs, but most changes involve
at least one effect of one log2 unit or less. Only four reversals exceed one log2 unit
in both definitions and pass the default adjusted-p threshold.
The [DGE review](10-dge.md) gives those actual examples and their expression fractions.

The DGE calculation describes clusters derived from the same expression data.
Its per-cell p-values do not establish population-level treatment or disease effects.
The detailed review separates cluster-selection bias and biological replication
from the correctness of the numerical formulas.

The graph is connected, although the UMAP plot shows apparent islands.
Its trustworthiness score is 0.95937, but exact neighbor overlap is much lower.
Leiden uses the original graph, so projection distortion does not directly change clusters.
The [UMAP review](07-umap.md) separates full-graph checks from sampled embedding metrics.
It also records one documentation mismatch: stored distance rows contain 15 nonself
entries, while the documented count is 14. Fuzzy memberships correctly use 14.

## Method and evidence

Each scientific review uses a dedicated Codex executor in a visible Herdr pane.
The supervisor checks the numerical evidence and source review before it accepts the work.
The executor pane closes only after that review.

The runner executes original notebook cells and preserves their scientific parameters.
It saves a checkpoint after each stage and a plot after each plotting cell.
Independent audit scripts compare these outputs with separate numerical calculations.
No production Scanpy code changes form part of this audit.

- [Execution stages and reproduction](README.md)
- [Source and dataset provenance](evidence/provenance.json)
- [Pinned package environment](evidence/environment.txt)

Large local artifacts are in `/home/fdr/scanpy-tutorial-audit`.
The dataset DOI is [10.6084/m9.figshare.22716739.v1](https://doi.org/10.6084/m9.figshare.22716739.v1).
The source notebook is [clustering.ipynb](../../docs/tutorials/basics/clustering.ipynb).
