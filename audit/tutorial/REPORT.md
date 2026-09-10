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

The HVG result supports the implementation of the stated method on this dataset.
It does not establish that each selected gene is a reliable biological marker.
The pooled-sample method selects a different list, with 1,445 genes in common.
See the [HVG review](05-hvg.md) for the formulas, source paths, tolerances, and limits.

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
