# Scientific audit of the Scanpy tutorial

This audit examines the scientific calculations in the preprocessing and clustering tutorial.
The baseline is `a6f1a2d2` and `docs/tutorials/basics/clustering.ipynb`.
The supervisor reviews dedicated Codex executors in visible Herdr panes.
The review starts with highly variable genes, differential expression, and UMAP.
Normalization, principal components, clustering, and doublet detection support that review.
Data import and basic QC only prepare the dataset.

| Stage | Notebook cells | Scope |
| --- | --- | --- |
| 01 | 1, 2, 4, 5 | Public data download and import |
| 02 | 9, 10, 12, 14, 16 | Quality metrics and filtering |
| 03 | 18 | Doublet detection |
| 04 | 22, 23 | Normalization and log transformation |
| 05 | 25, 26 | Highly variable genes |
| 06 | 28, 30, 32 | Principal components |
| 07 | 34, 36, 38 | Neighbors and UMAP |
| 08 | 41, 42, 45, 46, 52, 54, 55 | Leiden clustering and QC review |
| 09 | 59, 60, 62, 63 | Marker plots and manual annotation |
| 10 | 67, 69, 72, 73 | Differential expression and result plots |

The table describes execution dependencies, not the order of scientific review.
The runner saves an AnnData checkpoint after each stage and plots after each cell.
Data and large artifacts stay outside the Git repository.
Each scientific report records the executed commands, independent checks, source review,
and limits. The final report will distinguish implementation errors from statistical
assumptions and tutorial portability problems.

Run a stage with the preceding checkpoint present:

```sh
python scripts/tutorial_audit/run_stage.py 1 --output /home/fdr/scanpy-tutorial-audit
```

The runner selects a noninteractive plotting backend and eight workers.
It preserves the notebook calculations, parameter choices, and default random seeds.

Initial import and QC work used internal agents before the user clarified the required
Herdr workflow. Those agents stopped. Their preparation artifacts remain available,
but do not count as the requested scientific review.
