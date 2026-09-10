# Highly variable genes: full tutorial audit

The calculation selects the expected 2,000 genes. Independent calculations match every selection flag, batch count, and intersection flag. This audit found no implementation bug in the executed HVG calculation.

The main scientific limitation is the interpretation of batch selection. Only 628 selected genes qualify in both samples. The other 1,372 qualify in one sample. The selection also includes 184 genes detected in ten cells or fewer.

## Execution and scope

The source baseline is `a6f1a2d2`. The working HEAD at execution was `76c35e894c09235cf4e3a117d125e6319038131d`. Its changes from the baseline contain audit preparation only. Production HVG code matches the baseline.

The original dataset contains 17,125 cells. The full post-QC checkpoint contains 17,041 cells and 23,427 genes. This audit uses every post-QC cell and gene, without subsampling.

Original notebook cells 18, 22, 23, 25, and 26 completed successfully. No prerequisite failure or harness change occurred. Scrublet annotates doublets but does not remove them at this stage. Its successful execution does not establish the scientific correctness of Scrublet. The incomplete `check_03.py` was preliminary context only.

The notebook SHA-256 is `327c31fa8a746b15332620008c63b53470db3eefde3ec65f82bac047cc395851`. Package versions and numerical results are in [evidence/05.json](evidence/05.json).

Exact commands, from `/home/fdr/scanpy`:

```sh
mkdir -p /home/fdr/scanpy-tutorial-audit/logs
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/run_stage.py 3 --output /home/fdr/scanpy-tutorial-audit > /home/fdr/scanpy-tutorial-audit/logs/03-execution.log 2>&1
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/run_stage.py 4 --output /home/fdr/scanpy-tutorial-audit > /home/fdr/scanpy-tutorial-audit/logs/04-execution.log 2>&1
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/run_stage.py 5 --output /home/fdr/scanpy-tutorial-audit > /home/fdr/scanpy-tutorial-audit/logs/05-execution.log 2>&1
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_05.py --output /home/fdr/scanpy-tutorial-audit > /home/fdr/scanpy-tutorial-audit/logs/05-check.log 2>&1
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_05.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_05.py
```

## Calculation traced

Cell 25 calls `sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")`.

The default ScanpyV1 preset resolves `flavor` to `seurat`. It reads `.X`, uses 20 bins, writes `.var`, and keeps all genes in the matrix. `n_top_genes=2000` disables the mean and dispersion cutoffs. `span` and `check_values` do not apply to this flavor.

Cell 23 divides each cell's counts by its library size and multiplies by the global median library size, 5,853. It then applies the natural `log1p` transformation. HVG calculation applies `expm1` to recover normalized expression. Thus the variance calculation uses normalized expression before the logarithm, including implicit zeros.

For each sample and expressed gene, the calculation uses:

```text
mu = sum(x) / n
variance = sum((x - mu)**2) / (n - 1)
m = log1p(mu)
d = log(variance / mu)
z = (d - mean(d within the mean bin)) / sample_sd(d within the mean bin)
```

Scanpy computes variance through second moments with the same `n - 1` correction. The independent script uses centered squared deviations, including the contribution from implicit zeros.

`pd.cut(m, bins=20)` creates equal-width, right-closed bins. Its lower endpoint extends by 0.1% of the range to include the minimum. Bin standard deviations use `ddof=1`. A singleton bin gets normalized dispersion 1 through the source's special handling. The actual data contain no nonfinite normalized dispersions among expressed genes.

Within each sample, Scanpy selects every gene at or above the 2,000th largest normalized dispersion. Equal scores at this boundary can produce more than 2,000 per-sample HVGs. Neither sample has a boundary tie in this run. The script records minimum ranks for inspection. This flavor does not calculate or combine median ranks.

Genes absent from a sample do not participate in that sample's bins. Scanpy restores them with zero metrics and a false selection flag. It then averages each metric equally across samples, including these restored zeros. These means are averages of sample log means, not pooled expression means. Pandas excludes NaNs from these averages, but this run has no such expressed-gene values.

The final ordering uses decreasing HVG batch count, then decreasing mean normalized dispersion. Scanpy marks the first 2,000 genes and restores the original gene order. Exact ties in both sort keys inherit the grouped gene order in this implementation. Such ties do not determine the final boundary here.

Source: [`_highly_variable_genes.py`](../../src/scanpy/preprocessing/_highly_variable_genes.py), especially the single-batch calculation, bin normalization, `_subset_genes`, and batch aggregation. The Seurat primary documentation describes this binned dispersion method and its equal-width bins. Its current `mean.var.plot` interface uses cutoffs, so Scanpy's explicit top-2,000 extension is not an exact parameter match. [Seurat FindVariableFeatures](https://satijalab.org/seurat/reference/findvariablefeatures)

Seurat's integration feature selection breaks batch-count ties with median ranks. The executed Scanpy flavor instead uses mean normalized dispersion. [Seurat SelectIntegrationFeatures](https://satijalab.org/seurat/reference/selectintegrationfeatures)

## Numerical evidence

The expected results use only NumPy, SciPy, and pandas calculations. Public Scanpy calls supply additional per-sample comparison results, without supplying any expected metric.

| Quantity | s1d1 | s1d3 |
| --- | ---: | ---: |
| Cells | 8,713 | 8,328 |
| Expressed genes | 23,183 | 23,124 |
| Absent genes | 244 | 303 |
| Per-sample HVGs | 2,000 | 2,000 |
| Normalized dispersion cutoff | 0.9857832353 | 1.0285669035 |
| Selection disagreements | 0 | 0 |
| Largest normalized dispersion error | 9.375e-8 | 9.578e-8 |

Against the final checkpoint, maximum absolute errors are zero for `means`, `5.059e-8` for `dispersions`, and `4.316e-7` for `dispersions_norm`. The assertions use `rtol=1e-5, atol=1e-6`. Boolean flags and batch counts require exact equality.

The sparse input and `expm1` values use float32. Sparse CSR second-moment products use input precision before float64 accumulation. Scanpy stores final normalized dispersion in float32 after selection. Means and log dispersions remain float64. These precision choices explain the small numerical differences without changing selection.

A separate calculation starts from counts and performs normalization entirely in float64. Recovered float32 expression differs by at most 0.002148 expression units. All per-sample and final selection flags still match. Final normalized dispersion differs by at most `4.931e-7`.

The final boundary selects `TRERF1` at 0.7939602809 and excludes `AC139768.1` at 0.7939050194. Both qualify in one sample. Their score gap, approximately `5.526e-5`, exceeds the observed numerical errors.

The audit also checks that HVG calculation preserves `.X`, counts, cell order, gene order, observation annotations, and existing gene annotations.

## Scientifically consequential method choices

- Batch selection favors shared variability but does not require it. There are 628 shared HVGs and 2,744 genes that qualify in one sample. The final list contains all 628 shared HVGs and 1,372 of the others. It includes 26 genes absent from one sample.
- Batch selection changes the features substantially. A pooled calculation shares 1,445 genes with the tutorial selection. Selection by mean normalized dispersion alone shares 1,793. These comparisons change only the independent audit calculation, not the notebook.
- Mean cutoffs do not apply. Among selected genes, 999 have reported mean at most 0.0125 and six have mean at least 3. Detection ranges from three to 17,032 cells. The 184 genes detected in ten cells or fewer can represent rare biology or unstable sparse observations. This audit does not distinguish those explanations.
- Most genes occupy the first mean bin: 21,622 in s1d1 and 21,814 in s1d3. High-expression bins contain few genes. The bin structure therefore gives uneven reference populations for dispersion normalization.
- Singleton scores are conventional values, not measured within-bin z-scores. In s1d1, `HBB`, `MALAT1`, and `MT-ND2` receive score 1 and pass its cutoff. In s1d3, singleton genes receive 1 but fail its higher cutoff. The final list includes `HBB`, `MALAT1`, `MT-ND2`, and `HBA2`, but excludes `MT-ATP6`.
- Each sample receives equal weight despite different cell counts. Equal-weight and cell-weighted averages of sample log means differ by at most 0.01470. Neither quantity is the pooled log mean.

The next notebook PCA call automatically uses the HVG mask. Therefore these choices determine which features enter PCA and subsequent neighborhood calculations. HVG selection does not remove batch effects from expression values. It also supplies no significance test or false-discovery guarantee. The measured feature differences do not establish a particular change in clusters or biological conclusions.

## Plot inspection and limits

I inspected `/home/fdr/scanpy-tutorial-audit/plots/05-cell-26-1.png`. Both panels contain the expected gene clouds and selected-gene overlay. The normalized panel shows selected and unselected genes with overlapping scores, consistent with the priority for batch counts.

The axis labels omit the transformations and equal-sample averaging. The x-axis represents averaged `log1p` means. The unnormalized y-axis represents averaged log variance-to-mean ratios. Thus negative values there are valid log dispersions. The figure also hides batch counts and sparse detection, so it cannot establish that every selected gene varies in both samples.

No confirmed HVG implementation bug emerged. The method choices above are observed behavior. Sensitivity to bin counts, cell resampling, alternative HVG flavors, and later cluster assignments remains untested. General edge cases involving zero variance, empty batches, or exact final ties are outside this dataset audit.

The reusable script is [`check_05.py`](../../scripts/tutorial_audit/check_05.py). Small evidence is in [05.json](evidence/05.json). Complete per-gene metrics and minimum per-sample ranks are external:

- `/home/fdr/scanpy-tutorial-audit/metrics/05-independent-gene-metrics.csv`
- `/home/fdr/scanpy-tutorial-audit/metrics/05-s1d1-gene-metrics.csv`
- `/home/fdr/scanpy-tutorial-audit/metrics/05-s1d3-gene-metrics.csv`
- `/home/fdr/scanpy-tutorial-audit/metrics/05-selected-genes.txt`

The selected-gene list follows the checkpoint gene order. Its SHA-256 is `53e3d2d500f15644602a5cff71e537c3d65a764563730ebe51a3e9b76fb867ca`.
