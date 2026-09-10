# Differential expression: full tutorial audit

Every tested gene and group agrees with independent calculations. This audit found no implementation bug in the executed DGE path. The NK-cell description for group 7 is supported. Statistical choices substantially change significance counts and effect sizes.

## Scope and execution

The source baseline is `a6f1a2d2`. Production code and the notebook match that baseline. The audit HEAD was `47b4094c711c09fec20c9ffe0a20f788b96dba0b`.

Original stages 06–10 completed in order from checkpoint05. All original DGE cells, including both plots and the group-7 lookup, passed. There were no notebook runtime failures or parameter substitutions. DGE cell 67 took 4.725 seconds. This audit uses all 17,041 retained cells, all 23,427 genes, and all 17 groups: 398,259 comparisons. Group sizes range from 61 to 3,951 cells. The original dataset contains 17,125 cells before tutorial filtering. No audit subsampling occurred.

The reusable script is [`check_10.py`](../../scripts/tutorial_audit/check_10.py). Full results are in [`10.json`](evidence/10.json). Per-group sensitivity, marker expression, sample composition, and every underflowed gene/group pair are in the adjacent `10-*.csv` files.

Commands, from `/home/fdr/scanpy`:

```sh
for stage in 6 7 8 9 10; do
  NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/run_stage.py "$stage" --output /home/fdr/scanpy-tutorial-audit > "/home/fdr/scanpy-tutorial-audit/logs/$(printf '%02d' "$stage").log" 2>&1
done
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_10.py > /home/fdr/scanpy-tutorial-audit/logs/10-check.log 2>&1
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_10.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_10.py
```

Large checkpoints, execution logs, and original plots remain under `/home/fdr/scanpy-tutorial-audit`. The JSON includes execution metadata and package versions. Notebook SHA-256: `327c31fa8a746b15332620008c63b53470db3eefde3ec65f82bac047cc395851`.

## Exact statistical path

Cell 67 calls:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
```

The ScanpyV1 defaults select `tie_correct=False`, `mean_in_log_space=True`, `mask_var=None`, `reference="rest"`, and Benjamini–Hochberg correction. `.raw` is absent, so the call reads float32 `.X`. The `counts` layer does not enter DGE. `.X` contains natural-log `log1p` expression after library normalization to 5,853 counts. The 2,000-HVG annotation does not restrict the tested gene family.

For each gene, the independent script ranks all cells with SciPy average ranks. It sums ranks over the exact group mask and uses every other cell as the reference:

```text
N = 17041
n = cells in the group
z = (group_rank_sum - n*(N+1)/2) / sqrt(n*(N-n)*(N+1)/12)
p = 2*normal_survival(abs(z))
q[sorted_i] = min(1, min over j>=i of p[sorted_j]*23427/j)
```

There is no continuity correction. This is an unpaired rank-sum calculation with a two-sided normal approximation. SciPy documents the missing tie adjustment in its equivalent `ranksums` approximation. [SciPy ranksums](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ranksums.html)

BH operates separately over all 23,427 genes within each group. It does not adjust across the 17 groups. An additional `n_genes=5` call agrees with full-family BH, which excludes correction after truncation. The default output contains every gene in descending signed-score order. Top genes therefore favor positive scores, not absolute effects or smallest two-sided p-values. Equal scores need no unique gene order.

The default reported effect is:

```text
log2fc = log2((expm1(mean(log1p(expression_group))) + 1e-9)
           / (expm1(mean(log1p(expression_rest))) + 1e-9))
```

This compares shifted geometric-style means. It is not the ratio of arithmetic expression means.

Source trace: [`_RankGenes` and `_build_stats_dataframe`](../../src/scanpy/tools/_rank_genes_groups.py), [`ScanpyV1 defaults`](../../src/scanpy/_settings/presets.py), and [`rank_genes_groups_df`](../../src/scanpy/get/get.py). The previous logistic sign, layer-cache, and t-test findings do not occur in this notebook call.

## Numerical results

The expected calculations use NumPy/SciPy directly, in chunks of 128 genes. They do not call Scanpy rank or aggregation helpers.

| Quantity | Maximum absolute error |
| --- | ---: |
| Default scores | 3.792e-6 |
| Default p-values | 0 |
| Default BH values | 2.221e-16 |
| Default log2 fold changes | 9.537e-7 |
| Tie-adjusted scores | 3.813e-6 |
| Tie-adjusted p-values | 4.441e-16 |
| Tie-adjusted BH values | 8.882e-16 |
| Arithmetic log2 fold changes | 1.007e-6 |

Every comparison passes `rtol=2e-6`. Scores and fold changes also allow `atol=2e-6`. P-values and BH values use zero absolute tolerance. Scanpy stores scores and fold changes as float32. The remaining differences fit storage precision and arithmetic order.

All groups have complete unique gene names, descending scores, correct group/rest counts, and exact `get` field alignment. The arithmetic sensitivity preserves every Wilcoxon score.

No tested score, effect size, p-value, or BH value is NaN or infinite. Default normal tails return zero for 1,112 comparisons. BH preserves those zeros. Tie correction increases both zero counts to 3,205. [`10-underflow.csv`](evidence/10-underflow.csv) identifies every pair and gives a finite log-tail calculation. Default zero p-values have log10 tails from −1,619.147 to −310.333. Some remain representable through a log-tail calculation, so this is also a normal-survival implementation limit. Zero is not an exact probability or proof of biological certainty.

## Sensitivity on the actual data

Tie correction multiplies the null variance by `1 - sum(t**3-t)/(N**3-N)`, where `t` is each tied-value multiplicity. Its median coefficient is 0.040803. Its range is 0.000528–0.999999997. No gene is constant. Sparse expression makes this adjustment substantial.

| Result | Default | Tie correction |
| --- | ---: | ---: |
| Gene/group pairs with BH < 0.05 | 80,557 | 165,700 |
| Group-7 pairs with BH < 0.05 | 2,482 | 6,098 |
| Zero p-values | 1,112 | 3,205 |

Tie correction adds 85,143 significant pairs. Top-20 overlap ranges from 2 to 20 genes per group. Group 7 retains nine. Its adjusted top five are KLRF1, TRDC, SH2D1B, KLRD1, and CLIC3. The adjustment changes gene rankings, not just displayed p-values. It does not fix cluster selection or donor dependence.

The arithmetic sensitivity uses `mean_in_log_space=False` with Wilcoxon. No t-test runs. Across all pairs, 25,882 effect signs reverse. Of these, 25,873 have at least one absolute log2 effect at most one. Only nine exceed one in both definitions. Four of those also pass default BH < 0.05. Thus the overall reversal count mainly describes small or asymmetric effects. It does not describe 25,882 strong reversals.

All four material, significant examples are:

| Group | Gene | Default log2FC | Arithmetic log2FC | Default BH | Expressing: group | Expressing: rest |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 2 | IGHG1 | 1.384 | −1.259 | 0.01382 | 10.57% | 4.54% |
| 3 | IGHM | 1.053 | −1.052 | 9.324e-19 | 91.67% | 21.96% |
| 11 | IGHG4 | 1.133 | −1.013 | 0.007783 | 11.75% | 5.23% |
| 11 | IGHG1 | 1.145 | −1.644 | 0.004469 | 11.42% | 4.52% |

The group can contain more expressing cells while its arithmetic mean remains lower. For group-2 IGHG1, arithmetic means are 0.266 versus 0.636 normalized counts. For group-3 IGHM, they are 1.342 versus 2.783. Different detection frequencies and expression magnitudes explain why the effect definitions can disagree. These are sensitivity findings, not bugs. Group 7 has no reversal above one log2 unit in both definitions.

Exact examples and both expression fractions are in [`10-foldchange-reversals.csv`](evidence/10-foldchange-reversals.csv). Per-group threshold counts are in [`10-sensitivity.csv`](evidence/10-sensitivity.csv).

Per-group median absolute log2-effect differences range from 0.172 to 0.688. The largest difference is 4.478. Among default-significant pairs, positive effects above one log2 unit decrease from 21,558 to 16,317.

For group 7, NKG7 changes from 6.546 to 4.616 log2 units. GNLY changes from 7.379 to 5.503. Both remain strong positive markers. These effect definitions measure different quantities. The disagreement is a scientific sensitivity, not evidence of incorrect arithmetic.

## Group 7 and plots

Group 7 contains 459 cells: 166 from s1d1 and 293 from s1d3. Its original top five are NKG7, GNLY, KLRD1, PRF1, and CST7.

NKG7 appears in 100% of group-7 cells. GNLY, KLRD1, and PRF1 appear in 94.8%, 96.1%, and 95.0%. CD3D and TRAC appear in 13.3% and 10.9%, with negative group/rest effects. MS4A1, CD79A, and CD14 appear in fewer than 3%. This supports an NK-enriched annotation. Primary NK profiling uses the same NK-lineage markers. [Yang et al., 2019](https://www.nature.com/articles/s41467-019-11947-7)

The annotation does not prove purity or a particular NK subtype. TRDC appears in 74.9% of cells, and tie correction places it second. CD3E appears in 37.5%. TRDC is not exclusive to gamma-delta T cells. Primary single-cell work also detects it in NK cells. [Kazer et al., 2020](https://www.nature.com/articles/s41591-020-0799-2)

Additional cell-level or protein evidence is necessary to distinguish admixture, ambient expression, and related lymphocyte populations. This audit makes no stronger identity claim.

Both original plots were inspected. The UMAP panels use the expected five genes and show enrichment around group 7. The dotplot includes all 17 groups and five selected genes per group, with dendrogram reordering. All underlying dot fractions agree exactly with independent expression fractions. Scaled mean colors differ by at most 1.616e-7.

Dot color is the mean log expression, scaled separately from zero to one for each gene. Dot size is the fraction with expression above zero. Neither reports p-values or fold changes. The color scale cannot compare absolute expression between genes. Broad ribosomal and mitochondrial top markers in other groups show why top-rank selection alone does not establish cell identity.

## Scientific interpretation beyond arithmetic

The same expression data determine HVGs, PCA, neighbors, clusters, and cluster-marker tests. Gene reuse creates selection bias. BH arithmetic cannot repair p-values from this circular procedure. These outputs support exploratory cluster descriptions, not a calibrated claim of 5% false discoveries among biological markers. Primary research demonstrates this post-clustering inference problem. [Grabski et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC11282907/)

The rank-sum call also pools individual cells without a sample or donor term. The two sample labels contain 8,713 and 8,328 cells. Their proportions vary substantially across groups. For example, group 16 contains 57 s1d1 cells and eight s1d3 cells. The source dataset includes donor and site variation. [Lance et al., 2022](https://proceedings.mlr.press/v176/lance22a/lance22a.pdf)

Cells from one biological sample do not supply independent donor replication. This matters separately from circular cluster selection. Studies of single-cell DGE show false discoveries when methods ignore biological replication. [Zimmerman et al., 2021](https://www.nature.com/articles/s41467-021-21038-1), [Squair et al., 2021](https://www.nature.com/articles/s41467-021-25960-2)

Donor-aware pseudobulk or mixed models address replication for suitable prespecified contrasts. They do not automatically remove cluster-selection bias. This tutorial supplies no condition contrast or donor-level inferential analysis. Its marker calculation must not become a population-level treatment or disease claim.
