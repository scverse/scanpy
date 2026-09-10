# Normalization: full tutorial audit

Normalization passes the independent numerical checks. All 26,539,065 nonzero values agree within float32 accuracy. The audit found no calculation error in the tutorial path.

The audit uses baseline `a6f1a2d2`, notebook cells 22–23, and checkpoints `03.h5ad` and `04.h5ad`. The checkpoints reside in `/home/fdr/scanpy-tutorial-audit/checkpoints`. They contain 17,041 cells and 23,427 genes. [The script](../../scripts/tutorial_audit/check_04.py) records the calculations and package versions in [the evidence](evidence/04.json).

The [normalization source](../../src/scanpy/preprocessing/_normalization.py) uses every retained gene for each cell's denominator. The default `target_sum=None` selects the median positive row total across both samples: **5,853 counts**. The tutorial does not request separate sample targets or `exclude_highly_expressed=True`.

For raw count `c_ij`, library size `L_i`, and common target `T`, the independent float64 calculation is:

```text
L_i = sum_j(c_ij)
T = median(L_i[L_i > 0]) = 5853
normalized_ij = c_ij * T / L_i
logged_ij = ln(1 + normalized_ij)
```

The [log1p source](../../src/scanpy/preprocessing/_simple.py) uses NumPy's natural logarithm. The saved metadata equals `{'base': None}`. A complete Scanpy replay matches the checkpoint bit for bit.

The input, normalized values, logged values, and saved count layer use float32. All stored nonzero input counts are positive integers. Row totals range from 500 to 150,731, within the exact integer range of float32. Scanpy stores `L_i / T` as float32 and divides counts by that factor. Intermediate rounding explains the small difference from the float64 reference.

| Check | Maximum absolute error | Maximum relative error |
| --- | ---: | ---: |
| Every normalized nonzero value | 3.79e-4 | 1.11e-7 |
| Every logged nonzero value | 5.61e-7 | 1.98e-7 |
| Every normalized row total versus 5,853 | 4.80e-4 | 8.19e-8 |

Normalized values differ by at most one float32 representable step from the rounded reference. Logged values differ by at most two steps. Logged root mean square error is 3.68e-8. These errors have negligible size relative to the sample-target differences measured next.

Cell order, gene order, observation annotations, and gene annotations match. The count layer exactly preserves the input values and sparse structure. All 372,680,442 implicit zeros remain zero. There are no zero-library cells. Source review shows that the median excludes zero totals. The [division helper](../../src/scanpy/_utils/__init__.py) substitutes one for zero divisors, so zero rows remain zero. This dataset does not exercise that branch.

Stored QC totals exceed current row totals by up to 66 counts. Their median difference is zero. Normalization correctly recomputes denominators from the current matrix. This audit also uses current row totals for every association.

| Sample | Cells | Median library size | Median multiplier `5853 / L_i` | Separate target versus common target |
| --- | ---: | ---: | ---: | ---: |
| s1d1 | 8,713 | 5,643 | 1.0372 | −3.59% |
| s1d3 | 8,328 | 6,004.5 | 0.9748 | +2.59% |

Separate sample targets change each sample's normalized scale. Across all nonzero entries, their mean logged changes are −0.01881 and +0.01358, respectively. Maximum absolute changes are 0.03653 and 0.02555. These are parameter-choice effects, not numerical errors. The common target follows the tutorial and assigns equal values to equal count proportions across samples.

Normalization removes differences in prelog row totals. It does not remove all associations with library size or sample identity.

| Association with current library size | Pooled Spearman | s1d1 Spearman | s1d3 Spearman |
| --- | ---: | ---: | ---: |
| Sum of logged expression per cell | 0.195 | 0.312 | 0.067 |
| Number of detected genes | 0.556 | 0.663 | 0.437 |

A single observed count maps to logged values from 0.03810 to 2.54207 across these cells. Multiplicative scaling preserves zero patterns, while the nonlinear logarithm does not preserve equal row totals. These measurements demonstrate remaining depth dependence. They do not identify a technical cause or separate depth, cell type, and sample composition.

Total-count scaling treats each library size as a multiplicative exposure. It represents relative expression and cannot recover absolute transcript abundance. More counts in dominant genes increase the denominator and reduce the relative values of other genes. Biological changes in total RNA or composition can therefore affect the normalized result.

The optional exclusion rule flags a gene globally if its count exceeds 5% of any cell's total. On these data, 36 genes qualify. They include hemoglobin, immunoglobulin, mitochondrial, S100A8, S100A9, and MALAT1 genes. Their combined median count fractions are 25.35% in s1d1 and 16.66% in s1d3. Their 95th-percentile fractions exceed 91% in both samples. This concentration shows why composition matters for the denominator.

An additional audit replay enables exclusion and compares every nonzero value with independent arithmetic. The remaining positive totals have median 3,887. The replay agrees within float32 accuracy. Exclusion changes denominators but still scales all genes. Only the nonexcluded genes sum to the target. The tutorial keeps exclusion disabled, so this sensitivity check does not describe its saved output.

The audit does not establish that normalization removes batch effects or selects an optimal scientific model. It establishes arithmetic correctness and quantifies the remaining associations and alternative-target effects.

Reproduce the numerical checks and configured Ruff checks from the repository root:

```sh
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_04.py --output=/home/fdr/scanpy-tutorial-audit
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_04.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_04.py
```
