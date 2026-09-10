# PCA: full tutorial audit

The PCA calculation passes independent numerical checks. This audit found no implementation failure in the executed PCA path. All 17,041 cells contribute to the calculation.

PCA retains 61.03% of the variance across 2,000 selected genes. Hemoglobin expression strongly influences PC1. Mitochondrial fraction has its strongest association with PC5, which the tutorial scatter plots do not show. These associations do not establish a bug or a technical cause.

The audit uses baseline `a6f1a2d2`, checkpoint `05.h5ad`, and the exact original output `06.h5ad`. Both checkpoints reside under `/home/fdr/scanpy-tutorial-audit/checkpoints`. Original notebook cells 28, 30, and 32 completed successfully, as recorded in `logs/06.log`. The audit replay reproduces the saved scores exactly. Production source remains unchanged.

[The audit script](../../scripts/tutorial_audit/check_06.py) records all numerical results and package versions in [the evidence](evidence/06.json). The original plots are `plots/06-cell-30-1.png` and `plots/06-cell-32-1.png` under the same external directory.

The [PCA implementation](../../src/scanpy/preprocessing/_pca/__init__.py) resolves `sc.tl.pca(adata)` as follows:

- The default mask selects `.var['highly_variable']`: exactly 2,000 of 23,427 genes.
- Input is sparse CSR float32 `.X`, after total-count normalization and natural `log1p` transformation.
- PCA centers each gene implicitly. The tutorial does not scale genes to unit variance or regress covariates.
- The default solver is scikit-learn PCA with ARPACK. The default component count is 50.
- The [RNG compatibility wrapper](../../src/scanpy/_utils/random.py) supplies legacy `random_state=0` when the caller omits both RNG arguments. Explicit `rng=None` differs from omission.
- Scores occupy `.obsm['X_pca']`, with shape `(17041, 50)` and dtype float32.
- Loadings occupy `.varm['PCs']`, with shape `(23427, 50)`. Float64 padding stores the float32 solver loadings.
- `.uns['pca']` stores float32 variance and variance ratios. Its parameters record the mask and centering, but omit the solver and RNG.

Cell and gene order match between checkpoints. Expression, counts, and all observation and gene annotations remain unchanged. Every loading outside the selected genes equals zero. The selected-gene hash matches the [HVG audit](05-hvg.md).

The independent calculation uses float64 arithmetic:

```text
X = selected log-normalized expression
mu = column means of X
V = saved selected-gene loadings
T = X @ V - mu @ V
C = (X.T @ X - n * outer(mu, mu)) / (n - 1)
C @ V = V * variance
variance_ratio = variance / trace(C)
```

The independent SciPy `eigh` calculation diagonalizes the 2,000-by-2,000 covariance matrix. It does not call Scanpy or scikit-learn PCA. Sign alignment compares individual loadings. Principal angles compare subspaces without a sign assumption.

| Check | Measured error |
| --- | ---: |
| All reconstructed scores, maximum absolute | 3.13e-5 |
| All reconstructed scores, relative Frobenius norm | 9.34e-7 |
| Loading orthonormality, maximum absolute | 1.17e-6 |
| Score mean, maximum absolute | 2.34e-6 |
| Score variance versus stored variance, maximum relative | 1.11e-6 |
| Score covariance, maximum off-diagonal magnitude | 1.94e-6 |
| Covariance eigenvector residual, maximum relative | 1.75e-5 |
| Explained variance ratio, maximum absolute | 3.42e-9 |
| Independent eigenvalues, maximum relative | 1.88e-6 |
| Independent subspace, maximum principal angle | 1.16e-4 radians |
| Independent sign-aligned loadings, maximum absolute | 6.34e-5 |

These errors are consistent with the float32 solver output. The smallest adjacent relative eigenvalue gap among the first 51 components is 0.169%. The PC50-to-PC51 gap is 0.783%. The subspace comparison accounts for sensitivity near close eigenvalues.

Total selected-gene variance is 141.23958. PC1 explains 19.58%, PC2 explains 10.13%, and PC3 explains 5.74%. Cumulative fractions are 38.71% for four PCs, 48.60% for ten, 53.73% for twenty, and 61.03% for fifty. These fractions use selected-gene variance, not all-gene variance. The plot does not establish that fifty is the optimal component count.

The [default neighbors representation](../../src/scanpy/tools/_utils.py) uses all 50 saved PCs. An exact array comparison confirms this choice. The scores are not whitened, so higher-variance PCs contribute more to aggregate squared Euclidean distance.

Quantitative associations use the original stored QC annotations and all cells. Signs refer to this saved PCA orientation.

| Covariate | Strongest absolute Pearson association across 50 PCs | Spearman at that PC |
| --- | --- | ---: |
| Total counts | PC1: 0.469 | 0.186 |
| `log1p_total_counts` | PC6: 0.477 | 0.404 |
| Mitochondrial count percentage | PC5: 0.500 | 0.517 |
| Hemoglobin count percentage | PC1: 0.918 | 0.653 |
| Detected genes | PC6: 0.564 | 0.621 |
| Sample indicator, `s1d3=1` | PC16: 0.236 | 0.275 |

Normalization does not remove these associations. Library size, composition, and cell state can covary. Correlations alone cannot separate their causes. Sample identity explains 5.55% of PC16 variance and 1.07% of total retained score variance through sample means. This statistic does not measure local sample mixing or prove that batch effects are absent.

PC1 has its largest loadings on HBB, HBA2, HBA1, HBD, and MALAT1. HBB log expression correlates with PC1 at 0.942. HBA2 and HBA1 correlations are 0.925 and 0.922. PC2 emphasizes S100A9, S100A8, LYZ, and FCN1. PC3 emphasizes CD74, IGHM, and B-cell markers. PC4 includes NKG7 and GNLY. These patterns agree with marker programs in the tutorial, but PCA loadings alone do not assign cell labels.

Gene-group contributions retain cross-gene covariance. For a group `G`, the per-PC fraction is `sum(V[G] * (C @ V)[G]) / variance`. For exact eigenvectors, this equals the group's squared loading mass. Weighting these fractions by eigenvalues gives the retained-variance fraction. This attribution is not the effect of deleting genes and fitting PCA again.

| Selected gene group | Genes | Fraction of input variance | Fraction of retained variance |
| --- | ---: | ---: | ---: |
| Detected in at most ten cells | 184 | 0.0707% | 0.000670% |
| Absent from one sample | 26 | 0.00965% | 0.000110% |
| Mitochondrial | 2 | 2.043% | 3.150% |
| Hemoglobin | 7 | 12.101% | 19.463% |
| HVG `means > 3` | 6 | 13.169% | 21.197% |

These groups overlap. The six high-mean genes are MTRNR2L12, HBB, MALAT1, HBA2, HBA1, and MT-ND2. Their large contributions reflect covariance PCA without unit-variance scaling. Seven hemoglobin genes account for 51.02% of PC1 loading mass. The two mitochondrial genes account for 14.45% of PC5 loading mass.

The 184 rare genes contribute at most 0.00616% to any retained PC, with the maximum at PC50. Their unusual HVG selection has negligible aggregate influence on these saved PCs. This result does not establish sensitivity for individual rare cells or downstream neighbors.

The displayed PC1–PC4 plots omit the strongest mitochondrial association at PC5. A useful tutorial addition is a PC5/PC6 view with mitochondrial fraction and library size. The audit does not attribute the separate cluster-annotation failure to PCA.

Reproduction commands, from `/home/fdr/scanpy`:

```sh
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_06.py --output /home/fdr/scanpy-tutorial-audit > /home/fdr/scanpy-tutorial-audit/logs/06-check.log 2>&1
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_06.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_06.py
```

All numerical assertions and both configured Ruff checks pass.
