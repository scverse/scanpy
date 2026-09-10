# Leiden clustering and biological annotation audit

Leiden and marker dotplots reproduce correctly on the saved tutorial data. This audit found no algorithm or plot calculation bug in these paths.

The hardcoded annotation map fails on this run. It swaps the apparent B-cell and erythroid groups and leaves 1,743 cells unannotated. Retained cells with high mitochondrial percentages and predicted doublets also affect cluster interpretation.

## Scope and reproduction

The baseline is `a6f1a2d2`. Inputs contain 17,041 cells and 23,427 genes. The audit covers notebook cells 41–63 and the cluster-7 interpretation in cells 70–73. It uses checkpoints 07–09 and reads the saved cluster-7 DGE results from checkpoint 10.

The [audit script](../../scripts/tutorial_audit/check_08.py) writes [numerical evidence](evidence/08.json). The evidence includes all cluster sizes, sample counts, QC summaries, marker means, expression fractions, and resolution cross-tables.

```sh
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MPLBACKEND=Agg /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_08.py > /home/fdr/scanpy-tutorial-audit/logs/08-check.log 2>&1
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_08.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_08.py
```

The calculation and both Ruff checks passed. The environment uses Scanpy `1.14.0.dev6+ga6f1a2d23` and igraph `1.0.0`. No production code or unit tests changed.

## Confirmed annotation failure

Cell 62 maps numeric labels without reference to the current marker profiles. The saved coarse partition contains five labels, not the four labels in that map.

| Coarse cluster | Cells | Assigned annotation | Observed marker profile |
| --- | ---: | --- | --- |
| 0 | 9,763 | Lymphocytes | Broad T/NK signals, plus mixed subgroups |
| 1 | 1,507 | Monocytes | FCN1 and CD14 enrichment |
| 2 | 2,269 | Erythroid | MS4A1 and PAX5 enrichment supports a B-cell interpretation |
| 3 | 1,759 | B Cells | Strong HBA1/HBB with little MS4A1/PAX5 supports an erythroid interpretation |
| 4 | 1,743 | Missing | Strong erythroid markers with more MKI67 and GYPA |

These observations support broad interpretations. They do not establish biological ground truth or resolve every cell subtype.

The independent marker summaries make the swap clear. Means refer to log1p-normalized expression across all cells in each group.

| Marker | Cluster 2 mean / positive % | Cluster 3 mean / positive % | Cluster 4 mean / positive % |
| --- | ---: | ---: | ---: |
| MS4A1 | 1.496 / 81.8 | 0.004 / 1.0 | 0.033 / 4.8 |
| PAX5 | 0.737 / 62.8 | 0.003 / 0.7 | 0.016 / 3.3 |
| HBA1 | 1.176 / 84.8 | 6.142 / 100.0 | 5.438 / 100.0 |
| HBB | 2.805 / 99.9 | 8.204 / 100.0 | 7.431 / 100.0 |
| GYPA | 0.016 / 1.9 | 0.208 / 41.0 | 1.459 / 88.4 |
| MKI67 | 0.016 / 1.9 | 0.006 / 1.8 | 0.739 / 61.7 |

The two conflicting labels cover 4,028 cells (23.64%). Another 1,743 cells (10.23%) receive missing values. These are tutorial annotation failures, not Leiden membership or row-alignment failures.

The broad `Lymphocytes` label also hides heterogeneity. It includes all 84 cells of resolution-0.5 cluster 3 and all 106 cells of cluster 10. Their plotted profiles emphasize the notebook's dendritic and plasma markers, respectively. These small groups need separate review before biological use.

Numeric cluster labels do not encode cell identity. A portable tutorial must derive its annotations from the current profiles and cover every observed label. This audit does not replace the map with new biological ground truth.

## Leiden dispatch, graph, and objective

The [Leiden implementation](../../src/scanpy/tools/_leiden.py) explicitly dispatches these calls to igraph. It uses an undirected, weighted graph and sets `objective_function="modularity"`. Direct igraph defaults to CPM, so an unqualified direct call is not a valid parity comparison.

The saved connectivity matrix is exactly symmetric. It contains 383,780 nonzero entries, no self edges, and finite positive weights from approximately `6.95e-9` to 1.

The baseline [graph conversion](../../src/scanpy/_utils/__init__.py) creates an undirected edge for every nonzero matrix entry. Thus, the graph contains 191,890 pairs of parallel edges. Each pair represents one symmetric neighbor relationship.

This duplication doubles every edge weight in aggregate. It does not change generalized modularity for a fixed partition. The independent calculation uses the original symmetric matrix:

```text
S = sum_ij A_ij
K_c = sum_{i in c} sum_j A_ij
Q_gamma = sum_{i,j in same cluster} A_ij / S - gamma * sum_c (K_c / S)^2
```

The script also constructs a graph from only the upper triangle. Its fixed-partition objective agrees with the duplicated graph to floating-point precision. No selective edge weighting error appears here.

However, a fresh fit on the simple graph can follow a different stochastic path. Adjusted Rand indices against the saved fits are 0.81122, 1.0, 0.99851, and 0.96856 for the table order below. Equal objectives for fixed partitions do not imply identical optimization trajectories.

| Saved key | Resolution | Iterations | Clusters | Size range | Independent objective |
| --- | ---: | ---: | ---: | ---: | ---: |
| `leiden` | 1 | 2 | 25 | 61–2,798 | 0.8552573692601333 |
| `leiden_res_0.02` | 0.02 | -1 | 5 | 1,507–9,763 | 0.9886894814657039 |
| `leiden_res_0.50` | 0.5 | -1 | 17 | 61–3,951 | 0.9043706119400584 |
| `leiden_res_2.00` | 2 | -1 | 36 | 61–1,207 | 0.8013290868421141 |

Every independent objective agrees with the stored objective within `4e-16`. These scores use different resolutions and cannot rank biological accuracy across rows.

For all four calls, direct igraph and fresh Scanpy results match the saved labels exactly. Repeated direct calls also match exactly. The direct comparison reproduces the graph order, weights, objective, resolution, iterations, beta, and legacy RNG stream independently of Scanpy helpers.

The default seed is legacy `random_state=0`, through the [RNG wrapper](../../src/scanpy/_utils/random.py). An unrelated Python RNG with seed zero is not the same random stream. Exact reproduction here applies to this environment and graph order, not every package version or graph representation.

Categories follow natural numeric order, including labels above 9. Cell names remain unique and identical across checkpoints 07–09. Existing observation columns remain unchanged, and the graph is identical between checkpoints 07 and 08. No category ordering or row-alignment discrepancy appears.

## Iteration and resolution sensitivity

Cell 41 requests two iterations. Cell 52 omits `n_iterations`, so Scanpy passes `-1`. The initial result and later results therefore differ in both resolution and stopping rule.

At resolution 1 and seed 0, convergence produces 26 clusters instead of 25. The objective increases from 0.85525737 to 0.85611949. The adjusted Rand index is 0.97785.

Negative iterations stop after a stable iteration. They do not certify a global optimum. This distinction follows the [igraph API description](https://python.igraph.org/en/latest/api/igraph.community.html#_community_leiden) and matters for Scanpy's phrase “optimal clustering.”

Three additional seeds show the observed sensitivity:

| Resolution / iterations | Cluster counts for seeds 1, 2, 3 | Adjusted Rand index range against seed 0 |
| --- | --- | ---: |
| 1 / 2 | 27, 26, 25 | 0.79878–0.86998 |
| 0.02 / -1 | 5, 4, 5 | 0.95545–1.00000 |
| 0.5 / -1 | 17, 17, 18 | 0.92483–0.99276 |
| 2 / -1 | 36, 36, 35 | 0.83721–0.89831 |

Even the coarse cluster count changes with seed. The four-cluster solution at seed 2 has an objective only 0.00004194 below the saved five-cluster solution.

Resolution results are not strictly nested. For example, resolution-0.5 cluster 6 combines 310 cells from coarse cluster 0 and 71 from coarse cluster 1. A larger resolution does not guarantee a biological hierarchy.

## Marker dotplot calculations

I inspected the saved plots `09-cell-60-1.png` and `09-cell-63-1.png` under `/home/fdr/scanpy-tutorial-audit/plots`. Their visible profiles agree with the independent summaries.

The [dotplot implementation](../../src/scanpy/plotting/legacy/_dotplot.py) uses checkpoint 09 `.X`. This checkpoint has no `.raw`, and the notebook selects no layer. Colors use arithmetic means of log1p-normalized expression, including zeros. Dot sizes use the fraction with expression greater than zero.

For each gene, `standard_scale="var"` subtracts the minimum group mean and divides by the range across groups. Constant columns become zero. This transformation changes color values but preserves expression fractions.

| Comparison | Coarse plot maximum error | Resolution-0.5 plot maximum error |
| --- | ---: | ---: |
| Unscaled means | `2.06e-7` | `5.21e-7` |
| Scaled means | `9.19e-8` | `1.31e-7` |
| Positive fractions | 0 | 0 |

The small mean errors reflect float32 aggregation versus independent float64 calculations. Row order and all gene columns match, including repeated genes in the marker dictionary.

Colors show relative enrichment across the displayed groups. They do not compare absolute abundance across genes or across plots with different groupings. For example, HBB is positive in at least 99.9% of every coarse group, despite pale colors outside the erythroid groups.

This widespread signal does not prove erythroid identity in every group. Its cause remains unresolved by this audit. The notebook also lists negative markers, but the plot does not interpret their signs. Cell-type decisions still require positive and negative evidence together.

## Sample composition and retained QC signals

The saved data contain 8,713 `s1d1` cells and 8,328 `s1d3` cells. Every reported cluster contains cells from both samples. Some groups have substantial sample imbalance. This observation cannot distinguish technical effects from biological composition.

The table gives the resolution-0.5 partition used for marker analysis. The JSON contains equivalent summaries for all four partitions. Mitochondrial values refer to `pct_counts_mt`. The 20% threshold is descriptive, not a proposed universal filter.

| Cluster | Cells | s1d1 / s1d3 | Predicted doublets % | Mitochondrial mean / median % | Cells above 20% mitochondrial % |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 3951 | 1483 / 2468 | 0.13 | 6.86 / 6.12 | 1.9 |
| 1 | 1440 | 980 / 460 | 2.78 | 8.20 / 7.33 | 2.2 |
| 2 | 974 | 648 / 326 | 3.18 | 17.13 / 9.63 | 32.6 |
| 3 | 84 | 62 / 22 | 0.00 | 5.75 / 4.54 | 3.6 |
| 4 | 3045 | 1817 / 1228 | 0.23 | 33.52 / 27.95 | 77.3 |
| 5 | 1694 | 606 / 1088 | 0.06 | 0.41 / 0.33 | 0.0 |
| 6 | 381 | 310 / 71 | 0.26 | 35.67 / 36.46 | 80.3 |
| 7 | 459 | 166 / 293 | 0.44 | 8.62 / 7.71 | 1.7 |
| 8 | 386 | 109 / 277 | 0.78 | 7.45 / 5.62 | 6.2 |
| 9 | 1057 | 642 / 415 | 1.23 | 6.34 / 3.05 | 6.4 |
| 10 | 106 | 80 / 26 | 0.00 | 13.26 / 6.36 | 23.6 |
| 11 | 911 | 567 / 344 | 2.52 | 9.63 / 8.45 | 4.4 |
| 12 | 1300 | 467 / 833 | 0.08 | 8.51 / 7.46 | 2.9 |
| 13 | 624 | 404 / 220 | 9.94 | 2.83 / 1.60 | 1.3 |
| 14 | 61 | 42 / 19 | 37.70 | 4.87 / 4.42 | 1.6 |
| 15 | 503 | 273 / 230 | 0.20 | 11.24 / 6.41 | 13.9 |
| 16 | 65 | 57 / 8 | 0.00 | 1.14 / 0.49 | 0.0 |

All 213 predicted doublets remain (1.25% overall). Cluster 14 contains 23/61 predicted doublets and shows both erythroid and monocyte markers. Its FCN1 and CD14 positive fractions are 98.4% and 82.0%. Its GYPA fraction is 83.6%.

This combination requires review before a pure cell-type label. Predicted doublets are model calls, not known doublets. Neither the marker mixture nor these calls alone establish the mechanism.

Cluster 4 contains 3,045 cells with median mitochondrial percentage 27.95%. Cluster 6 has median 36.46%. Their QC profiles support caution before a distinct biological subtype interpretation.

The initial partition also contains strong QC-associated groups. Initial cluster 10 has median mitochondrial percentage 59.78%. Initial clusters 17 and 19 contain 32.88% and 37.70% predicted doublets. These numeric labels refer to `leiden`, not `leiden_res_0.50`.

## Hardcoded DGE cluster 7

The cluster-7 example is supported in this saved run. Its first five DGE genes are NKG7, GNLY, KLRD1, PRF1, and CST7. Independent means for NKG7 and GNLY are 3.420 and 3.534. Their positive fractions are 100% and 94.8%.

These profiles support the notebook's tentative NK interpretation. They do not establish cluster purity. TRBC2 is positive in 62.1% of these cells, so a categorical claim needs broader marker review.

The hardcoded number remains a portability assumption. It refers specifically to resolution-0.5 cluster 7. Initial `leiden` cluster 7 instead has 1,894 cells and median mitochondrial percentage 26.26%. The tutorial must identify the partition and current marker evidence whenever it assigns biological meaning to a numeric label.
