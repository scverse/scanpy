# Neighbor graph and UMAP audit

The executed graph agrees with the UMAP membership formula. This audit found no numerical implementation bug in graph construction or embedding output. It found a distance-storage documentation mismatch and substantial distortion of individual neighborhoods in two dimensions.

The approximate graph recovers 98.40% of exact neighbors in the evaluation sample. The embedding retains only 14.95% of exact neighbors. These measurements describe different operations. Neither establishes biological correctness.

## Scope and reproduction

The source baseline is `a6f1a2d2`. The audit reads the full saved tutorial checkpoints `06.h5ad` and `07.h5ad`. Both contain 17,041 cells and 23,427 genes. The original stage log records successful notebook cells 34, 36, and 38. The neighbor and UMAP calls took 20.669 and 4.862 seconds. This audit did not rerun that stage.

Artifacts reside under `/home/fdr/scanpy-tutorial-audit`:

- Inputs: `checkpoints/06.h5ad` and `checkpoints/07.h5ad`.
- Original execution: `logs/07.log` and `metrics/07-execution.json`.
- Inspected plot: `plots/07-cell-38-1.png`.
- Audit output: `logs/07-check.log`.

The [reusable audit script](../../scripts/tutorial_audit/check_07.py) writes [numerical evidence](evidence/07.json). The evidence includes package versions, query indices, cell names, and a row-name hash. The script uses no Scanpy neighbor helper or UMAP formula function for its independent calculations.

```sh
NUMBA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 /tmp/scanpy-de-audit-env/bin/python scripts/tutorial_audit/check_07.py > /home/fdr/scanpy-tutorial-audit/logs/07-check.log
/tmp/scanpy-audit-tools/bin/ruff check scripts/tutorial_audit/check_07.py
/tmp/scanpy-audit-tools/bin/ruff format --check scripts/tutorial_audit/check_07.py
```

No production code or unit tests changed.

## Representation, search, and self treatment

Cell 34 calls `sc.pp.neighbors(adata)`. The [representation selector](../../src/scanpy/tools/_utils.py) uses the existing `X_pca` because the gene count exceeds `settings.N_PCS`. With `n_pcs=None`, it uses all existing components. This checkpoint contains 50 PCs. The function does not choose a component count from the variance plot. It does not whiten the PCs or balance samples.

The defaults are `n_neighbors=15`, `metric="euclidean"`, `method="umap"`, and legacy `random_state=0`. At this cell count, the [backend selector](../../src/scanpy/neighbors/_common.py) chooses PyNNDescent, with 12 trees and 14 iterations. The execution harness sets eight jobs.

The saved distance matrix contains 16 stored entries per row: one zero self entry and 15 positive nonself distances. The [neighbor extraction](../../src/scanpy/neighbors/_common.py) retains self plus the nearest 14 nonself entries for fuzzy graph construction. Self membership is zero. The fifteenth nonself distance stays in storage but does not contribute directly to that row's fuzzy memberships.

This conflicts with the [public distance documentation](../../src/scanpy/neighbors/__init__.py), which promises `n_neighbors-1` nonzero entries per row. This is a documentation/storage-contract mismatch in the executed backend. The graph itself uses the intended 14 nonself neighbors. Downstream readers must not assume that every stored distance contributes an outgoing fuzzy edge.

## Independent actual-data results

Exact Euclidean search uses 1,024 deterministic query cells against **all 17,041 cells**. NumPy seed 703 selects 512 cells without replacement from each sample. Float64 SciPy `cdist` supplies all candidate distances. Self exclusion uses row identity. Evaluation uses 14 nonself neighbors and stable distance sorting.

The reported means give equal weight to the two sample strata. They are sample estimates, not exhaustive recall or embedding scores. Exact search and embedding rank calculations use the same queries. Distance accuracy, graph structure, formula reconstruction, and mixing use the full graph.

| Measurement | Result |
|---|---:|
| Mean approximate-neighbor recall | 0.983956 |
| Queries with all 14 exact neighbors | 866 / 1,024 |
| Minimum query recall | 0.571429 |
| 5th percentile query recall | 0.928571 |
| Maximum excess in approximate neighbor radius | 1.9706% |
| True ranks of substituted neighbors | 15–22 |
| Maximum absolute distance error, all stored entries | 1.1901e-6 |
| Maximum relative nonself distance error | 1.6691e-7 |
| Mean exact-neighbor overlap in UMAP | 0.149484 |
| Mean query trustworthiness | 0.959372 |
| Mean query continuity | 0.982313 |

Distance errors match float32 precision. Neighbor omissions are approximation behavior, not incorrect distance arithmetic. Recall differs slightly by sample: 0.985491 for `s1d1` and 0.982422 for `s1d3`.

For each query, trustworthiness penalizes UMAP neighbors whose PCA rank exceeds 14. Its normalization is `14 * (2*N - 3*14 - 1)`. Continuity reverses the two spaces. The script averages these per-query contributions with `N=17041`. It does not compute ranks within a reduced reference dataset.

High trustworthiness does not mean exact local preservation. Its rank penalty is relative to the full population. Mean exact overlap is about 2.09 of 14 neighbors. Median overlap is one neighbor. The 5th percentile is zero.

## Fuzzy membership and objective

The [Scanpy wrapper](../../src/scanpy/neighbors/_connectivity.py) passes the extracted neighbors to `umap.fuzzy_simplicial_set`. Installed umap-learn 0.5.12 uses `local_connectivity=1` and `set_op_mix_ratio=1`.

For this dataset, every nonself distance is positive. The local offset `rho_i` is the nearest nonself distance. The scale `sigma_i` solves the following equation, subject to the reference minimum scale:

```text
sum over 14 nonself neighbors exp(-max(0, d_ij-rho_i)/sigma_i) = log2(15)
p_ij = exp(-max(0, d_ij-rho_i)/sigma_i)
w_ij = p_ij + p_ji - p_ij*p_ji
```

The independent script solves this equation in float64 and reconstructs the full fuzzy union. Its support matches exactly, and its maximum weight difference is 3.55e-6. The small difference reflects the reference binary-search tolerance and float32 arithmetic. This agrees with the fuzzy union in the [primary UMAP paper](https://arxiv.org/abs/1802.03426) and the [authors' explanation](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html).

The saved graph has 383,780 directed entries, exact symmetry, zero diagonal, and one connected component. All weights are finite, between 6.95e-9 and 1. Degrees range from 14 to 338, with median 19. Symmetrization adds incoming neighbors, so degrees need not equal 14.

Cell 36 calls `sc.tl.umap(adata)`. The [UMAP wrapper](../../src/scanpy/tools/_umap.py) uses two dimensions, spectral initialization, 200 epochs, `min_dist=0.5`, and `spread=1`. It sets `alpha=1`, `gamma=1`, and `negative_sample_rate=5`. Saved curve parameters are `a=0.5830300203` and `b=1.3341669924`.

The low-dimensional membership is `q_ij = 1/(1+a*||y_i-y_j||^(2*b))`. The method targets fuzzy cross-entropy through sampled attractive edges and negative samples. The installed `umap/layouts.py` implements these updates with gradient clipping. It does not evaluate every pair's full loss at each step. The [primary method](https://arxiv.org/abs/1802.03426) describes the objective and sampling procedure.

The installed embedding routine removes weights below `max_weight/200` from its private graph copy before optimization. This removes 480 directed entries, but the graph still has one connected component. The saved connectivity matrix remains intact.

The [legacy RNG decorator](../../src/scanpy/_utils/random.py) supplies seed 0 when the tutorial omits both RNG arguments. Thus the apparent `rng=None` signature does not imply unseeded tutorial execution. The wrapper disables parallel embedding optimization for this legacy seed. This audit checks the saved realization, not coordinate reproducibility across versions or platforms.

## Sample mixing, identity, and interpretation

All observation names, observation metadata, variable names, and PCA values match checkpoint 06 in the same order. Observation names are unique. All 17,041 embedding rows are distinct and finite, with shape `(17041, 2)`.

The graph contains 8,713 `s1d1` cells and 8,328 `s1d3` cells. Cross-sample edges account for 30.74% of edges and 29.06% of graph weight. Mean per-cell cross-sample weight fractions are 29.39% and 29.04%. Random labels with these sample sizes imply about 48.87% and 51.13%, respectively.

This comparison is descriptive, not a batch-effect significance test. Samples mix less than random labels. Cell-type composition, biological differences, and technical effects can all contribute. These scores cannot separate those causes. Neither neighbors nor UMAP applies batch correction.

The inspected plot shows overlapping sample colors, local sample enrichment, narrow bridges, and apparent islands. The saved graph is nevertheless connected. Visible gaps do not establish disconnected biological populations. Plot distances between islands, their area, and their orientation do not measure biological distance or abundance. UMAP's local distance rescaling also prevents direct density interpretation, as the [authors explain](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html).

The measured loss of exact neighbors limits conclusions about individual adjacent points. Marker expression and analyses in the original representation remain necessary for biological interpretation. No embedding-quality score establishes cell identity, a lineage, or a valid cluster boundary.
