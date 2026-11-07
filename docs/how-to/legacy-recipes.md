# Migrating from legacy preprocessing recipes

`sc.pp.recipe_seurat`, `sc.pp.recipe_zheng17`, and `sc.pp.recipe_weinreb17` have been removed.
Compose individual preprocessing functions to make the transformations and their parameters explicit.
{func}`scanpy.experimental.pp.recipe_pearson_residuals` is unaffected.

The examples below start from non-logarithmized counts and reproduce the final expression matrices of the former Seurat and Zheng recipes with their default arguments and the `ScanpyV1` preset.
They subset into a new `AnnData` object; start with `adata = adata.copy()` if you also need to preserve the input object.
Unlike the recipes, calling {func}`~scanpy.pp.log1p` directly on `.X` also records `.uns["log1p"]`.
Nondefault arguments, especially `log=False`, require separate consideration of the original transformation order.

## Seurat

```python
import scanpy as sc

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
hvg = sc.pp.highly_variable_genes(
    adata, flavor="seurat", min_mean=0.0125, max_mean=3, min_disp=0.5, out=True
)
adata = adata[:, hvg["highly_variable"]].copy()
sc.pp.scale(adata, zero_center=True, max_value=10)
```

## Zheng et al. (2017)

```python
import scanpy as sc

sc.pp.filter_genes(adata, min_counts=1)
normalized = sc.pp.normalize_total(adata, out=True)
adata.X = normalized.normalized
adata.obs["n_counts_all"] = normalized.norm_factor
sc.pp.log1p(adata)
hvg = sc.pp.highly_variable_genes(
    adata, flavor="cell_ranger", n_top_genes=1000, out=True
)
adata = adata[:, hvg["highly_variable"]].copy()
sc.pp.scale(adata, zero_center=True)
```

The normalization target is the median nonzero cell total before highly variable gene selection.
`n_counts_all` retains the former recipe's normalization factors (cell totals divided by that target), rather than raw counts.
Gene selection uses the log-normalized matrix, matching the former implementation.
With the default `log=True`, the recipe's second normalization was discarded when it restored its saved log-transformed layer; do not add a second normalization to this workflow.
To inspect gene dispersion, call {func}`~scanpy.pl.highly_variable_genes` on `hvg` before subsetting (`log=True` for Zheng, `log=False` for the default Seurat recipe).

## Weinreb et al. (2017)

There is no numerically equivalent replacement using only the standard preprocessing functions.
The former recipe optionally logged **before** normalization, excluded genes exceeding 5% of a cell's total, and scaled by the mean total over all genes.
It selected genes by mean and coefficient of variation, then computed PCA using population standard deviations plus `0.0001` for standardization.
In contrast, {func}`~scanpy.pp.scale` uses sample standard deviations, and {func}`~scanpy.pp.highly_variable_genes` uses different selection criteria.
For new analyses, explicitly choose {func}`~scanpy.pp.normalize_total`, {func}`~scanpy.pp.log1p`, {func}`~scanpy.pp.highly_variable_genes`, {func}`~scanpy.pp.scale`, and {func}`~scanpy.pp.pca` parameters appropriate to the data; this is a change of analysis, not an equivalent substitution.
To reproduce an existing Weinreb analysis exactly, retain the Scanpy version and parameters used for that analysis.
