# Quickstart

Every scanpy analysis works on an {class}`~anndata.AnnData` object: {mod}`scanpy.pp` preprocesses it, {mod}`scanpy.tl` computes on it and {mod}`scanpy.pl` plots it.
The example below clusters 3k PBMCs and finds the marker genes of each cluster.

```python
import scanpy as sc

adata = sc.datasets.pbmc3k()

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.tl.rank_genes_groups(adata, "leiden")

sc.pl.umap(adata, color="leiden")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=3)
```

Most functions store their results in `adata` instead of returning them, e.g. the clusters in `adata.obs["leiden"]`.
Leiden clustering needs igraph, which `pip install 'scanpy[leiden]'` installs, see {doc}`installation`.

{doc}`tutorials/basics/clustering` walks through the same steps with quality control and doublet detection, {doc}`usage-principles` explains the conventions behind the API and {doc}`api/index` lists every function.
