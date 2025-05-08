## Plotting: `pl`

```{eval-rst}
.. module:: scanpy.pl
```

```{eval-rst}
.. currentmodule:: scanpy
```

The plotting module {mod}`scanpy.pl` largely parallels the `tl.*` and a few of the `pp.*` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

See {doc}`/tutorials/plotting/core` for an overview of how to use these functions.

```{note}
See the {ref}`settings` section for all important plotting configurations.
```

(pl-generic)=

### Generic

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   pl.scatter
   pl.heatmap
   pl.dotplot
   pl.tracksplot
   pl.violin
   pl.stacked_violin
   pl.matrixplot
   pl.clustermap
   pl.ranking
   pl.dendrogram

```

### Classes

These classes allow fine tuning of visual parameters.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/classes

    pl.DotPlot
    pl.MatrixPlot
    pl.StackedViolin

```

### Preprocessing

Methods for visualizing quality control and results of preprocessing functions.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.highest_expr_genes
   pl.filter_genes_dispersion
   pl.highly_variable_genes
   pl.scrublet_score_distribution

```

### Tools

Methods that extract and visualize tool-specific annotation in an
{class}`~anndata.AnnData` object.  For any method in module `tl`, there is
a method with the same name in `pl`.

#### PCA

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.pca
   pl.pca_loadings
   pl.pca_variance_ratio
   pl.pca_overview
```

(pl-embeddings)=

#### Embeddings

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.tsne
   pl.umap
   pl.diffmap
   pl.draw_graph
   pl.spatial
   pl.embedding
```

Compute densities on embeddings.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.embedding_density
```

#### Branching trajectories and pseudotime, clustering

Visualize clusters using one of the embedding methods passing `color='louvain'`.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.dpt_groups_pseudotime
   pl.dpt_timeseries
   pl.paga
   pl.paga_path
   pl.paga_compare
```

#### Marker genes

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.rank_genes_groups
   pl.rank_genes_groups_violin
   pl.rank_genes_groups_stacked_violin
   pl.rank_genes_groups_heatmap
   pl.rank_genes_groups_dotplot
   pl.rank_genes_groups_matrixplot
   pl.rank_genes_groups_tracksplot
```

#### Simulations

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.sim

```
