# Plotting: `pl`

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

## Generic

Functions with both backends:

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/
   :template: function-dual

   pl.scatter
   pl.heatmap
   pl.dotplot
   pl.tracksplot
   pl.violin
   pl.stacked_violin
   pl.matrixplot
   pl.ranking

```

Helper for computing dot areas, e.g. for {func}`~scanpy.pl.dotplot`:

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   pl.dot_area

```

Legacy (matplotlib) only:

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   pl.clustermap
   pl.dendrogram

```

## Classes

These classes allow fine tuning of visual parameters.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/classes

   pl.DotPlot
   pl.MatrixPlot
   pl.StackedViolin

```

## Preprocessing

Methods for visualizing quality control and results of preprocessing functions.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: function-dual

   pl.highest_expr_genes
   pl.highly_variable_genes
   pl.scrublet_score_distribution

```

## Tools

Methods that extract and visualize tool-specific annotation in an
{class}`~anndata.AnnData` object.  For any method in module `tl`, there is
a method with the same name in `pl`.

### PCA

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

### Embeddings

Functions with both backends:

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: function-dual

   pl.umap
   pl.tsne
   pl.diffmap
   pl.draw_graph
```

Compute densities on embeddings.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: function-dual

   pl.embedding_density
```

Legacy (matplotlib) only:

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.spatial
   pl.embedding
```

### Branching trajectories and pseudotime, clustering

Visualize clusters using one of the embedding methods passing e.g. `color='leiden'`.

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

Visualize hierarchical clustering results as a heatmap.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.correlation_matrix
```

### Marker genes

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

### Simulations

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.sim

```
