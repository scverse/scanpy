# Tools: `tl`

```{eval-rst}
.. module:: scanpy.tl
```

```{eval-rst}
.. currentmodule:: scanpy
```

Any transformation of the data matrix that is not *preprocessing*. In contrast to a *preprocessing* function, a *tool* usually adds an easily interpretable annotation to the data matrix, which can then be visualized with a corresponding plotting function.

## Embeddings

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.pca
   tl.tsne
   tl.umap
   tl.draw_graph
   tl.diffmap
```

Compute densities on embeddings.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.embedding_density
```

## Clustering and trajectory inference

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.leiden
   tl.dendrogram
   tl.dpt
   tl.paga
```

(data-integration)=

## Data integration

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.ingest
```

## Marker genes

| Question | Use |
| --- | --- |
| Which genes distinguish a cluster or cell type from the others? | {mod}`scanpy.tl.markers` |
| Which genes change between conditions (e.g. treated vs. control) within a cell type? | {func}`scanpy.tl.de.deseq2`, see {doc}`/tutorials/basics/differential-expression` |

Rank genes that distinguish groups of cells (e.g. clusters), for example to annotate them.
These fast tests treat every cell as an independent observation,
so their p-values are only useful to order genes {cite:p}`Squair2021`.
To plot the results, pass them to e.g. {func}`scanpy.pl.rank_genes_groups_dotplot` as `results=`.

```{eval-rst}
.. module:: scanpy.tl.markers
.. currentmodule:: scanpy
```

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.markers.wilcoxon
   tl.markers.ttest
   tl.markers.logreg
   tl.marker_gene_overlap
```

## Differential expression

Test for expression changes between conditions (e.g. treated vs. control),
using samples (e.g. donors) rather than cells as replicates.

```{eval-rst}
.. module:: scanpy.tl.de
.. currentmodule:: scanpy
```

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.de.deseq2
```

## Gene scores, Cell cycle

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.score_genes
   tl.score_genes_cell_cycle
```
