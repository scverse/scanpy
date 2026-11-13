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

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.rank_genes_groups
   tl.filter_rank_genes_groups
   tl.marker_gene_overlap
```

## Gene scores, Cell cycle

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.score_genes
   tl.score_genes_cell_cycle
```

## Output types

Functions with several outputs return a `…Result` if passed `out=True`,
and take a `…Out` mapping as `out` to write each output to a place of choice.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.CellCycleResult
   tl.CellCycleOut
   tl.DptResult
   tl.DptOut
```
