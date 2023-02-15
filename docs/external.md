# External API

```{eval-rst}
.. module:: scanpy.external
```

```{eval-rst}
.. currentmodule:: scanpy.external
```

```{note}
More tools that integrate well with scanpy and anndata can be found on the {doc}`ecosystem page <../ecosystem>`.
```

Import Scanpy's wrappers to external tools as:

```
import scanpy.external as sce
```

If you'd like to include a tool here, consider making a pull request ({doc}`instructions <../dev/external-tools>`).
If the tool already uses `scanpy` or `anndata`, it may fit better in the {doc}`ecosystem page <../ecosystem>`.

## Preprocessing: PP

### Data integration

```{eval-rst}
.. autosummary::
   :toctree: generated/

   pp.bbknn
   pp.harmony_integrate
   pp.mnn_correct
   pp.scanorama_integrate

```

### Sample demultiplexing, Doublet detection

```{eval-rst}
.. autosummary::
   :toctree: generated/

   pp.scrublet
   pp.scrublet_simulate_doublets
   pl.scrublet_score_distribution
   pp.hashsolo
```

### Imputation

Note that the fundamental limitations of imputation are still under [debate](https://github.com/scverse/scanpy/issues/189).

```{eval-rst}
.. autosummary::
   :toctree: generated/

   pp.dca
   pp.magic

```

## Tools: TL

### Embeddings

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.phate
   tl.palantir
   tl.trimap
   tl.sam
```

### Clustering and trajectory inference

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.phenograph
   tl.harmony_timeseries
   tl.wishbone
   tl.palantir
   tl.palantir_results
```

### Gene scores, Cell cycle

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.sandbag
   tl.cyclone

```

## Plotting: PL

```{eval-rst}
.. autosummary::
   :toctree: generated/

   pl.phate
   pl.trimap
   pl.sam
   pl.wishbone_marker_trajectory
```

## Exporting

```{eval-rst}
.. autosummary::
   :toctree: generated/

   exporting.spring_project
   exporting.cellbrowser
```
