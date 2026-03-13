# Preprocessing: PP

```{eval-rst}
.. module:: scanpy.external.pp
.. currentmodule:: scanpy.external
```

Previously found here, but now part of scanpy’s main API:
- {func}`scanpy.pp.harmony_integrate`
- {func}`scanpy.pp.scrublet`
- {func}`scanpy.pp.scrublet_simulate_doublets`

(external-data-integration)=

## Data integration

```{eval-rst}
.. autosummary::
   :toctree: ../generated/

   pp.bbknn
   pp.mnn_correct
   pp.scanorama_integrate
```

## Sample demultiplexing

```{eval-rst}
.. autosummary::
   :toctree: ../generated/

   pp.hashsolo
```

## Imputation

Note that the fundamental limitations of imputation are still under [debate](https://github.com/scverse/scanpy/issues/189).

```{eval-rst}
.. autosummary::
   :toctree: ../generated/

   pp.dca
   pp.magic
```
