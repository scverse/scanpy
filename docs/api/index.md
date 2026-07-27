# API

Import Scanpy as:

```
import scanpy as sc
```

```{note}
Additional functionality is available across the broader [scverse ecosystem](https://scverse.org/packages/#ecosystem), with some tools wrapped in the {mod}`scanpy.external` module.
```

(array-support)=
## Array type support

Different APIs have different levels of support for array types,
and this page lists the supported array types for each function
(⚡ indicates support of the type as chunk in a dask {class}`~dask.array.Array`):

```{eval-rst}
.. array-support:: all
```

```{toctree}
:maxdepth: 2
:hidden:

preprocessing
tools
plotting
io
get
queries
metrics
experimental
classes
settings
datasets
deprecated
```
