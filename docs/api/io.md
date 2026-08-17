(reading)=
(reading-and-writing)=

# Reading and Writing

```{eval-rst}
.. currentmodule:: scanpy
```

Write {class}`~anndata.AnnData` objects using its {doc}`writing <anndata:api>` methods

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   io.write
```

```{note}
For reading annotation use {ref}`pandas.read_… <pandas:io>`
and add it to your {class}`~anndata.AnnData` object. The following read functions are
intended for the numeric data in the data matrix `X`.
```

Read common file formats using

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   io.read
```

Read 10x formatted hdf5 files and directories containing `.mtx` files using

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   io.read_10x_h5
   io.read_10x_mtx
```

Read other formats using functions borrowed from {mod}`anndata`

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   io.read_h5ad
   io.read_csv
   io.read_excel
   io.read_hdf
   io.read_mtx
   io.read_text
   io.read_umi_tools

```
