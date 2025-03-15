# Usage Principles

Import Scanpy as:

```
import scanpy as sc
```

## Workflow

The typical workflow consists of subsequent calls of data analysis tools
in `sc.tl`, e.g.:

```
sc.tl.umap(adata, **tool_params)  # embed a neighborhood graph of the data using UMAP
```

where `adata` is an {class}`~anndata.AnnData` object.
Each of these calls adds annotation to an expression matrix *X*,
which stores *n_obs* observations (cells) of *n_vars* variables (genes).
For each tool, there typically is an associated plotting function in `sc.pl`:

```
sc.pl.umap(adata, **plotting_params)
```

If you pass `show=False`, a {class}`~matplotlib.axes.Axes` instance is returned
and you have all of matplotlib's detailed configuration possibilities.

To facilitate writing memory-efficient pipelines, by default,
Scanpy tools operate *inplace* on `adata` and return `None` –
this also allows to easily transition to [out-of-memory pipelines].
If you want to return a copy of the {class}`~anndata.AnnData` object
and leave the passed `adata` unchanged, pass `copy=True` or `inplace=False`.

## AnnData

Scanpy is based on {mod}`anndata`, which provides the {class}`~anndata.AnnData` class.

```{image} https://falexwolf.de/img/scanpy/anndata.svg
:width: 300px
```

At the most basic level, an {class}`~anndata.AnnData` object `adata` stores
a data matrix `adata.X`, annotation of observations
`adata.obs` and variables `adata.var` as `pd.DataFrame` and unstructured
annotation `adata.uns` as `dict`. Names of observations and
variables can be accessed via `adata.obs_names` and `adata.var_names`,
respectively. {class}`~anndata.AnnData` objects can be sliced like
dataframes, for example, `adata_subset = adata[:, list_of_gene_names]`.
For more, see this [blog post].

To read a data file to an {class}`~anndata.AnnData` object, call:

```
adata = sc.read(filename)
```

to initialize an {class}`~anndata.AnnData` object. Possibly add further annotation using, e.g., `pd.read_csv`:

```
import pandas as pd
anno = pd.read_csv(filename_sample_annotation)
adata.obs['cell_groups'] = anno['cell_groups']  # categorical annotation of type pandas.Categorical
adata.obs['time'] = anno['time']                # numerical annotation of type float
# alternatively, you could also set the whole dataframe
# adata.obs = anno
```

To write, use:

```
adata.write_h5ad(filename)
adata.write_zarr(filename)
adata.write_csvs(filename)
adata.write_loom(filename)
```

[blog post]: https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/
[matplotlib]: https://matplotlib.org/
[out-of-memory pipelines]: https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/
[seaborn]: https://seaborn.pydata.org/
