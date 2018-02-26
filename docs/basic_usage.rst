Basic Usage
-----------

Import the Scanpy API as::

    import scanpy.api as sc

Workflow
^^^^^^^^

The typical workflow consists of subsequent calls of data analysis tools
in ``sc.tl``, e.g.::

    sc.tl.louvain(adata, **tool_params)  # cluster cells using Louvain clustering

where ``adata`` is an :class:`~scanpy.api.AnnData` object. Each of these calls adds annotation to an expression matrix *X*, which stores *n_obs* observations of *n_vars* gene expression variables. For each tool, there is at least one associated plotting function in ``sc.pl``, which retrieves and plots the added annotation::

    sc.pl.louvain(adata, **plotting_params)

If you pass ``show=False``, a `matplotlib.Axes <https://matplotlib.org/api/axes_api.html>`_ instance is returned and you have all of matplotlib's detailed configuration possibilities.

To facilitate writing memory-efficient pipelines, by default, Scanpy tools operate *inplace* on ``adata`` and return ``None`` - this also allows to easily transition to `out-of-memory pipelines <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_. If you want to return a copy of the :class:`~scanpy.api.AnnData` object and leave the passed ``adata`` unchanged, pass ``copy=True``.

    
AnnData
^^^^^^^

Scanpy is based on `anndata <http://anndata.readthedocs.io>`_, which provides
the :class:`~scanpy.api.AnnData` class.

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

At the most basic level, an :class:`~scanpy.api.AnnData` object ``adata`` stores
a data matrix (``adata.X``), dataframe-like annotation of observations
(``adata.obs``) and variables (``adata.var``) and unstructured dict-like
annotation (``adata.uns``). Values can be retrieved and appended via
``adata.obs['key1']`` and ``adata.var['key2']``. Names of observations and
variables can be accessed via ``adata.obs_names`` and ``adata.var_names``,
respectively. :class:`~scanpy.api.AnnData` objects can be sliced like
dataframes, for example, ``adata_subset = adata[:, list_of_gene_names]``.
For more, see this `blog post <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_.
         
To read a data file to an :class:`~scanpy.api.AnnData` object, call::

    adata = sc.read(filename)

to initialize an :class:`~scanpy.api.AnnData` object. Possibly add further annotation using, e.g., ``pd.read_csv``::

    import pandas as pd 
    anno = pd.read_csv(filename_sample_annotation)
    adata.obs['cell_groups'] = anno['cell_groups']  # categorical annotation of type pandas.Categorical
    adata.obs['time'] = anno['time']                # numerical annotation of type float
    # alternatively, you could also set the whole dataframe
    # adata.obs = anno

To write, use::

    adata.write(filename)
    adata.write_csvs(filename)
    adata.write_loom(filename)    


.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/
