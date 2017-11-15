Basic Usage Principles
----------------------

Import the Scanpy API as::

    import scanpy.api as sc

The typical workflow consists of subsequent calls of data analysis tools
of the form::

    sc.tl.louvain(adata, **params)

where ``adata`` is an :class:`~scanpy.api.AnnData` object and ``params`` are optional parameters. Each of these calls adds annotation to an expression matrix *X*, which stores *n* *d*-dimensional gene expression measurements. To facilitate writing memory-efficient pipelines, by default, Scanpy tools operate *inplace* on ``adata`` and return ``None``. If you want to copy the :class:`~scanpy.api.AnnData` object, pass the ``copy`` argument::

    adata_copy = sc.tl.louvain(adata, copy=True, **params)

    
AnnData objects
^^^^^^^^^^^^^^^

An :class:`~scanpy.api.AnnData` object ``adata`` stores a data matrix
(``adata.data``), dataframe-like sample (``adata.smp``) and variable
(``adata.var``) annotation and unstructured dict-like annotation
(``adata.uns``).

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

Values can be retrieved and appended via ``adata.smp['key1']`` and
``adata.var['key2']``. Sample and variable names can be accessed via
``adata.smp_names`` and ``adata.var_names``,
respectively. :class:`~scanpy.api.AnnData` objects can be sliced like
dataframes, for example, ``adata_subset = adata[:, list_of_gene_names]``. The AnnData
class is similar to R's ExpressionSet [Huber15]_. For more, see :class:`~scanpy.api.AnnData`.
    

Reading and writing data files and AnnData objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To read, call::

    adata = sc.read(filename)

to initialize an :class:`~scanpy.api.AnnData` object. Possibly add further annotation using, e.g., ``pd.read_csv``::

    import pandas as pd 
    anno = pd.read_csv(filename_sample_annotation)
    adata.smp['cell_groups'] = anno['cell_groups']  # categorical annotation of type str or int
    adata.smp['time'] = anno['time']                # numerical annotation of type float

or set a whole dataframe::

    adata.smp = pd.read_csv(filename_sample_annotation)

To write, use::

    sc.write(filename, adata)

to save the :class:`~scanpy.api.AnnData` as a collection of data arrays to a file in a platform and language-independent way. Reading foresees filenames with extensions *h5*, *xlsx*, *mtx*, *txt*, *csv* and others. Writing foresees writing *h5*, *csv* and *txt*. For more, see :func:`~scanpy.api.write` and :func:`~scanpy.api.read`.

Plotting
^^^^^^^^

For each tool, there is an associated plotting function::

    sc.pl.tool(adata)

that retrieves and plots annotation in ``adata`` that has been added by ``sc.tl.tool(adata)``. Scanpy's plotting module can be viewed similar to Seaborn_: an extension of matplotlib_ that allows visualizing operations on AnnData objects with one-line commands. Detailed configuration has to be done via matplotlib functions, which is easy as Scanpy's plotting functions accept and return a ``Matplotlib.Axes`` object.

.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/
