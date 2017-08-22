Basic Usage Principles
----------------------

The typical workflow consists of subsequent calls of data analysis tools
of the form::

    sc.tl.louvain(adata, **params)

where ``adata`` is an ``AnnData`` object and ``params`` are optional parameters. Each of these calls adds annotation to an expression matrix *X*, which stores *n* *d*-dimensional gene expression measurements. To facilitate writing memory-efficient pipelines, by default, Scanpy tools operate *inplace* on ``adata`` and return ``None``. If you want to copy the ``AnnData`` object, pass the ``copy`` argument::

    adata_copy = sc.tl.louvain(adata, copy=True, **params)

Reading and writing data files and AnnData objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One usually calls::

    adata = sc.read(filename)

to initialize an AnnData object, possibly adds further annotation using, e.g., ``np.genfromtxt`` or ``pd.read_csv``::

    annotation = pd.read_csv(filename_annotation)
    adata.smp['cell_groups'] = annotation['cell_groups']  # categorical annotation of type str or int
    adata.smp['time'] = annotation['time']                # numerical annotation of type float

and uses::

    sc.write(filename, adata)

to save the ``adata`` as a collection of data arrays to a file in a platform and language-independent way. Reading foresees filenames with extensions *h5*, *xlsx*, *mtx*, *txt*, *csv* and others. Writing foresees writing *h5*, *csv* and *txt*.

AnnData objects
^^^^^^^^^^^^^^^

An ``AnnData`` object stores an array-like data matrix as ``adata.X``, dataframe-like sample annotation as ``adata.smp``, dataframe-like variable annotation as ``adata.var`` and additional unstructured dict-like annotation as ``adata.add``. While ``adata.add`` is a conventional dictionary, ``adata.smp`` and ``adata.var`` are instances of a low-level Pandas dataframe-like class.

Values can be retrieved and appended via ``adata.smp[key]`` and ``adata.var[key]``. Sample and variable names can be accessed via ``adata.smp_names`` and ``adata.var_names``, respectively. AnnData objects can be sliced like Pandas dataframes, for example, ``adata = adata[:, list_of_gene_names]``. The AnnData class is similar to R's ExpressionSet [Huber15]_ the latter though is not implemented for sparse data.

Plotting
^^^^^^^^

For each tool, there is an associated plotting function::

    sc.pl.tool(adata)

that retrieves and plots annotation in ``adata`` that has been added by ``sc.tl.tool(adata)``. Scanpy's plotting module can be viewed similar to Seaborn_: an extension of matplotlib_ that allows visualizing operations on AnnData objects with one-line commands. Detailed configuration has to be done via matplotlib functions, which is easy as Scanpy's plotting functions accept and return a ``Matplotlib.Axes`` object.

.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/
