def bbknn(adata, batch_key='batch', copy=False, **kwargs):
    """\
    Batch balanced kNN [Park18]_.

    Batch balanced kNN alters the kNN procedure to identify each
    cell's top neighbours in each batch separately instead of the
    entire cell pool with no accounting for batch. Aligns batches in a
    quick and lightweight manner.

    For use in the scanpy workflow as an alternative to
    :func:`~scanpy.pp.neighbors`.

    .. note::

        This is just a wrapper of :func:`bbknn.bbknn`: up to date docstring, more information
        and bug reports `here <https://github.com/Teichlab/bbknn>`__.

    Params
    ------
    adata : ``AnnData``
        Needs the PCA computed and stored in ``adata.obsm["X_pca"]``.
    batch_key : ``str``, optional (default: "batch")
        ``adata.obs`` column name discriminating between your batches.
    neighbors_within_batch : ``int``, optional (default: 3)
        How many top neighbours to report for each batch; total number of neighbours
        will be this number times the number of batches.
    n_pcs : ``int``, optional (default: 50)
        How many principal components to use in the analysis.
    trim : ``int`` or ``None``, optional (default: ``None``)
        Trim the neighbours of each cell to these many top connectivities. May help with
        population independence and improve the tidiness of clustering. The lower the value the
        more independent the individual populations, at the cost of more conserved batch effect.
        If ``None``, sets the parameter value automatically to 10 times the total number of
        neighbours for each cell. Set to 0 to skip.
    approx : ``bool``, optional (default: ``True``)
        If ``True``, use annoy's approximate neighbour finding. This results in a quicker run time
        for large datasets while also potentially increasing the degree of batch correction.
    n_trees : ``int``, optional (default: 10)
        Only used when ``approx=True``. The number of trees to construct in the annoy forest.
        More trees give higher precision when querying, at the cost of increased run time and
        resource intensity.
    use_faiss : ``bool``, optional (default: ``True``)
        If ``approx=False`` and the metric is "euclidean", use the faiss package to compute
        nearest neighbours if installed. This improves performance at a minor cost to numerical
        precision as faiss operates on float32.
    metric : ``str`` or ``sklearn.neighbors.DistanceMetric``, optional (default: "angular")
        What distance metric to use. If using ``approx=True``, the options are "angular",
        "euclidean", "manhattan" and "hamming". Otherwise, the options are "euclidean",
        a member of the ``sklearn.neighbors.KDTree.valid_metrics`` list, or parameterised
        ``sklearn.neighbors.DistanceMetric`` `objects
        <https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.DistanceMetric.html>`_:

        >>> from sklearn import neighbors
        >>> neighbors.KDTree.valid_metrics
        ['p', 'chebyshev', 'cityblock', 'minkowski', 'infinity', 'l2', 'euclidean', 'manhattan', 'l1']
        >>> pass_this_as_metric = neighbors.DistanceMetric.get_metric('minkowski',p=3)
    set_op_mix_ratio : ``float``, optional (default: 1)
        UMAP connectivity computation parameter, float between 0 and 1, controlling the
        blend between a connectivity matrix formed exclusively from mutual nearest neighbour
        pairs (0) and a union of all observed neighbour relationships with the mutual pairs
        emphasised (1)
    local_connectivity : ``int``, optional (default: 1)
        UMAP connectivity computation parameter, how many nearest neighbors of each cell
        are assumed to be fully connected (and given a connectivity value of 1)
    copy : ``bool``, optional (default: ``False``)
        If ``True``, return a copy instead of writing to the supplied adata.

    Returns
    -------
    The `adata` with the batch-corrected graph.
    """
    params = locals()  # Has to be first
    kwargs = params.pop('kwargs')
    try:
        from bbknn import bbknn
    except ImportError:
        raise ImportError('Please install bbknn: `pip install bbknn`.')
    return bbknn(**params, **kwargs)
