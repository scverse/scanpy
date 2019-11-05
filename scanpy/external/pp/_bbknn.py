from typing import Union, Optional

from anndata import AnnData
import sklearn

from ..._utils import lazy_import


# Import this lazily so we don’t slowly import sklearn.stats just for annotation
lazy_import("sklearn.neighbors")
del lazy_import


def bbknn(
    adata: AnnData,
    batch_key: str = 'batch',
    approx: bool = True,
    metric: Union[str, 'sklearn.neighbors.DistanceMetric'] = 'angular',
    copy: bool = False,
    *,
    n_pcs: int = 50,
    trim: Optional[int] = None,
    n_trees: int = 10,
    use_faiss: bool = True,
    set_op_mix_ratio: float = 1.0,
    local_connectivity: int = 1,
    **kwargs,
) -> AnnData:
    """\
    Batch balanced kNN [Polanski19]_.

    Batch balanced kNN alters the kNN procedure to identify each
    cell's top neighbours in each batch separately instead of the
    entire cell pool with no accounting for batch. Aligns batches in a
    quick and lightweight manner.

    For use in the scanpy workflow as an alternative to
    :func:`~scanpy.pp.neighbors`.

    .. note::

        This is just a wrapper of :func:`bbknn.bbknn`: up to date docstring,
        more information and bug reports there.

    Params
    ------
    adata
        Needs the PCA computed and stored in `adata.obsm["X_pca"]`.
    batch_key
        `adata.obs` column name discriminating between your batches.
    approx
        If `True`, use annoy's approximate neighbour finding.
        This results in a quicker run time for large datasets while also
        potentially increasing the degree of batch correction.
    metric
        What distance metric to use. If using `approx=True`, the options are
        `'angular'`, `'euclidean'`, `'manhattan'`, and `'hamming'`.
        Otherwise, the options are `"euclidean"`,
        an element of :class:`sklearn.neighbors.KDTree`’s `valid_metrics`,
        or parameterised :class:`sklearn.neighbors.DistanceMetric` objects:

        >>> from sklearn import neighbors
        >>> neighbors.KDTree.valid_metrics
        ['p', 'chebyshev', 'cityblock', 'minkowski', 'infinity', 'l2', 'euclidean', 'manhattan', 'l1']
        >>> pass_this_as_metric = neighbors.DistanceMetric.get_metric('minkowski',p=3)
    copy
        If `True`, return a copy instead of writing to the supplied adata.
    neighbors_within_batch
        How many top neighbours to report for each batch; total number of neighbours
        will be this number times the number of batches.
    n_pcs
        How many principal components to use in the analysis.
    trim
        Trim the neighbours of each cell to these many top connectivities.
        May help with population independence and improve the tidiness of clustering.
        The lower the value the more independent the individual populations,
        at the cost of more conserved batch effect. If `None`,
        sets the parameter value automatically to 10 times the total number of
        neighbours for each cell. Set to 0 to skip.
    n_trees
        Only used when `approx=True`.
        The number of trees to construct in the annoy forest.
        More trees give higher precision when querying,
        at the cost of increased run time and resource intensity.
    use_faiss
        If `approx=False` and the metric is `"euclidean"`,
        use the `faiss` package to compute nearest neighbours if installed.
        This improves performance at a minor cost to numerical
        precision as `faiss` operates on 32 bit floats.
    set_op_mix_ratio
        UMAP connectivity computation parameter, float between 0 and 1,
        controlling the blend between a connectivity matrix formed exclusively
        from mutual nearest neighbour pairs (0) and a union of all observed
        neighbour relationships with the mutual pairs emphasised (1)
    local_connectivity
        UMAP connectivity computation parameter,
        how many nearest neighbors per cell are assumed to be fully connected
        (and given a connectivity value of 1)

    Returns
    -------
    The `adata` with the batch-corrected graph.
    """
    try:
        from bbknn import bbknn
    except ImportError:
        raise ImportError('Please install bbknn: `pip install bbknn`.')
    return bbknn(
        adata=adata,
        batch_key=batch_key,
        approx=approx,
        metric=metric,
        copy=copy,
        n_pcs=n_pcs,
        trim=trim,
        n_trees=n_trees,
        use_faiss=use_faiss,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
        **kwargs,
    )
