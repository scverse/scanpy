from __future__ import annotations

from typing import TYPE_CHECKING

from ..._compat import old_positionals
from ..._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata import AnnData
    from sklearn.metrics import DistanceMetric


@old_positionals("batch_key", "use_rep", "approx", "use_annoy", "metric", "copy")
@doctest_needs("bbknn")
def bbknn(  # noqa: PLR0913
    adata: AnnData,
    *,
    batch_key: str = "batch",
    use_rep: str = "X_pca",
    approx: bool = True,
    use_annoy: bool = True,
    metric: str | Callable | DistanceMetric = "euclidean",
    copy: bool = False,
    neighbors_within_batch: int = 3,
    n_pcs: int = 50,
    trim: int | None = None,
    annoy_n_trees: int = 10,
    pynndescent_n_neighbors: int = 30,
    pynndescent_random_state: int = 0,
    use_faiss: bool = True,
    set_op_mix_ratio: float = 1.0,
    local_connectivity: int = 1,
    **kwargs,
) -> AnnData | None:
    """Batch balanced kNN :cite:p:`Polanski2019`.

    Batch balanced kNN alters the kNN procedure to identify each cell's top neighbours in
    each batch separately instead of the entire cell pool with no accounting for batch.
    The nearest neighbours for each batch are then merged to create a final list of
    neighbours for the cell. Aligns batches in a quick and lightweight manner.

    For use in the scanpy workflow as an alternative to :func:`~scanpy.pp.neighbors`.

    .. note::

        This is just a wrapper of :func:`bbknn.bbknn`: up to date docstring,
        more information and bug reports there.

    Params
    ------
    adata
        Needs the PCA computed and stored in `adata.obsm["X_pca"]`.
    batch_key
        `adata.obs` column name discriminating between your batches.
    use_rep
        The dimensionality reduction in `.obsm` to use for neighbour detection. Defaults to PCA.
    approx
        If `True`, use approximate neighbour finding - annoy or PyNNDescent. This results
        in a quicker run time for large datasets while also potentially increasing the degree of
        batch correction.
    use_annoy
        Only used when `approx=True`. If `True`, will use annoy for neighbour finding. If
        `False`, will use pyNNDescent instead.
    metric
        What distance metric to use. The options depend on the choice of neighbour algorithm.

        "euclidean", the default, is always available.

        Annoy supports "angular", "manhattan" and "hamming".

        PyNNDescent supports metrics listed in `pynndescent.distances.named_distances`
        and custom functions, including compiled Numba code.

        >>> import pynndescent
        >>> pynndescent.distances.named_distances.keys()  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        dict_keys(['euclidean', 'l2', 'sqeuclidean', 'manhattan', 'taxicab', 'l1', 'chebyshev', 'linfinity',
        'linfty', 'linf', 'minkowski', 'seuclidean', 'standardised_euclidean', 'wminkowski', ...])

        KDTree supports members of :class:`sklearn.neighbors.KDTree`’s ``valid_metrics`` list, or parameterised
        :class:`~sklearn.metrics.DistanceMetric` objects:

        >>> import sklearn.neighbors
        >>> sklearn.neighbors.KDTree.valid_metrics
        ['euclidean', 'l2', 'minkowski', 'p', 'manhattan', 'cityblock', 'l1', 'chebyshev', 'infinity']

        .. note:: check the relevant documentation for up-to-date lists.
    copy
        If `True`, return a copy instead of writing to the supplied adata.
    neighbors_within_batch
        How many top neighbours to report for each batch; total number of neighbours in
        the initial k-nearest-neighbours computation will be this number times the number
        of batches. This then serves as the basis for the construction of a symmetrical
        matrix of connectivities.
    n_pcs
        How many dimensions (in case of PCA, principal components) to use in the analysis.
    trim
        Trim the neighbours of each cell to these many top connectivities. May help with
        population independence and improve the tidiness of clustering. The lower the value the
        more independent the individual populations, at the cost of more conserved batch effect.
        If `None`, sets the parameter value automatically to 10 times `neighbors_within_batch`
        times the number of batches. Set to 0 to skip.
    annoy_n_trees
        Only used with annoy neighbour identification. The number of trees to construct in the
        annoy forest. More trees give higher precision when querying, at the cost of increased
        run time and resource intensity.
    pynndescent_n_neighbors
        Only used with pyNNDescent neighbour identification. The number of neighbours to include
        in the approximate neighbour graph. More neighbours give higher precision when querying,
        at the cost of increased run time and resource intensity.
    pynndescent_random_state
        Only used with pyNNDescent neighbour identification. The RNG seed to use when creating
        the graph.
    use_faiss
        If `approx=False` and the metric is "euclidean", use the faiss package to compute
        nearest neighbours if installed. This improves performance at a minor cost to numerical
        precision as faiss operates on float32.
    set_op_mix_ratio
        UMAP connectivity computation parameter, float between 0 and 1, controlling the
        blend between a connectivity matrix formed exclusively from mutual nearest neighbour
        pairs (0) and a union of all observed neighbour relationships with the mutual pairs
        emphasised (1)
    local_connectivity
        UMAP connectivity computation parameter, how many nearest neighbors of each cell
        are assumed to be fully connected (and given a connectivity value of 1)

    Returns
    -------
    The `adata` with the batch-corrected graph.

    """
    try:
        from bbknn import bbknn
    except ImportError as e:
        msg = "Please install bbknn: `pip install bbknn`."
        raise ImportError(msg) from e
    return bbknn(
        adata=adata,
        batch_key=batch_key,
        use_rep=use_rep,
        approx=approx,
        use_annoy=use_annoy,
        metric=metric,
        copy=copy,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=trim,
        annoy_n_trees=annoy_n_trees,
        pynndescent_n_neighbors=pynndescent_n_neighbors,
        pynndescent_random_state=pynndescent_random_state,
        use_faiss=use_faiss,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
        **kwargs,
    )
