"""\
Perform clustering using PhenoGraph
"""
from typing import Union, Tuple, Optional, Type, Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import spmatrix

from ..._compat import Literal
from ...tools._leiden import MutableVertexPartition
from ... import logging as logg


def phenograph(
    adata: Union[AnnData, np.ndarray, spmatrix],
    clustering_algo: Optional[Literal['louvain', 'leiden']] = 'louvain',
    k: int = 30,
    directed: bool = False,
    prune: bool = False,
    min_cluster_size: int = 10,
    jaccard: bool = True,
    primary_metric: Literal[
        'euclidean',
        'manhattan',
        'correlation',
        'cosine',
    ] = "euclidean",
    n_jobs: int = -1,
    q_tol: float = 1e-3,
    louvain_time_limit: int = 2000,
    nn_method: Literal['kdtree', 'brute'] = 'kdtree',
    partition_type: Optional[Type[MutableVertexPartition]] = None,
    resolution_parameter: float = 1,
    n_iterations: int = -1,
    use_weights: bool = True,
    seed: Optional[int] = None,
    copy: bool = False,
    **kargs: Any,
) -> Tuple[Optional[np.ndarray], spmatrix, Optional[float]]:
    """\
    PhenoGraph clustering [Levine15]_.

    **PhenoGraph** is a clustering method designed for high-dimensional single-cell
    data. It works by creating a graph ("network") representing phenotypic similarities
    between cells and then identifying communities in this graph. It supports both
    Louvain_ and Leiden_ algorithms for community detection.

    .. _Louvain: https://louvain-igraph.readthedocs.io/en/latest/

    .. _Leiden: https://leidenalg.readthedocs.io/en/latest/reference.html

    .. note::
       More information and bug reports `here
       <https://github.com/dpeerlab/PhenoGraph>`__.

    Parameters
    ----------
    adata
        AnnData, or Array of data to cluster, or sparse matrix of k-nearest neighbor
        graph. If ndarray, n-by-d array of n cells in d dimensions. if sparse matrix,
        n-by-n adjacency matrix.
    clustering_algo
        Choose between `'Louvain'` or `'Leiden'` algorithm for clustering.
    k
        Number of nearest neighbors to use in first step of graph construction.
    directed
        Whether to use a symmetric (default) or asymmetric (`'directed'`) graph.
        The graph construction process produces a directed graph, which is symmetrized
        by one of two methods (see `prune` below).
    prune
        `prune=False`, symmetrize by taking the average between the graph and its
        transpose. `prune=True`, symmetrize by taking the product between the graph
        and its transpose.
    min_cluster_size
        Cells that end up in a cluster smaller than min_cluster_size are considered
        outliers and are assigned to -1 in the cluster labels.
    jaccard
        If `True`, use Jaccard metric between k-neighborhoods to build graph. If
        `False`, use a Gaussian kernel.
    primary_metric
        Distance metric to define nearest neighbors. Note that performance will be
        slower for correlation and cosine.
    n_jobs
        Nearest Neighbors and Jaccard coefficients will be computed in parallel using
        n_jobs. If 1 is given, no parallelism is used. If set to -1, all CPUs are used.
        For n_jobs below -1, `n_cpus + 1 + n_jobs` are used.
    q_tol
        Tolerance, i.e. precision, for monitoring modularity optimization.
    louvain_time_limit
        Maximum number of seconds to run modularity optimization. If exceeded the best
        result so far is returned.
    nn_method
        Whether to use brute force or kdtree for nearest neighbor search.
        For very large high-dimensional data sets, brute force, with parallel
        computation, performs faster than kdtree.
    partition_type
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`. For the
        available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    resolution_parameter
        A parameter value controlling the coarseness of the clustering in Leiden. Higher
        values lead to more clusters. Set to `None` if overriding `partition_type` to
        one that does not accept a `resolution_parameter`.
    n_iterations
        Number of iterations to run the Leiden algorithm. If the number of iterations is
        negative, the Leiden algorithm is run until an iteration in which there was no
        improvement.
    use_weights
        Use vertices in the Leiden computation.
    seed
        Leiden initialization of the optimization.
    copy
        Return a copy or write to `adata`.
    kargs
        Additional arguments passed to :func:`~leidenalg.find_partition` and the
        constructor of the `partition_type`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields:

    **communities** - :class:`~numpy.ndarray` (:attr:`~anndata.AnnData.obs`, dtype `int`)
        integer array of community assignments for each row in data.

    **graph** - :class:`~scipy.sparse.spmatrix` (:attr:`~anndata.AnnData.obsp`, dtype `float`)
        the graph that was used for clustering.

    **Q** - `float` (:attr:`~anndata.AnnData.uns`, dtype `float`)
        the modularity score for communities on graph.

    Example
    -------
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> import numpy as np
    >>> import pandas as pd

    With annotated data as input:

    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.normalize_per_cell(adata)

    Then do PCA:

    >>> sc.tl.pca(adata, n_comps=100)

    Compute phenograph clusters:

    **Louvain** community detection

    >>> sce.tl.phenograph(adata, clustering_algo="louvain", k=30)

    **Leiden** community detection

    >>> sce.tl.phenograph(adata, clustering_algo="leiden", k=30)

    Return only `Graph` object

    >>> sce.tl.phenograph(adata, clustering_algo=None, k=30)

    Now to show phenograph on tSNE (for example):

    Compute tSNE:

    >>> sc.tl.tsne(adata, random_state=7)

    Plot phenograph clusters on tSNE:

    >>> sc.pl.tsne(
    ...     adata, color = ["pheno_louvain", "pheno_leiden"], s = 100,
    ...     palette = sc.pl.palettes.vega_20_scanpy, legend_fontsize = 10
    ... )

    Cluster and cluster centroids for input Numpy ndarray

    >>> df = np.random.rand(1000, 40)
    >>> dframe = pd.DataFrame(df)
    >>> dframe.index, dframe.columns = (map(str, dframe.index), map(str, dframe.columns))
    >>> adata = AnnData(dframe)
    >>> sc.tl.pca(adata, n_comps=20)
    >>> sce.tl.phenograph(adata, clustering_algo="leiden", k=50)
    >>> sc.tl.tsne(adata, random_state=1)
    >>> sc.pl.tsne(
    ...     adata, color=['pheno_leiden'], s=100,
    ...     palette=sc.pl.palettes.vega_20_scanpy, legend_fontsize=10
    ... )
    """
    start = logg.info("PhenoGraph clustering")

    try:
        import phenograph

        assert phenograph.__version__ >= "1.5.3"
    except (ImportError, AssertionError, AttributeError):
        raise ImportError(
            "please install the latest release of phenograph:\n\t"
            "pip install -U PhenoGraph"
        )

    if isinstance(adata, AnnData):
        try:
            data = adata.obsm["X_pca"]
        except KeyError:
            raise KeyError("Please run `sc.tl.pca` on `adata` and try again!")
    else:
        data = adata
        copy = True

    comm_key = (
        "pheno_{}".format(clustering_algo)
        if clustering_algo in ["louvain", "leiden"]
        else ''
    )
    ig_key = "pheno_{}_ig".format("jaccard" if jaccard else "gaussian")
    q_key = "pheno_{}_q".format("jaccard" if jaccard else "gaussian")

    communities, graph, Q = phenograph.cluster(
        data=data,
        clustering_algo=clustering_algo,
        k=k,
        directed=directed,
        prune=prune,
        min_cluster_size=min_cluster_size,
        jaccard=jaccard,
        primary_metric=primary_metric,
        n_jobs=n_jobs,
        q_tol=q_tol,
        louvain_time_limit=louvain_time_limit,
        nn_method=nn_method,
        partition_type=partition_type,
        resolution_parameter=resolution_parameter,
        n_iterations=n_iterations,
        use_weights=use_weights,
        seed=seed,
        **kargs,
    )

    logg.info("    finished", time=start)

    if copy:
        return communities, graph, Q
    else:
        adata.obsp[ig_key] = graph
        if comm_key:
            adata.obs[comm_key] = pd.Categorical(communities)
        if Q:
            adata.uns[q_key] = Q
