"""\
Embed high-dimensional data using TriMap
"""
from typing import Optional, Union

from anndata import AnnData
import scipy.sparse as scp

from ..._settings import settings
from ..._compat import Literal
from ... import logging as logg


def trimap(
    adata: AnnData,
    n_components: int = 2,
    n_inliers: int = 10,
    n_outliers: int = 5,
    n_random: int = 5,
    metric: Literal['angular', 'euclidean', 'hamming', 'manhattan'] = 'euclidean',
    weight_adj: float = 500.0,
    lr: float = 1000.0,
    n_iters: int = 400,
    verbose: Union[bool, int, None] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    TriMap: Large-scale Dimensionality Reduction Using Triplets [Amid19]_.

    TriMap is a dimensionality reduction method that uses triplet constraints
    to form a low-dimensional embedding of a set of points. The triplet
    constraints are of the form "point i is closer to point j than point k".
    The triplets are sampled from the high-dimensional representation of the
    points and a weighting scheme is used to reflect the importance of each
    triplet.

    TriMap provides a significantly better global view of the data than the
    other dimensionality reduction methods such t-SNE, LargeVis, and UMAP.
    The global structure includes relative distances of the clusters, multiple
    scales in the data, and the existence of possible outliers. We define a
    global score to quantify the quality of an embedding in reflecting the
    global structure of the data.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_components
        Number of dimensions of the embedding.
    n_inliers
        Number of inlier points for triplet constraints.
    n_outliers
        Number of outlier points for triplet constraints.
    n_random
        Number of random triplet constraints per point.
    metric
        Distance measure: 'angular', 'euclidean', 'hamming', 'manhattan'.
    weight_adj
        Adjusting the weights using a non-linear transformation.
    lr
        Learning rate.
    n_iters
        Number of iterations.
    verbose
        If `True`, print the progress report.
        If `None`, `sc.settings.verbosity` is used.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_trimap** : :class:`~numpy.ndarray`, (:attr:`~anndata.AnnData.obsm`, shape=(n_samples, n_components), dtype `float`)
        TriMap coordinates of data.

    Example
    -------
    
    >>> import scanpy as sc
    >>> import scanpy.external as sce
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> pbmc = sce.tl.trimap(pbmc, copy=True)
    >>> sce.pl.trimap(pbmc, color=['bulk_labels'], s=10)
    """

    try:
        from trimap import TRIMAP
    except ImportError:
        raise ImportError('\nplease install trimap: \n\n\tsudo pip install trimap')
    adata = adata.copy() if copy else adata
    start = logg.info('computing TriMap')
    adata = adata.copy() if copy else adata
    verbosity = settings.verbosity if verbose is None else verbose
    verbose = verbosity if isinstance(verbosity, bool) else verbosity > 0

    if 'X_pca' in adata.obsm:
        n_dim_pca = adata.obsm['X_pca'].shape[1]
        X = adata.obsm['X_pca'][:, : min(n_dim_pca, 100)]
    else:
        X = adata.X
        if scp.issparse(X):
            raise ValueError(
                'trimap currently does not support sparse matrices. Please'
                'use a dense matrix or apply pca first.'
            )
        logg.warning('`X_pca` not found. Run `sc.pp.pca` first for speedup.')
    X_trimap = TRIMAP(
        n_dims=n_components,
        n_inliers=n_inliers,
        n_outliers=n_outliers,
        n_random=n_random,
        lr=lr,
        distance=metric,
        weight_adj=weight_adj,
        n_iters=n_iters,
        verbose=verbose,
    ).fit_transform(X)
    adata.obsm['X_trimap'] = X_trimap
    logg.info(
        '    finished',
        time=start,
        deep="added\n    'X_trimap', TriMap coordinates (adata.obsm)",
    )
    return adata if copy else None
