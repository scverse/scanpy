from .. import settings
from .. import logging as logg


def umap(
        adata,
        n_neighbors=15,
        n_components=2,
        min_dist=0.1,
        metric='euclidean',
        alpha=1.0,
        init='spectral',
        local_connectivity=1.0,
        random_state=None,
        copy=False,
        umap_kwargs={}):
    """UMAP [McInnes18]_.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique for dimension reduction which is suitable for visualization of
    high-dimensional data. It is competitive with tSNE yet, it is faster and it
    arguably preserves more of the global structure.  We use the implementation
    of `umap-learn <https://github.com/lmcinnes/umap>`_ [McInnes18]_.

    *Note:* In contrast to other embeddings within Scanpy, UMAP does *not* use a
    PCA-reduced version of the data. We might change this behavior in the
    future.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    n_neighbors : `float`, optional (default: 15)
        The size of local neighborhood (in terms of number of neighboring
        sample points) used for manifold approximation. Larger values
        result in more global views of the manifold, while smaller
        values result in more local data being preserved. In general
        values should be in the range 2 to 100.
    n_components : `int`, optional (default: 2)
        The dimension of the space to embed into. This defaults to 2 to
        provide easy visualization, but can reasonably be set to any
        integer value in the range 2 to 100.
    min_dist : `float`, optional (default: 0.1)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points
        on the manifold are drawn closer together, while larger values will
        result on a more even dispersal of points. The value should be set
        relative to the ``spread`` value, which determines the scale at which
        embedded points will be spread out.
    metric : `string` or `function`, optional (default: 'euclidean')
        The metric to use to compute distances in high dimensional space.
        If a string is passed it must match a valid predefined metric. If
        a general metric is required a function that takes two 1d arrays and
        returns a float can be provided. For performance purposes it is
        required that this be a numba jit'd function. Valid string metrics
        include:
            * euclidean
            * manhattan
            * chebyshev
            * minkowski
            * canberra
            * braycurtis
            * mahalanobis
            * wminkowski
            * seuclidean
            * cosine
            * correlation
            * haversine
            * hamming
            * jaccard
            * dice
            * russelrao
            * kulsinski
            * rogerstanimoto
            * sokalmichener
            * sokalsneath
            * yule
        Metrics that take arguments (such as minkowski, mahalanobis etc.)
        can have arguments passed via the metric_kwds dictionary. At this
        time care must be taken and dictionary elements must be ordered
        appropriately; this will hopefully be fixed in the future.
    alpha : `float`, optional (default: 1.0)
        The initial learning rate for the embedding optimization.
    init : `string` or `np.array`, optional (default: 'spectral')
        How to initialize the low dimensional embedding.
        Options are:
            * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
            * 'random': assign initial embedding positions at random.
            * A numpy array of initial embedding positions.
    local_connectivity : `int`, optional (default: 1)
        The local connectivity required -- i.e. the number of nearest
        neighbors that should be assumed to be connected at a local level.
        The higher this value the more connected the manifold becomes
        locally. In practice this should be not more than the local intrinsic
        dimension of the manifold.
    random_state : `int`, `RandomState instance` or `None`, optional (default: `None`)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.
    umap_kwargs : `dict`, optional (default: `{}`)
        Additional keyword arguments for UMAP class constructor from umap-learn
        package.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_umap : `np.ndarray` (`adata.obs`, dtype `float`)
        UMAP coordinates of data.
    """
    try:
        import umap
    except ImportError:
        logg.error('UMAP visualization requires umap-learn package, however '
                   'it is not found. Follow instructions on GitHub to install '
                   'umap-learn: https://github.com/lmcinnes/umap#installing')

    logg.info('computing UMAP', r=True)
    adata = adata.copy() if copy else adata

    # params for umap-learn
    params_umap = {'n_neighbors': n_neighbors,
                   'n_components': n_components,
                   'min_dist': min_dist,
                   'metric': metric,
                   'alpha': alpha,
                   'init': init,
                   'local_connectivity': local_connectivity,
                   'random_state': random_state,
                   'verbose': max(0, settings.verbosity-3),
                   **umap_kwargs
                   }

    um = umap.UMAP(**params_umap)
    X_umap = um.fit_transform(adata.X)

    # update AnnData instance
    adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'X_umap\', UMAP coordinates (adata.obs)')

    return adata if copy else None
