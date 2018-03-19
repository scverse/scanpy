from ..utils import doc_params
from ._utils import choose_representation, doc_use_rep
from .. import settings
from .. import logging as logg


@doc_params(use_rep=doc_use_rep)
def umap(
        adata,
        use_rep=None,
        n_neighbors=15,
        n_components=2,
        min_dist=0.5,
        spread=1.0,
        metric='euclidean',
        n_epochs=None,
        alpha=1.0,
        gamma=1.0,
        negative_sample_rate=5,
        init='spectral',
        local_connectivity=1.0,
        random_state=0,
        a=None,
        b=None,
        copy=False,
        umap_kwargs={}):
    """\
    UMAP [McInnes18]_.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique for dimension reduction which is suitable for visualization of
    high-dimensional data. It is competitive with tSNE yet, it is faster and it
    arguably preserves more of the global structure.  We use the implementation
    of `umap-learn <https://github.com/lmcinnes/umap>`_ [McInnes18]_.


    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    {use_rep}
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
    min_dist : `float`, optional (default: 0.5)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points
        on the manifold are drawn closer together, while larger values will
        result on a more even dispersal of points. The value should be set
        relative to the ``spread`` value, which determines the scale at which
        embedded points will be spread out.
    spread: `float` (optional, default 1.0)
        The effective scale of embedded points. In combination with ``min_dist``
        this determines how clustered/clumped the embedded points are.
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
    negative_sample_rate: int (optional, default 5)
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
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
    gamma: float (optional, default 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    a: `float` (optional, default `None`)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and
        ``spread``.
    b: `float` (optional, default `None`)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and
        ``spread``.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    X_umap : `adata.obsm`
        UMAP coordinates of data.
    """
    from ..neighbors import umap
    adata = adata.copy() if copy else adata
    logg.info('computing UMAP', r=True)
    if 'neighbors' in adata.uns:
        if not adata.uns['neighbors']['params']['umap']:
            logg.warn('neighbors have not been computed using umap')
        from sklearn.utils import check_random_state
        random_state = check_random_state(random_state)
        from ..neighbors.umap.umap_ import find_ab_params, simplicial_set_embedding
        if n_epochs is None: n_epochs = 0
        if a is None or b is None:
            a, b = find_ab_params(spread, min_dist)
        else:
            a = a
            b = b
        X_umap = simplicial_set_embedding(
            adata.uns['neighbors']['connectivities'].tocoo(),
            n_components,
            alpha,
            a,
            b,
            gamma,
            negative_sample_rate,
            n_epochs,
            init,
            random_state,
            max(0, settings.verbosity-3))
    else:
        X = choose_representation(adata, use_rep=use_rep)
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
        X_umap = um.fit_transform(X)
    adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n'
              '    \'X_umap\', UMAP coordinates (adata.obsm)')
    return adata if copy else None
