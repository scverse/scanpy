from ._utils import get_init_pos_from_paga
from .._settings import settings
from .. import logging as logg
from ..logging import (
    _settings_verbosity_greater_or_equal_than,
    _VERBOSITY_LEVELS_FROM_STRINGS,
)

def umap(
    adata,
    min_dist=0.5,
    spread=1.0,
    n_components=2,
    maxiter=None,
    alpha=1.0,
    gamma=1.0,
    negative_sample_rate=5,
    init_pos='spectral',
    random_state=0,
    a=None,
    b=None,
    copy=False,
):
    """Embed the neighborhood graph using UMAP [McInnes18]_.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.  We use the
    implementation of `umap-learn <https://github.com/lmcinnes/umap>`__
    [McInnes18]_. For a few comparisons of UMAP with tSNE, see this `preprint
    <https://doi.org/10.1101/298430>`__.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    min_dist : `float`, optional (default: 0.5)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread : `float` (optional, default 1.0)
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components : `int`, optional (default: 2)
        The number of dimensions of the embedding.
    maxiter : `int`, optional (default: `None`)
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha : `float`, optional (default: 1.0)
        The initial learning rate for the embedding optimization.
    gamma : `float` (optional, default 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate : `int` (optional, default 5)
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos : `string` or `np.array`, optional (default: 'spectral')
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP.
        Options are:

        * Any key for `adata.obsm`.
        * 'paga': positions from :func:`~scanpy.api.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state : `int`, `RandomState` or `None`, optional (default: 0)
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a : `float` (optional, default `None`)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b : `float` (optional, default `None`)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_umap** : `adata.obsm` field
        UMAP coordinates of data.
    """
    adata = adata.copy() if copy else adata
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'Did not find \'neighbors/connectivities\'. Run `sc.pp.neighbors` first.')
    logg.info('computing UMAP', r=True)
    if ('params' not in adata.uns['neighbors']
        or adata.uns['neighbors']['params']['method'] != 'umap'):
        logg.warn('neighbors/connectivities have not been computed using umap')
    from ..neighbors.umap.umap_ import find_ab_params, simplicial_set_embedding
    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    else:
        a = a
        b = b
    if init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif init_pos == 'paga':
        init_coords = get_init_pos_from_paga(adata, random_state=random_state)
    else:
        init_coords = init_pos
    from sklearn.utils import check_random_state
    random_state = check_random_state(random_state)
    n_epochs = maxiter
    verbosity = _VERBOSITY_LEVELS_FROM_STRINGS.get(settings.verbosity, settings.verbosity)
    X_umap = simplicial_set_embedding(
        adata.uns['neighbors']['connectivities'].tocoo(),
        n_components,
        alpha,
        a,
        b,
        gamma,
        negative_sample_rate,
        n_epochs,
        init_coords,
        random_state,
        max(0, verbosity-3))
    adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    logg.info('    finished', time=True, end=' ' if _settings_verbosity_greater_or_equal_than(3) else '\n')
    logg.hint('added\n'
              '    \'X_umap\', UMAP coordinates (adata.obsm)')
    return adata if copy else None
