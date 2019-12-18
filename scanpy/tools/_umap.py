<<<<<<< HEAD
from ._utils import get_init_pos_from_paga
from .. import settings
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
        copy=False):
=======
from typing import Optional, Union

import numpy as np
from anndata import AnnData
from numpy.random.mtrand import RandomState
from sklearn.utils import check_random_state, check_array

from ._utils import get_init_pos_from_paga, _choose_representation
from .. import logging as logg
from .._settings import settings
from .._compat import Literal


_InitPos = Literal['paga', 'spectral', 'random']


def umap(
    adata: AnnData,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: Optional[int] = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: Union[_InitPos, np.ndarray, None] = 'spectral',
    random_state: Optional[Union[int, RandomState]] = 0,
    a: Optional[float] = None,
    b: Optional[float] = None,
    copy: bool = False,
    method: Literal['umap', 'rapids'] = 'umap'
) -> Optional[AnnData]:
>>>>>>> upstream/master
    """\
    Embed the neighborhood graph using UMAP [McInnes18]_.

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
<<<<<<< HEAD
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    min_dist : `float`, optional (default: 0.5)
=======
    adata
        Annotated data matrix.
    min_dist
>>>>>>> upstream/master
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
<<<<<<< HEAD
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
=======
    spread
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components
        The number of dimensions of the embedding.
    maxiter
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha
        The initial learning rate for the embedding optimization.
    gamma
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Options are:

        * Any key for `adata.obsm`.
        * 'paga': positions from :func:`~scanpy.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state
>>>>>>> upstream/master
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
<<<<<<< HEAD
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
=======
    a
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    copy
        Return a copy instead of writing to adata.
    method
        Use the original 'umap' implementation, or 'rapids' (experimental, GPU only)
>>>>>>> upstream/master

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

<<<<<<< HEAD
    X_umap : `adata.obsm`
=======
    **X_umap** : `adata.obsm` field
>>>>>>> upstream/master
        UMAP coordinates of data.
    """
    adata = adata.copy() if copy else adata
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'Did not find \'neighbors/connectivities\'. Run `sc.pp.neighbors` first.')
<<<<<<< HEAD
    logg.info('computing UMAP', r=True)
    if ('params' not in adata.uns['neighbors']
        or adata.uns['neighbors']['params']['method'] != 'umap'):
        logg.warn('neighbors/connectivities have not been computed using umap')
    from ..neighbors.umap.umap_ import find_ab_params, simplicial_set_embedding
=======
    start = logg.info('computing UMAP')
    if ('params' not in adata.uns['neighbors']
        or adata.uns['neighbors']['params']['method'] != 'umap'):
        logg.warning('neighbors/connectivities have not been computed using umap')
    from umap.umap_ import find_ab_params, simplicial_set_embedding
>>>>>>> upstream/master
    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    else:
        a = a
        b = b
<<<<<<< HEAD
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
=======
    adata.uns['umap'] = {'params':{'a': a, 'b': b}}
    if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif isinstance(init_pos, str) and init_pos == 'paga':
        init_coords = get_init_pos_from_paga(adata, random_state=random_state)
    else:
        init_coords = init_pos  # Let umap handle it
    if hasattr(init_coords, "dtype"):
        init_coords = check_array(init_coords, dtype=np.float32, accept_sparse=False)

    random_state = check_random_state(random_state)
    neigh_params = adata.uns['neighbors']['params']
    X = _choose_representation(
        adata, neigh_params.get('use_rep', None), neigh_params.get('n_pcs', None), silent=True)
    if method == 'umap':
        # the data matrix X is really only used for determining the number of connected components
        # for the init condition in the UMAP embedding
        n_epochs = 0 if maxiter is None else maxiter
        X_umap = simplicial_set_embedding(
            X,
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
            neigh_params.get('metric', 'euclidean'),
            neigh_params.get('metric_kwds', {}),
            verbose=settings.verbosity > 3,
        )
    elif method == 'rapids':
        metric = neigh_params.get('metric', 'euclidean')
        if metric != 'euclidean':
            raise ValueError(
                f'`sc.pp.neighbors` was called with `metric` {metric!r}, '
                "but umap `method` 'rapids' only supports the 'euclidean' metric."
            )
        from cuml import UMAP
        n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
        n_epochs = 500 if maxiter is None else maxiter # 0 is not a valid value for rapids, unlike original umap
        X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
        umap = UMAP(
            n_neighbors=n_neighbors,
            n_components=n_components,
            n_epochs=n_epochs,
            learning_rate=alpha,
            init=init_pos,
            min_dist=min_dist,
            spread=spread,
            negative_sample_rate=negative_sample_rate,
            a=a,
            b=b,
            verbose=settings.verbosity > 3,
        )
        X_umap = umap.fit_transform(X_contiguous)
    adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    logg.info(
        '    finished',
        time=start,
        deep=(
            'added\n'
            "    'X_umap', UMAP coordinates (adata.obsm)"
        ),
    )
>>>>>>> upstream/master
    return adata if copy else None
