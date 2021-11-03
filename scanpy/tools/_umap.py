from typing import Optional, Union
import warnings

import numpy as np
from packaging import version
from anndata import AnnData
from sklearn.utils import check_random_state, check_array

from ._utils import get_init_pos_from_paga, _choose_representation
from .. import logging as logg
from .._settings import settings
from .._compat import Literal
from .._utils import AnyRandom, NeighborsView


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
    random_state: AnyRandom = 0,
    a: Optional[float] = None,
    b: Optional[float] = None,
    copy: bool = False,
    method: Literal['umap', 'rapids'] = 'umap',
    neighbors_key: Optional[str] = None,
) -> Optional[AnnData]:
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
    adata
        Annotated data matrix.
    min_dist
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
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
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState` or `Generator`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
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
    neighbors_key
        If not specified, umap looks .uns['neighbors'] for neighbors settings
        and .obsp['connectivities'] for connectivities
        (default storage places for pp.neighbors).
        If specified, umap looks .uns[neighbors_key] for neighbors settings and
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_umap** : `adata.obsm` field
        UMAP coordinates of data.
    """
    adata = adata.copy() if copy else adata

    if neighbors_key is None:
        neighbors_key = 'neighbors'

    if neighbors_key not in adata.uns:
        raise ValueError(
            f'Did not find .uns["{neighbors_key}"]. Run `sc.pp.neighbors` first.'
        )
    start = logg.info('computing UMAP')

    neighbors = NeighborsView(adata, neighbors_key)

    if 'params' not in neighbors or neighbors['params']['method'] != 'umap':
        logg.warning(
            f'.obsp["{neighbors["connectivities_key"]}"] have not been computed using umap'
        )

    # Compat for umap 0.4 -> 0.5
    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        import umap

    if version.parse(umap.__version__) >= version.parse("0.5.0"):

        def simplicial_set_embedding(*args, **kwargs):
            from umap.umap_ import simplicial_set_embedding

            X_umap, _ = simplicial_set_embedding(
                *args,
                densmap=False,
                densmap_kwds={},
                output_dens=False,
                **kwargs,
            )
            return X_umap

    else:
        from umap.umap_ import simplicial_set_embedding
    from umap.umap_ import find_ab_params

    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    else:
        a = a
        b = b
    adata.uns['umap'] = {'params': {'a': a, 'b': b}}
    if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif isinstance(init_pos, str) and init_pos == 'paga':
        init_coords = get_init_pos_from_paga(
            adata, random_state=random_state, neighbors_key=neighbors_key
        )
    else:
        init_coords = init_pos  # Let umap handle it
    if hasattr(init_coords, "dtype"):
        init_coords = check_array(init_coords, dtype=np.float32, accept_sparse=False)

    if random_state != 0:
        adata.uns['umap']['params']['random_state'] = random_state
    random_state = check_random_state(random_state)

    neigh_params = neighbors['params']
    X = _choose_representation(
        adata,
        neigh_params.get('use_rep', None),
        neigh_params.get('n_pcs', None),
        silent=True,
    )
    if method == 'umap':
        # the data matrix X is really only used for determining the number of connected components
        # for the init condition in the UMAP embedding
        default_epochs = 500 if neighbors['connectivities'].shape[0] <= 10000 else 200
        n_epochs = default_epochs if maxiter is None else maxiter
        X_umap = simplicial_set_embedding(
            X,
            neighbors['connectivities'].tocoo(),
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

        n_neighbors = neighbors['params']['n_neighbors']
        n_epochs = (
            500 if maxiter is None else maxiter
        )  # 0 is not a valid value for rapids, unlike original umap
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
            random_state=random_state,
        )
        X_umap = umap.fit_transform(X_contiguous)
    adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    logg.info(
        '    finished',
        time=start,
        deep=('added\n' "    'X_umap', UMAP coordinates (adata.obsm)"),
    )
    return adata if copy else None
