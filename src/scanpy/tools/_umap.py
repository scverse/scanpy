from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

import numpy as np
from sklearn.utils import check_array, check_random_state

from .. import logging as logg
from .._compat import old_positionals
from .._settings import settings
from .._utils import NeighborsView
from ._utils import _choose_representation, get_init_pos_from_paga

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

    from .._compat import _LegacyRandom
    from ._types import DensmapMethodKwds
    from .._utils.random import _LegacyRandom

    _InitPos = Literal["paga", "spectral", "random"]


@old_positionals(
    "min_dist",
    "spread",
    "n_components",
    "maxiter",
    "alpha",
    "gamma",
    "negative_sample_rate",
    "init_pos",
    "random_state",
    "a",
    "b",
    "copy",
    "method",
    "neighbors_key",
)
def umap(  # noqa: PLR0913, PLR0915
    adata: AnnData,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: int | None = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: _InitPos | np.ndarray | None = "spectral",
    random_state: _LegacyRandom = 0,
    a: float | None = None,
    b: float | None = None,
    method: Literal["umap", "rapids", "densmap"] = "umap",
    method_kwds: DensmapMethodKwds | None = None,
    key_added: str | None = None,
    neighbors_key: str = "neighbors",
    copy: bool = False,
) -> AnnData | None:
    r"""Embed the neighborhood graph using UMAP :cite:p:`McInnes2018`.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.
    We use the implementation of umap-learn_ :cite:p:`McInnes2018`.
    For a few comparisons of UMAP with tSNE, see :cite:t:`Becht2018`.

    .. _umap-learn: https://github.com/lmcinnes/umap

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
    method
        Chosen implementation.

        ``'umap'``
            Umap’s simplical set embedding.
        ``'densmap'``
            Umap’s simplical set embedding with densmap=True :cite:p:`Narayan2021`.
        ``'rapids'``
            GPU accelerated implementation.

            .. deprecated:: 1.10.0
                Use :func:`rapids_singlecell.tl.umap` instead.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_umap']` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['umap']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]`` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    neighbors_key
        Umap looks in
        :attr:`~anndata.AnnData.uns`\ ``[neighbors_key]`` for neighbors settings and
        :attr:`~anndata.AnnData.obsp`\ ``[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
    copy
        Return a copy instead of writing to adata.

    method_kwds
        Additional method parameters.

        If method is ``'densmap'``, the following parameters are available:

            ``dens_lambda`` : `float`, optional (default: 2.0)
                Controls the regularization weight of the density correlation term
                in densMAP. Higher values prioritize density preservation over the
                UMAP objective, and vice versa for values closer to zero. Setting this
                parameter to zero is equivalent to running the original UMAP algorithm.
            ``dens_frac`` : `float`, optional (default: 0.3)
                Controls the fraction of epochs (between 0 and 1) where the
                density-augmented objective is used in densMAP. The first
                (1 - dens_frac) fraction of epochs optimize the original UMAP objective
                before introducing the density correlation term.
            ``dens_var_shift`` : `float`, optional (default: 0.1)
                A small constant added to the variance of local radii in the
                embedding when calculating the density correlation objective to
                prevent numerical instability from dividing by a small number

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields unless method is 'densmap':

    `adata.obsm['X_umap' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        UMAP coordinates of data.
    `adata.uns['umap' | key_added]` : :class:`dict`
        UMAP parameters.

    When method is 'densmap', sets the following fields:

    `adata.obsm['X_densmap']` : :class:`numpy.ndarray` (dtype `float`)
        densMAP coordinates of data.
    `adata.uns['densmap']` : :class:`dict`
        densMAP parameters.

    """
    adata = adata.copy() if copy else adata

    key_obsm, key_uns = (
        (("X_densmap", "densmap") if method == "densmap" else ("X_umap", "umap"))
        if key_added is None
        else [key_added] * 2
    )
    method_name = "DensMAP" if method == "densmap" else "UMAP"

    if neighbors_key is None:  # backwards compat
        neighbors_key = "neighbors"
    if neighbors_key not in adata.uns:
        msg = f"Did not find .uns[{neighbors_key!r}]. Run `sc.pp.neighbors` first."
        raise ValueError(msg)

    start = logg.info("computing UMAP")

    neighbors = NeighborsView(adata, neighbors_key)

    if "params" not in neighbors or neighbors["params"]["method"] != "umap":
        logg.warning(
            f'.obsp["{neighbors["connectivities_key"]}"] have not been computed using umap'
        )

    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        import umap

    from umap.umap_ import find_ab_params, simplicial_set_embedding

    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)

    adata.uns[key_uns] = dict(params=dict(a=a, b=b))
    if isinstance(init_pos, str) and init_pos in adata.obsm:
        init_coords = adata.obsm[init_pos]
    elif isinstance(init_pos, str) and init_pos == "paga":
        init_coords = get_init_pos_from_paga(
            adata, random_state=random_state, neighbors_key=neighbors_key
        )
    else:
        init_coords = init_pos  # Let umap handle it
    if hasattr(init_coords, "dtype"):
        init_coords = check_array(init_coords, dtype=np.float32, accept_sparse=False)

    if random_state != 0:
        adata.uns[key_uns]["params"]["random_state"] = random_state
    random_state = check_random_state(random_state)

    neigh_params = neighbors["params"]
    x = _choose_representation(
        adata,
        use_rep=neigh_params.get("use_rep", None),
        n_pcs=neigh_params.get("n_pcs", None),
        silent=True,
    )

    if method_kwds is None:
        method_kwds = {}

    densmap_kwds = (
        {
            "graph_dists": neighbors["distances"],
            "n_neighbors": neigh_params.get("n_neighbors", 15),
            # Default params from umap package
            # Reference: https://github.com/lmcinnes/umap/blob/868e55cb614f361a0d31540c1f4a4b175136025c/umap/umap_.py#L1692
            # If user provided method_kwds, the user-provided values should
            # overwrite the default values specified above.
            "lambda": method_kwds.get("dens_lambda", 2.0),
            "frac": method_kwds.get("dens_frac", 0.3),
            "var_shift": method_kwds.get("dens_var_shift", 0.1),
        }
        if method == "densmap"
        else {}
    )
    if method == "densmap":
        adata.uns[key_uns]["params"].update(
            {
                "dens_lambda": densmap_kwds["lambda"],
                "dens_frac": densmap_kwds["frac"],
                "dens_var_shift": densmap_kwds["var_shift"],
            }
        )

    if method == "umap" or method == "densmap":
        # the data matrix X is really only used for determining the number of connected components
        # for the init condition in the UMAP embedding
        default_epochs = 500 if neighbors["connectivities"].shape[0] <= 10000 else 200
        n_epochs = default_epochs if maxiter is None else maxiter
        x_umap, _ = simplicial_set_embedding(
            data=x,
            graph=neighbors["connectivities"].tocoo(),
            n_components=n_components,
            initial_alpha=alpha,
            a=a,
            b=b,
            gamma=gamma,
            negative_sample_rate=negative_sample_rate,
            n_epochs=n_epochs,
            init=init_coords,
            random_state=random_state,
            metric=neigh_params.get("metric", "euclidean"),
            metric_kwds=neigh_params.get("metric_kwds", {}),
            densmap=(method == "densmap"),
            densmap_kwds=densmap_kwds,
            output_dens=False,
            verbose=settings.verbosity > 3,
        )
    elif method == "rapids":
        msg = (
            "`method='rapids'` is deprecated. "
            "Use `rapids_singlecell.tl.louvain` instead."
        )
        warnings.warn(msg, FutureWarning, stacklevel=2)
        metric = neigh_params.get("metric", "euclidean")
        if metric != "euclidean":
            msg = (
                f"`sc.pp.neighbors` was called with `metric` {metric!r}, "
                "but umap `method` 'rapids' only supports the 'euclidean' metric."
            )
            raise ValueError(msg)
        from cuml import UMAP

        n_neighbors = neighbors["params"]["n_neighbors"]
        n_epochs = (
            500 if maxiter is None else maxiter
        )  # 0 is not a valid value for rapids, unlike original umap
        x_contiguous = np.ascontiguousarray(x, dtype=np.float32)
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
        x_umap = umap.fit_transform(x_contiguous)
    adata.obsm[key_obsm] = x_umap  # annotate samples with UMAP coordinates
    logg.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {key_obsm!r}, {method_name} coordinates (adata.obsm)\n"
            f"    {key_uns!r}, {method_name} parameters (adata.uns)"
        ),
    )
    return adata if copy else None
