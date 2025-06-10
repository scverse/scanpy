from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from packaging.version import Version

from .. import logging as logg
from .._compat import old_positionals
from .._settings import settings
from .._utils import _doc_params, raise_not_implemented_error_if_backed_type
from ..neighbors._doc import doc_n_pcs, doc_use_rep
from ._utils import _choose_representation

if TYPE_CHECKING:
    from anndata import AnnData

    from .._utils.random import _LegacyRandom


@old_positionals(
    "use_rep",
    "perplexity",
    "early_exaggeration",
    "learning_rate",
    "random_state",
    "use_fast_tsne",
    "n_jobs",
    "copy",
)
@_doc_params(doc_n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def tsne(  # noqa: PLR0913
    adata: AnnData,
    n_pcs: int | None = None,
    *,
    use_rep: str | None = None,
    perplexity: float = 30,
    metric: str = "euclidean",
    early_exaggeration: float = 12,
    learning_rate: float = 1000,
    random_state: _LegacyRandom = 0,
    use_fast_tsne: bool = False,
    n_jobs: int | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    r"""t-SNE :cite:p:`vanDerMaaten2008,Amir2013,Pedregosa2011`.

    t-distributed stochastic neighborhood embedding (tSNE, :cite:t:`vanDerMaaten2008`) was
    proposed for visualizating single-cell data by :cite:t:`Amir2013`. Here, by default,
    we use the implementation of *scikit-learn* :cite:p:`Pedregosa2011`. You can achieve
    a huge speedup and better convergence if you install Multicore-tSNE_
    by :cite:t:`Ulyanov2016`, which will be automatically detected by Scanpy.

    .. _multicore-tsne: https://github.com/DmitryUlyanov/Multicore-TSNE

    Parameters
    ----------
    adata
        Annotated data matrix.
    {doc_n_pcs}
    {use_rep}
    perplexity
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    metric
        Distance metric calculate neighbors on.
    early_exaggeration
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state
        Change this to use different intial states for the optimization.
        If `None`, the initial state is not reproducible.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy.settings.n_jobs`.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_tsne']` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['tsne']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]`` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsm['X_tsne' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        tSNE coordinates of data.
    `adata.uns['tsne' | key_added]` : :class:`dict`
        tSNE parameters.

    """
    import sklearn

    start = logg.info("computing tSNE")
    adata = adata.copy() if copy else adata
    X = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
    raise_not_implemented_error_if_backed_type(X, "tsne")
    # params for sklearn
    n_jobs = settings.n_jobs if n_jobs is None else n_jobs
    params_sklearn = dict(
        perplexity=perplexity,
        random_state=random_state,
        verbose=settings.verbosity > 3,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_jobs=n_jobs,
        metric=metric,
    )
    if metric != "euclidean" and (Version(sklearn.__version__) < Version("1.3.0rc1")):
        params_sklearn["square_distances"] = True

    # Backwards compat handling: Remove in scanpy 1.9.0
    if n_jobs != 1 and not use_fast_tsne:
        warnings.warn(
            "In previous versions of scanpy, calling tsne with n_jobs > 1 would use "
            "MulticoreTSNE. Now this uses the scikit-learn version of TSNE by default. "
            "If you'd like the old behaviour (which is deprecated), pass "
            "'use_fast_tsne=True'. Note, MulticoreTSNE is not actually faster anymore.",
            UserWarning,
            stacklevel=2,
        )
    if use_fast_tsne:
        warnings.warn(
            "Argument `use_fast_tsne` is deprecated, and support for MulticoreTSNE "
            "will be dropped in a future version of scanpy.",
            FutureWarning,
            stacklevel=2,
        )

    # deal with different tSNE implementations
    if use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE

            tsne = TSNE(**params_sklearn)
            logg.info("    using the 'MulticoreTSNE' package by Ulyanov (2017)")
            # need to transform to float64 for MulticoreTSNE...
            X_tsne = tsne.fit_transform(X.astype("float64"))
        except ImportError:
            use_fast_tsne = False
            warnings.warn(
                "Could not import 'MulticoreTSNE'. Falling back to scikit-learn.",
                UserWarning,
                stacklevel=2,
            )
    if use_fast_tsne is False:  # In case MultiCore failed to import
        from sklearn.manifold import TSNE

        # unfortunately, sklearn does not allow to set a minimum number
        # of iterations for barnes-hut tSNE
        tsne = TSNE(**params_sklearn)
        logg.info("    using sklearn.manifold.TSNE")
        X_tsne = tsne.fit_transform(X)

    # update AnnData instance
    params = dict(
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        n_jobs=n_jobs,
        metric=metric,
        use_rep=use_rep,
    )
    key_uns, key_obsm = ("tsne", "X_tsne") if key_added is None else [key_added] * 2
    adata.obsm[key_obsm] = X_tsne  # annotate samples with tSNE coordinates
    adata.uns[key_uns] = dict(params={k: v for k, v in params.items() if v is not None})

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added\n"
            f"    {key_obsm!r}, tSNE coordinates (adata.obsm)\n"
            f"    {key_uns!r}, tSNE parameters (adata.uns)"
        ),
    )

    return adata if copy else None
