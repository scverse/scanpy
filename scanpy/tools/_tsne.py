from packaging import version
from typing import Optional, Union
import warnings

from anndata import AnnData

from .._utils import _doc_params, AnyRandom
from ..tools._utils import _choose_representation, doc_use_rep, doc_n_pcs
from .._settings import settings
from .. import logging as logg


@_doc_params(doc_n_pcs=doc_n_pcs, use_rep=doc_use_rep)
def tsne(
    adata: AnnData,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    perplexity: Union[float, int] = 30,
    early_exaggeration: Union[float, int] = 12,
    learning_rate: Union[float, int] = 1000,
    random_state: AnyRandom = 0,
    use_fast_tsne: bool = False,
    n_jobs: Optional[int] = None,
    copy: bool = False,
    *,
    metric: str = "euclidean",
) -> Optional[AnnData]:
    """\
    t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.

    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been
    proposed for visualizating single-cell data by [Amir13]_. Here, by default,
    we use the implementation of *scikit-learn* [Pedregosa11]_. You can achieve
    a huge speedup and better convergence if you install `Multicore-tSNE
    <https://github.com/DmitryUlyanov/Multicore-TSNE>`__ by [Ulyanov16]_, which
    will be automatically detected by Scanpy.

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
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_tsne** : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    """
    import sklearn

    start = logg.info('computing tSNE')
    adata = adata.copy() if copy else adata
    X = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)
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
    # square_distances will default to true in the future, we'll get ahead of the
    # warning for now
    if metric != "euclidean":
        sklearn_version = version.parse(sklearn.__version__)
        if sklearn_version >= version.parse("0.24.0"):
            params_sklearn["square_distances"] = True
        else:
            warnings.warn(
                "Results for non-euclidean metrics changed in sklearn 0.24.0, while "
                f"you are using {sklearn.__version__}.",
                UserWarning,
            )

    # Backwards compat handling: Remove in scanpy 1.9.0
    if n_jobs != 1 and not use_fast_tsne:
        warnings.warn(
            UserWarning(
                "In previous versions of scanpy, calling tsne with n_jobs > 1 would use "
                "MulticoreTSNE. Now this uses the scikit-learn version of TSNE by default. "
                "If you'd like the old behaviour (which is deprecated), pass "
                "'use_fast_tsne=True'. Note, MulticoreTSNE is not actually faster anymore."
            )
        )
    if use_fast_tsne:
        warnings.warn(
            FutureWarning(
                "Argument `use_fast_tsne` is deprecated, and support for MulticoreTSNE "
                "will be dropped in a future version of scanpy."
            )
        )

    # deal with different tSNE implementations
    if use_fast_tsne:
        try:
            from MulticoreTSNE import MulticoreTSNE as TSNE

            tsne = TSNE(**params_sklearn)
            logg.info("    using the 'MulticoreTSNE' package by Ulyanov (2017)")
            # need to transform to float64 for MulticoreTSNE...
            X_tsne = tsne.fit_transform(X.astype('float64'))
        except ImportError:
            use_fast_tsne = False
            warnings.warn(
                UserWarning(
                    "Could not import 'MulticoreTSNE'. Falling back to scikit-learn."
                )
            )
    if use_fast_tsne is False:  # In case MultiCore failed to import
        from sklearn.manifold import TSNE

        # unfortunately, sklearn does not allow to set a minimum number
        # of iterations for barnes-hut tSNE
        tsne = TSNE(**params_sklearn)
        logg.info('    using sklearn.manifold.TSNE')
        X_tsne = tsne.fit_transform(X)

    # update AnnData instance
    adata.obsm['X_tsne'] = X_tsne  # annotate samples with tSNE coordinates
    adata.uns["tsne"] = {
        "params": {
            k: v
            for k, v in {
                "perplexity": perplexity,
                "early_exaggeration": early_exaggeration,
                "learning_rate": learning_rate,
                "n_jobs": n_jobs,
                "metric": metric,
                "use_rep": use_rep,
            }.items()
            if v is not None
        }
    }

    logg.info(
        '    finished',
        time=start,
        deep="added\n    'X_tsne', tSNE coordinates (adata.obsm)",
    )

    return adata if copy else None
