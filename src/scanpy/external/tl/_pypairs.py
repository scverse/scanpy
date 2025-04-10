"""Calculate scores based on relative expression change of maker pairs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from packaging.version import Version

from ..._settings import settings
from ..._utils._doctests import doctest_needs

if TYPE_CHECKING:
    from collections.abc import Collection, Mapping

    import pandas as pd
    from anndata import AnnData

    Genes = Collection[str | int | bool]


@doctest_needs("pypairs")
def sandbag(
    adata: AnnData,
    annotation: Mapping[str, Genes] | None = None,
    *,
    fraction: float = 0.65,
    filter_genes: Genes | None = None,
    filter_samples: Genes | None = None,
) -> dict[str, list[tuple[str, str]]]:
    """Calculate marker pairs of genes :cite:p:`Scialdone2015,Fechtner2018`.

    Calculates the pairs of genes serving as marker pairs for each phase,
    based on a matrix of gene counts and an annotation of known phases.

    This reproduces the approach of :cite:t:`Scialdone2015` in the implementation of
    :cite:t:`Fechtner2018`.

    More information and bug reports `here
    <https://github.com/rfechtner/pypairs>`__.

    Parameters
    ----------
    adata
        The annotated data matrix.
    annotation
        Mapping from category to genes, e.g. `{'phase': [Gene1, ...]}`.
        Defaults to ``data.vars['category']``.
    fraction
        Fraction of cells per category where marker criteria must be satisfied.
    filter_genes
        Genes for sampling the reference set. Defaults to all genes.
    filter_samples
        Cells for sampling the reference set. Defaults to all samples.

    Returns
    -------
    A dict mapping from category to lists of marker pairs, e.g.:
    `{'Category_1': [(Gene_1, Gene_2), ...], ...}`.

    Examples
    --------
    >>> from scanpy.external.tl import sandbag
    >>> from pypairs import datasets
    >>> adata = datasets.leng15()
    >>> marker_pairs = sandbag(adata, fraction=0.5)

    """
    _check_import()
    from pypairs import settings as pp_settings
    from pypairs.pairs import sandbag

    pp_settings.verbosity = settings.verbosity
    pp_settings.n_jobs = settings.n_jobs
    pp_settings.writedir = settings.writedir
    pp_settings.cachedir = settings.cachedir
    pp_settings.logfile = settings.logfile

    return sandbag(
        data=adata,
        annotation=annotation,
        fraction=fraction,
        filter_genes=filter_genes,
        filter_samples=filter_samples,
    )


def cyclone(
    adata: AnnData,
    marker_pairs: Mapping[str, Collection[tuple[str, str]]] | None = None,
    *,
    iterations: int = 1000,
    min_iter: int = 100,
    min_pairs: int = 50,
) -> pd.DataFrame:
    """Assign scores and predicted class to observations :cite:p:`Scialdone2015` :cite:p:`Fechtner2018`.

    Calculates scores for each observation and each phase and assigns prediction
    based on marker pairs indentified by :func:`~scanpy.external.tl.sandbag`.

    This reproduces the approach of :cite:t:`Scialdone2015` in the implementation of
    :cite:t:`Fechtner2018`.

    Parameters
    ----------
    adata
        The annotated data matrix.
    marker_pairs
        Mapping of categories to lists of marker pairs.
        See :func:`~scanpy.external.tl.sandbag` output.
    iterations
        An integer scalar specifying the number of
        iterations for random sampling to obtain a cycle score.
    min_iter
        An integer scalar specifying the minimum number of iterations
        for score estimation.
    min_pairs
        An integer scalar specifying the minimum number of pairs
        for score estimation.

    Returns
    -------
    A :class:`~pandas.DataFrame` with samples as index and categories as columns
    with scores for each category for each sample and a additional column with
    the name of the max scoring category for each sample.

    If `marker_pairs` contains only the cell cycle categories G1, S and G2M an
    additional column `pypairs_cc_prediction` will be added.
    Where category S is assigned to samples where G1 and G2M score are < 0.5.

    """
    _check_import()
    from pypairs import settings as pp_settings
    from pypairs.pairs import cyclone

    pp_settings.verbosity = settings.verbosity
    pp_settings.n_jobs = settings.n_jobs
    pp_settings.writedir = settings.writedir
    pp_settings.cachedir = settings.cachedir
    pp_settings.logfile = settings.logfile

    return cyclone(
        data=adata,
        marker_pairs=marker_pairs,
        iterations=iterations,
        min_iter=min_iter,
        min_pairs=min_pairs,
    )


def _check_import():
    try:
        import pypairs
    except ImportError as e:
        msg = "You need to install the package `pypairs`."
        raise ImportError(msg) from e

    min_version = Version("3.0.9")
    if Version(pypairs.__version__) < min_version:
        msg = f"Please only use `pypairs` >= {min_version}"
        raise ImportError(msg)
