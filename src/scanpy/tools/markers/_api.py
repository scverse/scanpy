from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

from ... import logging as logg
from ..._docs import doc_mask
from ..._settings import Default, settings
from ..._utils import _doc_params, _numba_thread_limit, check_nonnegative_integers
from ._comparison import _Comparison
from ._results import _group_results, _marker_table
from ._scorers import _illico, _logreg_signed, _t_test, _wilcoxon
from ._stats import _group_stats

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Literal

    import pandas as pd
    from anndata import AnnData

    from ..._settings.presets import WilcoxonBackend
    from ...get.get import Mask
    from ._results import _CorrMethod
    from ._scorers import _TestResult
    from ._stats import _GroupStats


doc_marker_params = """\
adata
    Annotated data matrix. Expects logarithmized data.
groupby
    Key in `.obs` of the grouping (e.g. clusters) to find marker genes for.
groups
    Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which the comparison shall be restricted,
    or `'all'` (default), for all groups.
    If `reference='rest'`, cells of all other groups still form each group’s rest.
reference
    If `'rest'`, compare each group to the union of all other cells.
    If a group identifier, compare each group to this group.
layer
    Key from `adata.layers` to use instead of `adata.X`."""

doc_marker_params_tail = """\
n_top
    Number of top-ranked genes per group to return. Defaults to all genes.
    All genes are scored; this only limits the size of the returned table.
corr_method
    p-value correction method, applied per group across genes."""

doc_mean_in_log_space = r"""mean_in_log_space
    Whether to compute the log fold change as :math:`\log_2(e^{\operatorname{mean}(x)})` (`True`, faster)
    or as :math:`\log_2(\operatorname{mean}(e^x))` (`False`, accurate)."""

doc_returns = """\
:class:`pandas.DataFrame` with one row per group and gene, ordered by descending `score` within each group:

`group`, `reference`
    The group and what it was compared to (`'rest'` or the reference group).
`gene`
    The gene (`adata.var_names`).
`score`
    The method’s test statistic (see above).
`log_fc`
    log2 fold change of the group mean versus the reference mean.
`p_value`, `adj_p_value`
    Raw and multiple-testing corrected p-values.
`frac_group`, `frac_reference`
    Fraction of cells in the group / reference with non-zero expression."""

doc_note = """\
These are fast cell-level tests that rank genes distinguishing groups of cells, e.g. to annotate clusters.
Cells are treated as independent observations, so the p-values are only useful to order genes
:cite:p:`Squair2021`.
To test for differential expression between conditions, use :func:`scanpy.tl.de.deseq2`."""


def _run_markers(  # noqa: PLR0913
    adata: AnnData,
    groupby: str,
    scorer: Callable[[_Comparison, _GroupStats | None], Iterable[_TestResult]],
    *,
    groups: Literal["all"] | Iterable[str],
    reference: str,
    layer: str | None,
    mask: Mask | None,
    n_top: int | None,
    corr_method: _CorrMethod,
    mean_in_log_space: bool | Default,
    test_needs_var: bool,
    func_name: str,
) -> pd.DataFrame:
    if isinstance(mean_in_log_space, Default):
        mean_in_log_space = settings.preset.markers.mean_in_log_space
    start = logg.info("ranking marker genes")
    cmp = _Comparison.from_adata(
        adata,
        groupby,
        groups=groups,
        reference=reference,
        layer=layer,
        mask=mask,
        func_name=func_name,
    )
    if func_name != "logreg" and check_nonnegative_integers(cmp.x):
        logg.warning(
            f"It seems you use {func_name} on raw count data. "
            "Please logarithmize your data first."
        )
    with _numba_thread_limit(settings.n_jobs):
        if test_needs_var:
            test_stats = _group_stats(cmp, need_var=True, comp_pts=True)
            lfc_stats = (
                test_stats
                if mean_in_log_space
                else _group_stats(cmp, exponentiate_values=True)
            )
        else:
            lfc_stats = test_stats = _group_stats(
                cmp, exponentiate_values=not mean_in_log_space, comp_pts=True
            )
        group_results = _group_results(
            cmp,
            lfc_stats,
            scorer(cmp, test_stats),
            corr_method=corr_method,
            mean_in_log_space=mean_in_log_space,
        )
        df = _marker_table(cmp, test_stats, group_results, n_top=n_top)
    logg.info("    finished", time=start)
    return df


@_doc_params(
    params=doc_marker_params,
    mask=doc_mask("Select subset of genes to test.", dim="var"),
    params_tail=doc_marker_params_tail,
    mean_in_log_space=doc_mean_in_log_space,
    returns=doc_returns,
    note=doc_note,
)
def wilcoxon(  # noqa: PLR0913
    adata: AnnData,
    groupby: str,
    *,
    groups: Literal["all"] | Iterable[str] = "all",
    reference: str = "rest",
    layer: str | None = None,
    mask: Mask | None = None,
    n_top: int | None = None,
    corr_method: _CorrMethod = "benjamini-hochberg",
    mean_in_log_space: bool | Default = Default(
        preset=("markers", "mean_in_log_space")
    ),
    tie_correct: bool = False,
    backend: WilcoxonBackend | Default = Default(
        preset=("markers", "wilcoxon_backend")
    ),
) -> pd.DataFrame:
    """Rank marker genes per group with a Wilcoxon rank-sum test.

    .. array-support:: tl.markers.wilcoxon

    {note}

    The `score` is the z-score underlying the test’s normal approximation.

    Parameters
    ----------
    {params}
    {mask}
    {params_tail}
    {mean_in_log_space}
    tie_correct
        Use tie correction for the z-scores.
    backend
        Implementation: scanpy’s own (`'numba'`) or `illico <https://github.com/remydubois/illico>`__ (`'illico'`).

    Returns
    -------
    {returns}

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced().raw.to_adata()
    >>> markers = sc.tl.markers.wilcoxon(adata, "bulk_labels", n_top=5)

    """
    if isinstance(backend, Default):
        backend = settings.preset.markers.wilcoxon_backend
    match backend:
        case "numba":
            scorer = partial(_wilcoxon_scorer, tie_correct=tie_correct)
        case "illico":
            scorer = partial(_illico_scorer, tie_correct=tie_correct)
        case _:
            msg = f"`backend` must be 'numba' or 'illico', not {backend!r}."
            raise ValueError(msg)
    return _run_markers(
        adata,
        groupby,
        scorer,
        groups=groups,
        reference=reference,
        layer=layer,
        mask=mask,
        n_top=n_top,
        corr_method=corr_method,
        mean_in_log_space=mean_in_log_space,
        test_needs_var=False,
        func_name="wilcoxon",
    )


@_doc_params(
    params=doc_marker_params,
    mask=doc_mask("Select subset of genes to test.", dim="var"),
    params_tail=doc_marker_params_tail,
    mean_in_log_space=doc_mean_in_log_space,
    returns=doc_returns,
    note=doc_note,
)
def ttest(
    adata: AnnData,
    groupby: str,
    *,
    groups: Literal["all"] | Iterable[str] = "all",
    reference: str = "rest",
    layer: str | None = None,
    mask: Mask | None = None,
    n_top: int | None = None,
    corr_method: _CorrMethod = "benjamini-hochberg",
    mean_in_log_space: bool | Default = Default(
        preset=("markers", "mean_in_log_space")
    ),
    overestimate_var: bool = False,
) -> pd.DataFrame:
    """Rank marker genes per group with Welch’s t-test.

    .. array-support:: tl.markers.ttest

    {note}

    The `score` is the t-statistic, computed on the data as given.

    Parameters
    ----------
    {params}
    {mask}
    {params_tail}
    {mean_in_log_space}
        Only affects `log_fc`.
    overestimate_var
        Overestimate the variance of each group, by using the group’s size for the reference too.

    Returns
    -------
    {returns}

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced().raw.to_adata()
    >>> markers = sc.tl.markers.ttest(adata, "bulk_labels", n_top=5)

    """
    return _run_markers(
        adata,
        groupby,
        partial(_t_test_scorer, overestimate_var=overestimate_var),
        groups=groups,
        reference=reference,
        layer=layer,
        mask=mask,
        n_top=n_top,
        corr_method=corr_method,
        mean_in_log_space=mean_in_log_space,
        test_needs_var=True,
        func_name="ttest",
    )


@_doc_params(
    mask=doc_mask("Select subset of genes to use as features.", dim="var"),
    mean_in_log_space=doc_mean_in_log_space,
)
def logreg(
    adata: AnnData,
    groupby: str,
    *,
    groups: Literal["all"] | Iterable[str] = "all",
    layer: str | None = None,
    mask: Mask | None = None,
    n_top: int | None = None,
    mean_in_log_space: bool | Default = Default(
        preset=("markers", "mean_in_log_space")
    ),
    **kwds,
) -> pd.DataFrame:
    """Rank marker genes per group by logistic regression coefficients :cite:p:`Ntranos2019`.

    .. array-support:: tl.markers.logreg

    Fits one multi-variate classifier predicting `groupby` from all genes,
    and uses each group’s coefficients as `score`.
    There are no p-values.

    Parameters
    ----------
    adata
        Annotated data matrix. Expects logarithmized data.
    groupby
        Key in `.obs` of the grouping (e.g. clusters) to find marker genes for.
    groups
        Subset of groups to fit the classifier on, or `'all'` (default).
    layer
        Key from `adata.layers` to use instead of `adata.X`.
    {mask}
    n_top
        Number of top-ranked genes per group to return. Defaults to all genes.
    {mean_in_log_space}
        Only affects `log_fc`.
    kwds
        Passed to :class:`sklearn.linear_model.LogisticRegression`,
        e.g. `penalty='l1'` for a sparse solution.

    Returns
    -------
    :class:`pandas.DataFrame` with the same columns as :func:`~scanpy.tl.markers.wilcoxon`,
    without `p_value` and `adj_p_value`.
    `log_fc` and `frac_reference` compare each group to all other cells.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced().raw.to_adata()
    >>> markers = sc.tl.markers.logreg(adata, "bulk_labels", n_top=5, max_iter=1000)

    """
    return _run_markers(
        adata,
        groupby,
        partial(_logreg_scorer, **kwds),
        groups=groups,
        reference="rest",
        layer=layer,
        mask=mask,
        n_top=n_top,
        corr_method="benjamini-hochberg",
        mean_in_log_space=mean_in_log_space,
        test_needs_var=False,
        func_name="logreg",
    )


def _wilcoxon_scorer(
    cmp: _Comparison, _stats: _GroupStats | None, *, tie_correct: bool
) -> Iterable[_TestResult]:
    return _wilcoxon(cmp, tie_correct=tie_correct)


def _illico_scorer(
    cmp: _Comparison, _stats: _GroupStats | None, *, tie_correct: bool
) -> Iterable[_TestResult]:
    return _illico(cmp, tie_correct=tie_correct)


def _t_test_scorer(
    cmp: _Comparison, stats: _GroupStats | None, *, overestimate_var: bool
) -> Iterable[_TestResult]:
    return _t_test(cmp, stats, overestimate_var=overestimate_var)


def _logreg_scorer(
    cmp: _Comparison, _stats: _GroupStats | None, **kwds
) -> Iterable[_TestResult]:
    return _logreg_signed(cmp, **kwds)
