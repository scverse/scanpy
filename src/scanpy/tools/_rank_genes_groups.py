"""Rank genes according to differential expression."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scverse_misc import Deprecation, deprecated_arg

from .. import _utils
from .. import logging as logg
from .._compat import CSBase, warn
from .._docs import doc_mask
from .._settings import Default, Preset, settings
from .._settings.presets import DETest
from .._utils import (
    _doc_params,
    _numba_thread_limit,
    check_nonnegative_integers,
    get_literal_vals,
)
from ..get import _check_mask, _get_arr
from ..get.get import _mask_arg
from .markers._comparison import _Comparison, _expm1_func
from .markers._results import _group_results, _legacy_stats_frame
from .markers._scorers import _illico, _logreg, _t_test, _wilcoxon
from .markers._stats import _group_stats

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import NDArray

    from ..get.get import Mask
    from .markers._results import _CorrMethod
    from .markers._stats import _GroupStats


def _legacy_compute(
    cmp: _Comparison,
    method: DETest,
    *,
    pts: bool,
    corr_method: _CorrMethod,
    n_genes_user: int | None,
    rankby_abs: bool,
    tie_correct: bool,
    mean_in_log_space: bool,
    **kwds,
) -> tuple[_GroupStats | None, pd.DataFrame | None]:
    stats = None
    if method in {"t-test", "t-test_overestim_var"}:
        stats = _group_stats(
            cmp,
            exponentiate_values=not mean_in_log_space,
            need_var=True,
            comp_pts=pts,
        )
        results = _t_test(cmp, stats, overestimate_var=method == "t-test_overestim_var")
    elif method == "wilcoxon_illico":
        results = _illico(cmp, tie_correct=tie_correct)
        stats = _group_stats(
            cmp, exponentiate_values=not mean_in_log_space, comp_pts=pts
        )
    elif method == "wilcoxon":
        results = _wilcoxon(cmp, tie_correct=tie_correct)
        stats = _group_stats(
            cmp, exponentiate_values=not mean_in_log_space, comp_pts=pts
        )
    else:
        results = _logreg(cmp, **kwds)

    group_results = _group_results(
        cmp,
        stats,
        results,
        corr_method=corr_method,
        mean_in_log_space=mean_in_log_space,
    )
    frame = _legacy_stats_frame(
        cmp, group_results, n_genes_user=n_genes_user, rankby_abs=rankby_abs
    )
    return stats, frame


@_doc_params(
    mask=doc_mask("Select subset of genes to use in statistical tests.", dim="var")
)
@deprecated_arg("mask_var", Deprecation("1.13.0", "Use `mask` instead."))
def rank_genes_groups(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    groupby: str,
    *,
    mask: Mask | None = None,
    use_raw: bool | None = None,
    groups: Literal["all"] | Iterable[str] = "all",
    reference: str = "rest",
    n_genes: int | None = None,
    rankby_abs: bool = False,
    pts: bool = False,
    key_added: str | None = None,
    copy: bool = False,
    method: DETest | Default = Default(preset=("rank_genes_groups", "method")),
    corr_method: _CorrMethod = "benjamini-hochberg",
    tie_correct: bool = False,
    layer: str | None = None,
    mean_in_log_space: bool | Default = Default(
        preset=("rank_genes_groups", "mean_in_log_space")
    ),
    mask_var: Mask | None = None,
    **kwds,
) -> AnnData | None:
    r"""Rank genes for characterizing groups.

    Expects logarithmized data.

    .. array-support:: tl.rank_genes_groups

    ..  warning::

        Comparing between cells leads to highly inflated p-values,
        since cells are not independent observations :cite:p:`Squair2021`.
        Especially in single-cell data, consider instead to use more appropriate methods such as combining pseudobulking with :doc:`pydeseq2:index`.

        :func:`decoupler.pp.pseudobulk` or :func:`scanpy.get.aggregate` can be used to aggregate samples for pseudobulking.
        Ours is a bit more verbose, but supports :doc:`dask:index` arrays for improved performance.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    {mask}
    use_raw
        Use `raw` attribute of `adata` if present. The default behavior is to use `raw` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    groups
        Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
        shall be restricted, or `'all'` (default), for all groups. Note that if
        `reference='rest'` all groups will still be used as the reference, not
        just those specified in `groups`.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    n_genes
        The number of genes that appear in the returned tables.
        Defaults to all genes.
    method
        The default method is `'t-test'`,
        `'t-test_overestim_var'` overestimates variance of each group,
        `'wilcoxon'` uses Wilcoxon rank-sum,
        `'logreg'` uses logistic regression. See :cite:t:`Ntranos2019`,
        `here <https://github.com/scverse/scanpy/issues/95>`__ and `here
        <https://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,
        for why this is meaningful.
    corr_method
        p-value correction method.
        Used only for `'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'`.
    tie_correct
        Use tie correction for `'wilcoxon'` scores.
        Used only for `'wilcoxon'`.
    rankby_abs
        Rank genes by the absolute value of the score, not by the
        score. The returned scores are never the absolute values.
    pts
        Compute the fraction of cells expressing the genes.
    key_added
        The key in `adata.uns` information is saved to.
    copy
        Whether to copy `adata` or modify it inplace.
    mean_in_log_space
        Whether to do :math:`\log(\operatorname{{mean}}(e^x))` (`False`)
        or :math:`\log(e^{{\operatorname{{mean}}(x)}})` (`True`).
        The former is accurate, while the latter is a faster approximation
        that underestimates this accurate result in the presence of many outliers.
    kwds
        Are passed to test methods. Currently this affects only parameters that
        are passed to :class:`sklearn.linear_model.LogisticRegression`.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.uns['rank_genes_groups' | key_added]['names']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['scores']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the z-score
        underlying the computation of a p-value for each gene for each
        group. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['logfoldchanges']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the log2
        fold change for each gene for each group. Ordered according to
        scores. Only provided if method is 't-test' like.
        Note: if `mean_in_log_space=True`, this is an approximation calculated from mean-log values.
    `adata.uns['rank_genes_groups' | key_added]['pvals']` : structured :class:`numpy.ndarray` (dtype `float`)
        p-values.
    `adata.uns['rank_genes_groups' | key_added]['pvals_adj']` : structured :class:`numpy.ndarray` (dtype `float`)
        Corrected p-values.
    `adata.uns['rank_genes_groups' | key_added]['pts']` : :class:`pandas.DataFrame` (dtype `float`)
        Fraction of cells expressing the genes for each group.
    `adata.uns['rank_genes_groups' | key_added]['pts_rest']` : :class:`pandas.DataFrame` (dtype `float`)
        Only if `reference` is set to `'rest'`.
        Fraction of cells from the union of the rest of each group
        expressing the genes.

    Notes
    -----
    There are slight inconsistencies depending on whether sparse
    or dense data are passed. See `here <https://github.com/scverse/scanpy/blob/main/tests/test_rank_genes_groups.py>`__.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")
    >>> # to visualize the results
    >>> sc.pl.rank_genes_groups(adata)

    """
    mask = _mask_arg(mask, mask_var, dim="var")
    if isinstance(mean_in_log_space, Default):
        mean_in_log_space = settings.preset.rank_genes_groups.mean_in_log_space
    # If scanpy presets are used for v2, use illico - prevents the presets from showing the `wilcoxon_illico` method and allows us to silently replace `wilcoxon`'s implementation.
    if method is None or isinstance(method, Default):
        method = settings.preset.rank_genes_groups.method
        if settings.preset is Preset.ScanpyV2Preview:
            method = "wilcoxon_illico"
    # Otherwise, nudge people to use the presets.
    elif "illico" in method:
        msg = (
            "`wilcoxon_illico` flavor will be removed in scanpy 2.0 and be simply the new `wilcoxon` implementation."
            "To remove theis warning, you can locally do `with sc.settings.override(preset=sc.Preset.ScanpyV2Preview)`."
        )
        warn(msg, DeprecationWarning)

    mask_var = _check_mask(adata, mask, "var")

    if use_raw is None:
        use_raw = adata.raw is not None
    elif use_raw is True and adata.raw is None:
        msg = "Received `use_raw=True`, but `adata.raw` is empty."
        raise ValueError(msg)

    if "only_positive" in kwds:
        rankby_abs = not kwds.pop("only_positive")  # backwards compat

    start = logg.info("ranking genes")
    if method not in (avail_methods := get_literal_vals(DETest)):
        msg = f"Method must be one of {avail_methods}."
        raise ValueError(msg)

    avail_corr = {"benjamini-hochberg", "bonferroni"}
    if corr_method not in avail_corr:
        msg = f"Correction method must be one of {avail_corr}."
        raise ValueError(msg)

    adata = adata.copy() if copy else adata
    _utils.sanitize_anndata(adata)
    # for clarity, rename variable
    if groups == "all":
        groups_order = "all"
    elif isinstance(groups, str | int):
        msg = "Specify a sequence of groups"
        raise ValueError(msg)
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if reference != "rest" and reference not in set(groups_order):
            groups_order += [reference]
    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        cats = adata.obs[groupby].cat.categories.tolist()
        msg = f"reference = {reference} needs to be one of groupby = {cats}."
        raise ValueError(msg)

    if key_added is None:
        key_added = "rank_genes_groups"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = dict(
        groupby=groupby,
        reference=reference,
        method=method,
        use_raw=use_raw,
        layer=layer,
        corr_method=corr_method,
    )

    comparison = _Comparison.legacy(
        adata,
        groups_order,
        groupby,
        mask_var=mask_var,
        reference=reference,
        use_raw=use_raw,
        layer=layer,
    )

    if check_nonnegative_integers(comparison.x) and method != "logreg":
        logg.warning(
            "It seems you use rank_genes_groups on the raw count data. "
            "Please logarithmize your data before calling rank_genes_groups."
        )

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    # defaults to all genes
    if n_genes_user is None or n_genes_user > comparison.x.shape[1]:
        n_genes_user = comparison.x.shape[1]

    logg.debug(f"consider {groupby!r} groups:")
    logg.debug(f"with sizes: {np.count_nonzero(comparison.groups_masks_obs, axis=1)}")

    with _numba_thread_limit(settings.n_jobs if method == "wilcoxon" else None):
        stats, frame = _legacy_compute(
            comparison,
            method,
            pts=pts,
            corr_method=corr_method,
            n_genes_user=n_genes_user,
            rankby_abs=rankby_abs,
            tie_correct=tie_correct,
            mean_in_log_space=mean_in_log_space,
            **kwds,
        )

    groups_names = [str(name) for name in comparison.groups_order]
    if stats is not None and stats.pts is not None:
        adata.uns[key_added]["pts"] = pd.DataFrame(
            stats.pts.T, index=comparison.var_names, columns=groups_names
        )
    if stats is not None and stats.pts_rest is not None:
        adata.uns[key_added]["pts_rest"] = pd.DataFrame(
            stats.pts_rest.T, index=comparison.var_names, columns=groups_names
        )

    frame.columns = frame.columns.swaplevel()

    dtypes = {
        "names": "O",
        "scores": "float32",
        "logfoldchanges": "float32",
        "pvals": "float64",
        "pvals_adj": "float64",
    }

    for col in frame.columns.levels[0]:
        adata.uns[key_added][col] = frame[col].to_records(
            index=False, column_dtypes=dtypes[col]
        )

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added to `.uns[{key_added!r}]`\n"
            "    'names', sorted np.recarray to be indexed by group ids\n"
            "    'scores', sorted np.recarray to be indexed by group ids\n"
            + (
                "    'logfoldchanges', sorted np.recarray to be indexed by group ids\n"
                "    'pvals', sorted np.recarray to be indexed by group ids\n"
                "    'pvals_adj', sorted np.recarray to be indexed by group ids"
                if method in {"t-test", "t-test_overestim_var", "wilcoxon"}
                else ""
            )
        ),
    )
    return adata if copy else None


def _calc_frac(x: NDArray[np.number] | CSBase, /) -> NDArray[np.float64]:
    n_nonzero = (
        x.getnnz(axis=0) if isinstance(x, CSBase) else np.count_nonzero(x, axis=0)
    )
    return n_nonzero / x.shape[0]


def filter_rank_genes_groups(
    adata: AnnData,
    *,
    key: str | None = None,
    groupby: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    key_added: str = "rank_genes_groups_filtered",
    min_in_group_fraction: float = 0.25,
    min_fold_change: float = 1,
    max_out_group_fraction: float = 0.5,
    compare_abs: bool = False,
) -> None:
    """Filter out genes based on two criteria.

    1. log fold change and
    2. fraction of genes expressing the
       gene within and outside the `groupby` categories.

    See :func:`~scanpy.tl.rank_genes_groups`.

    Results are stored in `adata.uns[key_added]`
    (default: 'rank_genes_groups_filtered').

    To preserve the original structure of adata.uns['rank_genes_groups'],
    filtered genes are set to `NaN`.

    Parameters
    ----------
    adata
    key
    groupby
    layer
    use_raw
    key_added
    min_in_group_fraction
    min_fold_change
    max_out_group_fraction
    compare_abs
        If `True`, compare absolute values of log fold change with `min_fold_change`.

    Returns
    -------
    Same output as :func:`scanpy.tl.rank_genes_groups` but with filtered genes names set to `nan`.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")
    >>> sc.tl.filter_rank_genes_groups(adata, min_fold_change=3)
    >>> # visualize results
    >>> sc.pl.rank_genes_groups(adata, key="rank_genes_groups_filtered")
    >>> # visualize results using dotplot
    >>> sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered")

    """
    if key is None:
        key = "rank_genes_groups"

    if groupby is None:
        groupby = adata.uns[key]["params"]["groupby"]

    if use_raw is None:
        use_raw = adata.uns[key]["params"]["use_raw"] if layer is None else False

    x = _get_arr(adata, use_raw=use_raw, layer=layer)

    same_params = (
        adata.uns[key]["params"]["groupby"] == groupby
        and adata.uns[key]["params"]["reference"] == "rest"
        and adata.uns[key]["params"]["use_raw"] == use_raw
    )

    use_logfolds = same_params and "logfoldchanges" in adata.uns[key]
    use_fraction = same_params and "pts_rest" in adata.uns[key]

    # convert structured numpy array into DataFrame
    gene_names = pd.DataFrame(adata.uns[key]["names"])

    fraction_in_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )
    fraction_out_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )

    if use_logfolds:
        fold_change_matrix = pd.DataFrame(adata.uns[key]["logfoldchanges"])
    else:
        fold_change_matrix = pd.DataFrame(
            np.zeros(gene_names.shape),
            columns=gene_names.columns,
            index=gene_names.index,
        )

        expm1_func = _expm1_func(adata)

    logg.info(
        f"Filtering genes using: "
        f"min_in_group_fraction: {min_in_group_fraction} "
        f"min_fold_change: {min_fold_change}, "
        f"max_out_group_fraction: {max_out_group_fraction}"
    )

    for cluster in gene_names.columns:
        # iterate per column
        var_names = gene_names[cluster].array

        if not use_logfolds or not use_fraction:
            var_idx = (adata.raw if use_raw else adata).var_names.get_indexer(var_names)
            sub_x = x[:, var_idx]
            in_group = (adata.obs[groupby] == cluster).to_numpy()
            x_in = sub_x[in_group]
            x_out = sub_x[~in_group]

        if use_fraction:
            fraction_in_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]["pts"][cluster].loc[var_names].array
            )
            fraction_out_cluster_matrix.loc[:, cluster] = (
                adata.uns[key]["pts_rest"][cluster].loc[var_names].array
            )
        else:
            fraction_in_cluster_matrix.loc[:, cluster] = _calc_frac(x_in)
            fraction_out_cluster_matrix.loc[:, cluster] = _calc_frac(x_out)

        if not use_logfolds:
            # compute mean value
            mean_in_cluster = np.ravel(x_in.mean(0))
            mean_out_cluster = np.ravel(x_out.mean(0))
            # compute fold change
            fold_change_matrix.loc[:, cluster] = np.log2(
                (expm1_func(mean_in_cluster) + 1e-9)
                / (expm1_func(mean_out_cluster) + 1e-9)
            )

    if compare_abs:
        fold_change_matrix = fold_change_matrix.abs()
    # filter original_matrix
    gene_names = gene_names[
        (fraction_in_cluster_matrix > min_in_group_fraction)
        & (fraction_out_cluster_matrix < max_out_group_fraction)
        & (fold_change_matrix > min_fold_change)
    ]
    # create new structured array using 'key_added'.
    adata.uns[key_added] = adata.uns[key].copy()
    adata.uns[key_added]["names"] = gene_names.to_records(index=False)
