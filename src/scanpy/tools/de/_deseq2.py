from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from anndata import AnnData

from ... import logging as logg
from ..._compat import CSBase
from ..._settings import settings
from ..._utils import (
    check_nonnegative_integers,
    raise_not_implemented_error_if_backed_type,
)
from ..._utils._doctests import doctest_needs
from ...get import aggregate

if TYPE_CHECKING:
    from collections.abc import Sequence


@doctest_needs("pydeseq2")
def deseq2(
    adata: AnnData,
    condition: str,
    *,
    reference: str,
    sample_key: str,
    groupby: str | None = None,
    covariates: Sequence[str] = (),
    layer: str | None = None,
    min_cells: int = 10,
    n_jobs: int | None = None,
) -> pd.DataFrame:
    """Test for differential expression between conditions with pseudobulk DESeq2.

    Use this to find genes whose expression changes between conditions
    (e.g. treated vs. control, disease vs. healthy), using biological replicates.
    Counts are summed per sample (pseudobulk) with :func:`scanpy.get.aggregate`,
    then tested with :doc:`PyDESeq2 <pydeseq2:index>` :cite:p:`Love2014,Muzellec2023`.
    Treating cells instead of samples as replicates inflates p-values :cite:p:`Squair2021`.

    To find genes that distinguish clusters instead, use :mod:`scanpy.tl.markers`.

    Requires `pydeseq2`: `pip install 'scanpy[pydeseq2]'`.
    For more complex designs (interactions, mixed models, other methods), use PyDESeq2 or
    `pertpy <https://pertpy.readthedocs.io>`__ directly.

    Parameters
    ----------
    adata
        Annotated data matrix with raw counts in `adata.X` or `adata.layers[layer]`.
    condition
        Key in `.obs` of the condition to test, e.g. `'treatment'`.
    reference
        The baseline level of `condition`, e.g. `'control'`.
        Every other level is compared to it.
    sample_key
        Key in `.obs` of the biological replicate (e.g. donor or sample) each cell comes from.
    groupby
        Key in `.obs` of a cell grouping (e.g. cell type) to test separately.
        By default, all cells are tested together.
    covariates
        Keys in `.obs` of sample-level covariates to adjust for (e.g. batch).
        Each must be constant within a sample.
        To compare conditions within each donor (paired design), pass `sample_key` here
        and a `sample_key` that is shared between conditions.
    layer
        Key in `.layers` containing raw counts. Defaults to `adata.X`.
    min_cells
        Pseudobulk samples made of fewer cells are dropped.
    n_jobs
        Number of CPUs to use. Defaults to :attr:`scanpy.settings.n_jobs`.

    Returns
    -------
    :class:`pandas.DataFrame` with one row per tested group, condition and gene:

    `group`
        The `groupby` value (only if `groupby` is given).
    `condition`, `reference`
        The tested level of `condition` and the baseline it is compared to.
    `gene`
        The gene (`adata.var_names`).
    `base_mean`
        Mean of normalized counts across samples.
    `log_fc`, `log_fc_se`
        log2 fold change of `condition` versus `reference`, and its standard error.
    `statistic`, `p_value`, `adj_p_value`
        Wald statistic, p-value, and Benjamini–Hochberg adjusted p-value (per group and condition).

    Groups or conditions with fewer than two samples on either side are skipped with a warning.

    Examples
    --------
    >>> import scanpy as sc
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> means = rng.uniform(0.5, 50, 200)
    >>> counts = rng.negative_binomial(5, 5 / (5 + means), (400, 200))
    >>> adata = sc.AnnData(counts.astype(np.float32))
    >>> adata.obs["donor"] = np.repeat([f"d{i}" for i in range(8)], 50)
    >>> adata.obs["treatment"] = np.tile(np.repeat(["ctrl", "stim"], 25), 8)
    >>> res = sc.tl.de.deseq2(
    ...     adata, "treatment", reference="ctrl", sample_key="donor", n_jobs=1
    ... )
    >>> res.columns.tolist()
    ['condition', 'reference', 'gene', 'base_mean', 'log_fc', 'log_fc_se', 'statistic', 'p_value', 'adj_p_value']

    """
    try:
        from pydeseq2.default_inference import DefaultInference
    except ImportError as e:
        e.add_note(
            "Please install `scanpy[pydeseq2]` (or `pydeseq2` directly) and try again."
        )
        raise

    keys = [condition, sample_key, *([] if groupby is None else [groupby]), *covariates]
    if missing := [k for k in keys if k not in adata.obs.columns]:
        msg = f"Could not find {missing} in `adata.obs`."
        raise KeyError(msg)
    if reference not in set(adata.obs[condition]):
        msg = f"`reference={reference!r}` is not a value of `adata.obs[{condition!r}]`."
        raise ValueError(msg)

    x = adata.X if layer is None else adata.layers[layer]
    raise_not_implemented_error_if_backed_type(x, "de.deseq2")
    if not check_nonnegative_integers(x):
        msg = (
            "`deseq2` needs raw counts (non-negative integers). "
            "If they are in a layer, pass e.g. `layer='counts'`."
        )
        raise ValueError(msg)

    by = [sample_key, condition, *([] if groupby is None else [groupby])]
    sample_cols = [sample_key, condition]
    covariate_table = _covariate_table(adata.obs, sample_cols, covariates)
    pseudobulk = aggregate(
        AnnData(X=x, obs=adata.obs[by].copy(), var=pd.DataFrame(index=adata.var_names)),
        by=by,
        func="sum",
    )
    pseudobulk = pseudobulk[pseudobulk.obs["n_obs_aggregated"] >= min_cells]
    inference = DefaultInference(n_cpus=settings.n_jobs if n_jobs is None else n_jobs)

    groups = (
        [(None, pseudobulk)]
        if groupby is None
        else [
            (group, pseudobulk[pseudobulk.obs[groupby] == group])
            for group in pseudobulk.obs[groupby].cat.categories
        ]
    )
    frames = [
        frame
        for group, pb in groups
        for frame in _test_group(
            pb,
            group,
            condition=condition,
            reference=reference,
            covariate_table=covariate_table,
            inference=inference,
        )
    ]
    if not frames:
        msg = "No group had at least two samples per condition."
        raise ValueError(msg)
    df = pd.concat(frames, ignore_index=True)
    if groupby is None:
        df = df.drop(columns="group")
    for col in ["group", "condition", "reference"]:
        if col in df:
            df[col] = df[col].astype("category")
    return df


def _test_group(
    pseudobulk: AnnData,
    group: str | None,
    *,
    condition: str,
    reference: str,
    covariate_table: pd.DataFrame,
    inference,
) -> list[pd.DataFrame]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    label = "all cells" if group is None else f"group {group!r}"
    n_samples = pseudobulk.obs[condition].value_counts()
    if n_samples.get(reference, 0) < 2:
        logg.warning(f"Skipping {label}: fewer than two samples for {reference!r}.")
        return []
    levels = [lvl for lvl, n in n_samples.items() if lvl != reference and n >= 2]
    if not levels:
        logg.warning(f"Skipping {label}: no other condition has two samples.")
        return []
    pb = pseudobulk[pseudobulk.obs[condition].isin([reference, *levels])]

    metadata = _design_metadata(pb.obs, condition, covariate_table)
    counts = pb.layers["sum"]
    counts = counts.toarray() if isinstance(counts, CSBase) else np.asarray(counts)
    dds = DeseqDataSet(
        counts=pd.DataFrame(
            counts.astype(np.int64), index=pb.obs_names, columns=pb.var_names
        ),
        metadata=metadata.set_index(pb.obs_names),
        design="~" + " + ".join([*metadata.columns[1:], "condition"]),
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    frames = []
    for level in levels:
        stats = DeseqStats(
            dds,
            contrast=["condition", str(level), str(reference)],
            inference=inference,
            quiet=True,
        )
        stats.summary()
        frames.append(_result_frame(stats.results_df, group, level, reference))
    return frames


def _covariate_table(
    obs: pd.DataFrame, sample_cols: list[str], covariates: Sequence[str]
) -> pd.DataFrame:
    covariates = list(covariates)
    grouped = obs[[*dict.fromkeys([*sample_cols, *covariates])]].groupby(
        sample_cols, observed=True
    )
    table = grouped.first()
    n_unique = grouped.nunique()
    if (varying := (n_unique > 1).any()).any():
        msg = (
            f"Covariates {varying[varying].index.tolist()} vary within samples. "
            "Covariates must be constant within each sample."
        )
        raise ValueError(msg)
    for col in sample_cols:
        if col in covariates:
            table[col] = table.index.get_level_values(col)
    return table[covariates]


def _design_metadata(
    pb_obs: pd.DataFrame, condition: str, covariate_table: pd.DataFrame
) -> pd.DataFrame:
    metadata = pd.DataFrame(
        {"condition": pd.Categorical(pb_obs[condition].astype(str))},
        index=pd.RangeIndex(len(pb_obs)),
    )
    samples = pd.MultiIndex.from_frame(pb_obs[covariate_table.index.names])
    for i, (_, values) in enumerate(covariate_table.reindex(samples).items()):
        col = values.reset_index(drop=True)
        if not pd.api.types.is_numeric_dtype(col):
            col = col.astype(str).astype("category")
        metadata[f"covariate_{i}"] = col
    return metadata


def _result_frame(
    results: pd.DataFrame, group: str | None, level: str, reference: str
) -> pd.DataFrame:
    return pd.DataFrame({
        "group": group,
        "condition": str(level),
        "reference": str(reference),
        "gene": results.index,
        "base_mean": results["baseMean"].to_numpy(),
        "log_fc": results["log2FoldChange"].to_numpy(),
        "log_fc_se": results["lfcSE"].to_numpy(),
        "statistic": results["stat"].to_numpy(),
        "p_value": results["pvalue"].to_numpy(),
        "adj_p_value": results["padj"].to_numpy(),
    })
