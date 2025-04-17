"""Calculate scores based on the expression of gene lists."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .. import logging as logg
from .._compat import CSBase, old_positionals
from .._utils import _check_use_raw, is_backed_type
from ..get import _get_obs_rep

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Sequence
    from typing import Literal

    from anndata import AnnData
    from numpy.typing import DTypeLike, NDArray

    from .._utils.random import _LegacyRandom

    try:
        _StrIdx = pd.Index[str]
    except TypeError:  # Sphinx
        _StrIdx = pd.Index
    _GetSubset = Callable[[_StrIdx], np.ndarray | CSBase]


def _sparse_nanmean(X: CSBase, axis: Literal[0, 1]) -> NDArray[np.float64]:
    """np.nanmean equivalent for sparse matrices."""
    if not isinstance(X, CSBase):
        msg = "X must be a compressed sparse matrix"
        raise TypeError(msg)

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype="float64")  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m


@old_positionals(
    "ctrl_size", "gene_pool", "n_bins", "score_name", "random_state", "copy", "use_raw"
)
def score_genes(  # noqa: PLR0913
    adata: AnnData,
    gene_list: Sequence[str] | pd.Index[str],
    *,
    ctrl_as_ref: bool = True,
    ctrl_size: int = 50,
    gene_pool: Sequence[str] | pd.Index[str] | None = None,
    n_bins: int = 25,
    score_name: str = "score",
    random_state: _LegacyRandom = 0,
    copy: bool = False,
    use_raw: bool | None = None,
    layer: str | None = None,
) -> AnnData | None:
    """Score a set of genes :cite:p:`Satija2015`.

    The score is the average expression of a set of genes after subtraction by
    the average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat :cite:p:`Satija2015` and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    ctrl_as_ref
        Allow the algorithm to use the control genes as reference.
        Will be changed to `False` in scanpy 2.0.
    ctrl_size
        Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in `.obs`.
    random_state
        The random seed for sampling.
    copy
        Copy `adata` or modify it inplace.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.

        .. versionchanged:: 1.4.5
           Default value changed from `False` to `None`.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following field:

    `adata.obs[score_name]` : :class:`numpy.ndarray` (dtype `float`)
        Scores of each cell.

    Examples
    --------
    See this `notebook <https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle>`__.

    """
    start = logg.info(f"computing score {score_name!r}")
    adata = adata.copy() if copy else adata
    use_raw = _check_use_raw(adata, use_raw, layer=layer)
    if is_backed_type(adata.X) and not use_raw:
        msg = f"score_genes is not implemented for matrices of type {type(adata.X)}"
        raise NotImplementedError(msg)

    if random_state is not None:
        np.random.seed(random_state)

    gene_list, gene_pool, get_subset = _check_score_genes_args(
        adata, gene_list, gene_pool, use_raw=use_raw, layer=layer
    )
    del use_raw, layer, random_state

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    control_genes = pd.Index([], dtype="string")
    for r_genes in _score_genes_bins(
        gene_list,
        gene_pool,
        ctrl_as_ref=ctrl_as_ref,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        get_subset=get_subset,
    ):
        control_genes = control_genes.union(r_genes)

    if len(control_genes) == 0:
        msg = "No control genes found in any cut."
        if ctrl_as_ref:
            msg += " Try setting `ctrl_as_ref=False`."
        raise RuntimeError(msg)

    means_list, means_control = (
        _nan_means(get_subset(genes), axis=1, dtype="float64")
        for genes in (gene_list, control_genes)
    )
    score = means_list - means_control

    adata.obs[score_name] = pd.Series(
        np.array(score).ravel(), index=adata.obs_names, dtype="float64"
    )

    logg.info(
        "    finished",
        time=start,
        deep=(
            "added\n"
            f"    {score_name!r}, score of gene set (adata.obs).\n"
            f"    {len(control_genes)} total control genes are used."
        ),
    )
    return adata if copy else None


def _check_score_genes_args(
    adata: AnnData,
    gene_list: pd.Index[str] | Sequence[str],
    gene_pool: pd.Index[str] | Sequence[str] | None,
    *,
    layer: str | None,
    use_raw: bool,
) -> tuple[pd.Index[str], pd.Index[str], _GetSubset]:
    """Restrict `gene_list` and `gene_pool` to present genes in `adata`.

    Also returns a function to get subset of `adata.X` based on a set of genes passed.
    """
    var_names = adata.raw.var_names if use_raw else adata.var_names
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    genes_to_ignore = gene_list.difference(var_names, sort=False)  # first get missing
    gene_list = gene_list.intersection(var_names)  # then restrict to present
    if len(genes_to_ignore) > 0:
        logg.warning(f"genes are not in var_names and ignored: {genes_to_ignore}")
    if len(gene_list) == 0:
        msg = "No valid genes were passed for scoring."
        raise ValueError(msg)

    if gene_pool is None:
        gene_pool = var_names.astype("string")
    else:
        gene_pool = pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        msg = "No valid genes were passed for reference set."
        raise ValueError(msg)

    def get_subset(genes: pd.Index[str]):
        x = _get_obs_rep(adata, use_raw=use_raw, layer=layer)
        if len(genes) == len(var_names):
            return x
        idx = var_names.get_indexer(genes)
        return x[:, idx]

    return gene_list, gene_pool, get_subset


def _score_genes_bins(
    gene_list: pd.Index[str],
    gene_pool: pd.Index[str],
    *,
    ctrl_as_ref: bool,
    ctrl_size: int,
    n_bins: int,
    get_subset: _GetSubset,
) -> Generator[pd.Index[str], None, None]:
    # average expression of genes
    obs_avg = pd.Series(_nan_means(get_subset(gene_pool), axis=0), index=gene_pool)
    # Sometimes (and I don’t know how) missing data may be there, with NaNs for missing entries
    obs_avg = obs_avg[np.isfinite(obs_avg)]

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    keep_ctrl_in_obs_cut = False if ctrl_as_ref else obs_cut.index.isin(gene_list)

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes: pd.Index[str] = obs_cut[(obs_cut == cut) & ~keep_ctrl_in_obs_cut].index
        if len(r_genes) == 0:
            msg = (
                f"No control genes for {cut=}. You might want to increase "
                f"gene_pool size (current size: {len(gene_pool)})"
            )
            logg.warning(msg)
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size).index
        if ctrl_as_ref:  # otherwise `r_genes` is already filtered
            r_genes = r_genes.difference(gene_list)
        yield r_genes


def _nan_means(
    x: np.ndarray | CSBase, *, axis: Literal[0, 1], dtype: DTypeLike | None = None
) -> NDArray[np.float64]:
    if isinstance(x, CSBase):
        return np.array(_sparse_nanmean(x, axis=axis)).flatten()
    return np.nanmean(x, axis=axis, dtype=dtype)


@old_positionals("s_genes", "g2m_genes", "copy")
def score_genes_cell_cycle(
    adata: AnnData,
    *,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    copy: bool = False,
    **kwargs,
) -> AnnData | None:
    """Score cell cycle genes :cite:p:`Satija2015`.

    Given two lists of genes associated to S phase and G2M phase, calculates
    scores and assigns a cell cycle phase (G1, S or G2M). See
    :func:`~scanpy.tl.score_genes` for more explanation.

    Parameters
    ----------
    adata
        The annotated data matrix.
    s_genes
        List of genes associated with S phase.
    g2m_genes
        List of genes associated with G2M phase.
    copy
        Copy `adata` or modify it inplace.
    **kwargs
        Are passed to :func:`~scanpy.tl.score_genes`. `ctrl_size` is not
        possible, as it's set as `min(len(s_genes), len(g2m_genes))`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obs['S_score']` : :class:`pandas.Series` (dtype `object`)
        The score for S phase for each cell.
    `adata.obs['G2M_score']` : :class:`pandas.Series` (dtype `object`)
        The score for G2M phase for each cell.
    `adata.obs['phase']` : :class:`pandas.Series` (dtype `object`)
        The cell cycle phase (`S`, `G2M` or `G1`) for each cell.

    See Also
    --------
    score_genes

    Examples
    --------
    See this `notebook <https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle>`__.

    """
    logg.info("calculating cell cycle phase")

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    for genes, name in [(s_genes, "S_score"), (g2m_genes, "G2M_score")]:
        score_genes(adata, genes, score_name=name, ctrl_size=ctrl_size, **kwargs)
    scores = adata.obs[["S_score", "G2M_score"]]

    # default phase is S
    phase = pd.Series("S", index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores["G2M_score"] > scores["S_score"]] = "G2M"

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = "G1"

    adata.obs["phase"] = phase
    logg.hint("    'phase', cell cycle phase (adata.obs)")
    return adata if copy else None
