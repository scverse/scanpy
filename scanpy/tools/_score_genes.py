"""Calculate scores based on the expression of gene lists."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from scanpy._utils import _check_use_raw

from .. import logging as logg
from .._compat import old_positionals

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

    from .._utils import AnyRandom


def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

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
def score_genes(
    adata: AnnData,
    gene_list: Sequence[str] | pd.Index[str],
    *,
    ctrl_size: int = 50,
    gene_pool: Sequence[str] | pd.Index[str] | None = None,
    n_bins: int = 25,
    score_name: str = "score",
    random_state: AnyRandom = 0,
    copy: bool = False,
    use_raw: bool | None = None,
) -> AnnData | None:
    """\
    Score a set of genes [Satija15]_.

    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
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
    use_raw = _check_use_raw(adata, use_raw)

    if random_state is not None:
        np.random.seed(random_state)

    var_names = adata.raw.var_names if use_raw else adata.var_names
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    genes_to_ignore = gene_list.difference(var_names, sort=False)  # first get missing
    gene_list = gene_list.intersection(var_names)  # then restrict to present
    if len(genes_to_ignore) > 0:
        logg.warning(f"genes are not in var_names and ignored: {genes_to_ignore}")
    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")

    if gene_pool is None:
        gene_pool = pd.Index(var_names, dtype="string")
    else:
        gene_pool = pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    _adata = adata.raw if use_raw else adata
    _adata_subset = (
        _adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
    )
    # average expression of genes
    if issparse(_adata_subset.X):
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(),
            index=gene_pool,
        )
    else:
        obs_avg = pd.Series(np.nanmean(_adata_subset.X, axis=0), index=gene_pool)

    # Sometimes (and I don't know how) missing data may be there, with nansfor
    obs_avg = obs_avg[np.isfinite(obs_avg)]

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = pd.Index([], dtype="string")

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes: pd.Index[str] = obs_cut[obs_cut == cut].index
        if ctrl_size < len(r_genes):
            r_genes = r_genes.to_series().sample(ctrl_size).index
        control_genes = control_genes.union(r_genes.difference(gene_list))

    X_list = _adata[:, gene_list].X
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
    else:
        X_list = np.nanmean(X_list, axis=1, dtype="float64")

    X_control = _adata[:, control_genes].X
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1, dtype="float64")

    score = X_list - X_control

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


@old_positionals("s_genes", "g2m_genes", "copy")
def score_genes_cell_cycle(
    adata: AnnData,
    *,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    copy: bool = False,
    **kwargs,
) -> AnnData | None:
    """\
    Score cell cycle genes [Satija15]_.

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

    See also
    --------
    score_genes

    Examples
    --------
    See this `notebook <https://github.com/scverse/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    logg.info("calculating cell cycle phase")

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    score_genes(
        adata, gene_list=s_genes, score_name="S_score", ctrl_size=ctrl_size, **kwargs
    )
    # add g2m-score
    score_genes(
        adata,
        gene_list=g2m_genes,
        score_name="G2M_score",
        ctrl_size=ctrl_size,
        **kwargs,
    )
    scores = adata.obs[["S_score", "G2M_score"]]

    # default phase is S
    phase = pd.Series("S", index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores.G2M_score > scores.S_score] = "G2M"

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = "G1"

    adata.obs["phase"] = phase
    logg.hint("    'phase', cell cycle phase (adata.obs)")
    return adata if copy else None
