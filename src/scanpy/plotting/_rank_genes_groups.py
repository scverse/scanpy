from __future__ import annotations

from typing import Literal, Sequence

import numpy as np
import pandas as pd

from anndata import AnnData

from ._anndata import heatmap, matrixplot, dotplot


ArrayLike = np.ndarray
ValuesToPlot = Literal["scores", "logfoldchanges", "pvals", "pvals_adj"]


def _extract_rgg_values(
    adata: AnnData,
    values_to_plot: ValuesToPlot,
    groups: Sequence[str] | str | None,
    n_genes: int,
):
    """Extract dataframe (groups × genes) for the selected rank_genes_groups metric."""
    rgg = adata.uns.get("rank_genes_groups", None)
    if rgg is None:
        raise ValueError("`adata.uns['rank_genes_groups']` not found.")

    groups_order = rgg["names"].dtype.names
    if isinstance(groups, str) and groups != "all":
        groups = [groups]
    elif groups is None or groups == "all":
        groups = list(groups_order)

    # gather top N genes per group
    selected_genes = []
    for g in groups:
        arr = rgg["names"][g][:n_genes]
        selected_genes.extend(arr)

    selected_genes = list(dict.fromkeys(selected_genes))  # deduplicate, preserve order

    # build dataframe values_df[group][gene]
    df = pd.DataFrame(index=groups, columns=selected_genes, dtype=float)

    for g in groups:
        metrics = rgg[values_to_plot][g]  # ndarray length = total ranked genes
        names = rgg["names"][g]

        # map each gene to its metric
        mapping = {gene: metrics[i] for i, gene in enumerate(names)}

        # fill row
        df.loc[g, :] = [mapping.get(gn, np.nan) for gn in selected_genes]

    return df, selected_genes, groups


# ------------------------------------------------------------------------------
# Matrixplot
# ------------------------------------------------------------------------------

def rank_genes_groups_matrixplot(
    adata: AnnData,
    *,
    values_to_plot: ValuesToPlot | None = None,
    groups: Sequence[str] | str | None = None,
    n_genes: int = 20,
    **kwargs,
    MatrixPlot wrapper for rank_genes_groups with DE metric selection.

    Example:
        sc.pl.rank_genes_groups_matrixplot(
            adata,
            values_to_plot="logfoldchanges",
            groups=["0","1"],
            n_genes=20,
        )
    """
    if values_to_plot is None:
        # default: plot expression of marker genes

        raise ValueError(
            "`values_to_plot` must be provided. Options: "
            "['scores', 'logfoldchanges', 'pvals', 'pvals_adj']"
        )

    values_df, genes, groups = _extract_rgg_values(
        adata, values_to_plot, groups, n_genes
    )

    return matrixplot(
        adata,
        var_names=genes,
        groupby=groups,
        values_df=values_df,
        **kwargs,
    )


# ------------------------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------------------------

def rank_genes_groups_heatmap(
    adata: AnnData,
    *,
    values_to_plot: ValuesToPlot | None = None,
    groups: Sequence[str] | str | None = None,
    n_genes: int = 20,
    **kwargs,
):
    """
    Heatmap wrapper for rank_genes_groups with DE metric selection.
    """
    if values_to_plot is None:
        raise ValueError(
            "`values_to_plot` must be provided. Options: "
            "['scores', 'logfoldchanges', 'pvals', 'pvals_adj']"
        )

    values_df, genes, groups = _extract_rgg_values(
        adata, values_to_plot, groups, n_genes
    )

    return heatmap(
        adata,
        var_names=genes,
        groupby=groups,
        values_df=values_df,
        **kwargs,
    )


# ------------------------------------------------------------------------------
# Dotplot (for completeness parity with the issue text)
# ------------------------------------------------------------------------------

def rank_genes_groups_dotplot(
    adata: AnnData,
    *,
    values_to_plot: ValuesToPlot | None = None,
    groups: Sequence[str] | str | None = None,
    n_genes: int = 20,
    **kwargs,
):
    """
    DotPlot wrapper for rank_genes_groups with DE metric selection.
    This adds parity with the existing sc.pl.rank_genes_groups_dotplot API.
    """
    if values_to_plot is None:
        raise ValueError(
            "`values_to_plot` must be provided. Options: "
            "['scores', 'logfoldchanges', 'pvals', 'pvals_adj']"
        )

    values_df, genes, groups = _extract_rgg_values(
        adata, values_to_plot, groups, n_genes
    )

    # DotPlot uses values_df as dot_color_df
    return dotplot(
        adata,
        var_names=genes,
        groupby=groups,
        dot_color_df=values_df,
        **kwargs,
    )

