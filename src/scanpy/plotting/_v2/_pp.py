from __future__ import annotations

from typing import TYPE_CHECKING

import holoviews as hv
from hv_anndata import A

from ...get._out import Prefix, _out_dests, _out_key
from ...preprocessing._scrublet import SCRUBLET_BASES, SCRUBLET_KINDS
from .._common import highest_expr_genes as _highest_expr_genes

if TYPE_CHECKING:
    from collections.abc import Mapping

    import pandas as pd
    from anndata import AnnData
    from anndata.acc import AdRef

    from ...get._types import MultiAcc, RepAcc, VecRef


def highest_expr_genes(
    adata: AnnData,
    n_top: int = 20,
    *,
    use: RepAcc | str | None = None,
    gene_symbols: AdRef | None = None,
) -> hv.BoxWhisker:
    """Get ``n_top`` genes by mean expression.

    Parameters
    ----------
    adata
        The AnnData object.
    n_top
        The number of top genes to plot.
    use
        The count matrix to use, e.g. `A.layers['counts']`.
    gene_symbols
        Gene symbols to label the genes with, e.g. `A.var['symbol']`.
        Uses `.var_names` by default.

    Returns
    -------
    A box-and-whisker plot.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.highest_expr_genes(adata, use=A.layers["counts"])

    """
    hxg = _highest_expr_genes(adata, n_top, use=use, gene_symbols=gene_symbols)
    hxg_melted = hxg.melt(var_name="gene", value_name="frac_pct")
    return hv.BoxWhisker(hxg_melted, ["gene"], ["frac_pct"]).opts(
        ylabel="% of total counts", invert_axes=True
    )


def _hvg_table(
    adata: AnnData, use: MultiAcc | str | Prefix | None
) -> tuple[pd.DataFrame, dict]:
    """Find what `pp.highly_variable_genes` wrote, given its `out` (`None` for the default)."""
    if use is None or isinstance(use, Prefix):
        prefix = "" if use is None else f"{use.prefix}_"
        cols = [c for c in adata.var.columns if c.startswith(prefix)]
        table = adata.var[cols].rename(columns=lambda c: c.removeprefix(prefix))
        return table, adata.uns[f"{prefix}hvg"]
    key = _out_key(use, None, kind="varm", base="hvg")
    return adata.varm[key], adata.uns[key]


def highly_variable_genes(
    adata: AnnData, *, use: MultiAcc | str | Prefix | None = None
) -> hv.Layout:
    """Plot dispersions used to identify highly variable genes.

    Parameters
    ----------
    adata
        The AnnData object.
    use
        Where :func:`~scanpy.pp.highly_variable_genes` wrote its results, i.e. the `out` passed to it:
        `None` for the default places, `sc.Prefix('p')` for the `p_`-prefixed columns,
        or an accessor like `A.varm['hvg']`.

    Returns
    -------
    A layout containing two :class:`~holoviews.Scatter` plots,
    one normalized and one not.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        sc.pp.highly_variable_genes(adata, layer="counts")  # TODO: this should be the default
        sc.pl.highly_variable_genes(adata)

    """
    table, info = _hvg_table(adata, use)
    d1, d2 = (
        ("variances", "variances_norm")
        if info["flavor"].startswith("seurat_v3")
        else ("dispersions", "dispersions_norm")
    )

    return hv.Layout([
        hv.Scatter(table, ["means"], [d, "highly_variable"]).opts(
            color="highly_variable",
            **(  # https://github.com/holoviz/holoviews/issues/6938
                dict(cmap={True: "black", False: "gray"})
                if hv.Store.current_backend != "plotly"
                else {}
            ),
            legend_labels={
                True: "highly variable",
                False: "not highly variable",
            },
            legend_position="bottom_right",
            xlabel="mean expression of genes",
            ylabel=f"{d1} of genes ({'' if 'norm' in d else 'not '}normalized)",
        )
        for d in [d2, d1]
    ])


def scrublet_score_distribution(
    adata: AnnData, *, use: Mapping[str, VecRef | str] | Prefix | None = None
) -> hv.Layout:
    """Plot the doublet score distribution.

    Plots the doublet score probability densities for observed transcriptomes
    and simulated doublets.

    Parameters
    ----------
    adata
        The AnnData object.
    use
        Where :func:`~scanpy.pp.scrublet` wrote its results, i.e. the `out` passed to it:
        `None` for the default places, `sc.Prefix('p')` for the `p_`-prefixed ones,
        or a mapping like `{'doublet_score': A.obs['score']}`.

    Returns
    -------
    Layout containing two histograms.

    Examples
    --------

    ..  holoviews::

        import scanpy as sc

        sc.settings.preset = sc.Preset.ScanpyV2Preview
        A = sc.pl.hv_init(FAKE_BACKEND)

        adata = sc.datasets.pbmc68k_reduced()
        adata_sim = sc.pp.scrublet_simulate_doublets(adata)
        sc.pp.scrublet(adata, adata_sim)
        sc.pl.scrublet_score_distribution(adata)

    """
    uns_key, dests = _out_dests(
        False if use is None else use,
        None,
        dests=lambda _: SCRUBLET_BASES,
        kinds=SCRUBLET_KINDS,
        bases=SCRUBLET_BASES,
        uns="scrublet",
        default="scrublet",
    )
    score = A.obs[dests.doublet_score["obs"]]
    labels = dict(
        xlabel="Doublet score",
        ylabel="Probability density",
    )

    observed = (
        hv
        .Dataset(adata, [], [score])
        .hist(score, adjoin=False)
        .opts(
            xlim=(0, 1),
            logy=True,
            ylim=(1, None),
            title="Observed transcriptomes",
            **labels,
        )
    )

    doublets_sim = (
        hv
        .Table(adata.uns[uns_key]["doublet_scores_sim"], "scores")
        .hist("scores", adjoin=False)
        .opts(xlim=(0, 1), shared_axes=False, title="Simulated doublets", **labels)
    )

    return hv.Layout([observed, doublets_sim]) * hv.VLine(
        adata.uns[uns_key]["threshold"]
    )
