from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .._compat import old_positionals
from ._utils import savefig_or_show

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    Scale = Literal["linear", "log", "symlog", "logit"] | str


@old_positionals(
    "scale_hist_obs", "scale_hist_sim", "figsize", "return_fig", "show", "save"
)
def scrublet_score_distribution(
    adata: AnnData,
    *,
    scale_hist_obs: Scale = "log",
    scale_hist_sim: Scale = "linear",
    figsize: tuple[float | int, float | int] = (8, 3),
    return_fig: bool = False,
    show: bool = True,
    save: str | bool | None = None,
) -> Figure | Sequence[tuple[Axes, Axes]] | tuple[Axes, Axes] | None:
    """Plot histogram of doublet scores for observed transcriptomes and simulated doublets.

    The histogram for simulated doublets is useful for determining the correct doublet
    score threshold.

    Scrublet must have been run previously with the input object.

    Parameters
    ----------
    adata
        An AnnData object resulting from :func:`~scanpy.pp.scrublet`.
    scale_hist_obs
        Set y axis scale transformation in matplotlib for the plot of observed transcriptomes
    scale_hist_sim
        Set y axis scale transformation in matplotlib for the plot of simulated doublets
    figsize
        width, height
    show
        Show the plot, do not return axis.
    save
        If :data:`True` or a :class:`str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.

    Returns
    -------
    If ``return_fig`` is True, a :class:`~matplotlib.figure.Figure`.
    If ``show==False`` a list of :class:`~matplotlib.axes.Axes`.

    See Also
    --------
    :func:`~scanpy.pp.scrublet`: Main way of running Scrublet, runs
        preprocessing, doublet simulation and calling.
    :func:`~scanpy.pp.scrublet_simulate_doublets`: Run Scrublet's doublet
        simulation separately for advanced usage.

    """
    if "scrublet" not in adata.uns:
        msg = "Please run scrublet before trying to generate the scrublet plot."
        raise ValueError(msg)

    # If batched_by is populated, then we know Scrublet was run over multiple batches

    if "batched_by" in adata.uns["scrublet"]:
        batched_by = adata.uns["scrublet"]["batched_by"]
        batches = adata.obs[batched_by].astype("category", copy=False)
        n_batches = len(batches.cat.categories)
        figsize = (figsize[0], figsize[1] * n_batches)
    else:
        batches = pd.Series(
            np.broadcast_to(0, adata.n_obs), dtype="category", index=adata.obs_names
        )
        n_batches = 1

    fig, axs = plt.subplots(n_batches, 2, figsize=figsize)

    for idx, (batch_key, sub_obs) in enumerate(
        adata.obs.groupby(batches, observed=True)
    ):
        obs_ax: Axes
        sim_ax: Axes
        # We'll need multiple rows if Scrublet was run in multiple batches
        if "batched_by" in adata.uns["scrublet"]:
            threshold = adata.uns["scrublet"]["batches"][batch_key].get(
                "threshold", None
            )
            doublet_scores_sim = adata.uns["scrublet"]["batches"][batch_key][
                "doublet_scores_sim"
            ]
            axis_lab_suffix = f" ({batch_key})"
            obs_ax = axs[idx][0]
            sim_ax = axs[idx][1]

        else:
            threshold = adata.uns["scrublet"].get("threshold", None)
            doublet_scores_sim = adata.uns["scrublet"]["doublet_scores_sim"]
            axis_lab_suffix = ""
            obs_ax = axs[0]
            sim_ax = axs[1]

        # Make the plots
        _plot_scores(
            obs_ax,
            sub_obs["doublet_score"],
            scale=scale_hist_obs,
            title=f"Observed transcriptomes {axis_lab_suffix}",
            threshold=threshold,
        )
        _plot_scores(
            sim_ax,
            doublet_scores_sim,
            scale=scale_hist_sim,
            title=f"Simulated doublets {axis_lab_suffix}",
            threshold=threshold,
        )

    fig.tight_layout()

    savefig_or_show("scrublet_score_distribution", show=show, save=save)
    if return_fig:
        return fig
    elif not show:
        return axs


def _plot_scores(
    ax: Axes,
    scores: np.ndarray,
    scale: Scale,
    title: str,
    threshold: float | None = None,
) -> None:
    ax.hist(
        scores,
        np.linspace(0, 1, 50),
        color="gray",
        linewidth=0,
        density=True,
    )
    ax.set_yscale(scale)
    yl = ax.get_ylim()
    ax.set_ylim(yl)

    if threshold is not None:
        ax.plot(threshold * np.ones(2), yl, c="black", linewidth=1)

    ax.set_title(title)
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Prob. density")
