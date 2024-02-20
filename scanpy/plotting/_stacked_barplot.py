"""Stacked barplots for analyzing compositions of cell groups.
"""
from typing import TYPE_CHECKING, Literal
from collections.abc import Sequence
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Colormap, ListedColormap, is_color_like
from matplotlib.axes import Axes
import seaborn as sns
import numpy as np
from pandas.api.types import CategoricalDtype
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from anndata import AnnData
from cycler import Cycler

from .._settings import settings
from ._utils import ColorLike
from . import _utils


def stacked_barplot(
    adata: AnnData, 
    *,
    obs_key: str,
    groupby: str,
    order: Sequence[str] | None = None,
    order_method: str | None = None,
    orderby: Literal["value", "dendrogram"] | None = None,
    plot_dendrogram: bool = False,
    palette: Cycler | ListedColormap | ColorLike | Sequence[ColorLike] | None = None,
    rotation: float | None = None,
    figsize: tuple[float, float] | None = None, 
    show: bool | None = None,
    save: bool | str | None = None,
    ax: Axes | None = None
) -> Axes | tuple[Axes, Axes] | None:
    """
    Stacked barplot.

    Makes a stacked barplot for showing compositions of groups of cells (e.g., fraction of each cluster 
    originating from each batch).

    Parameters
    ----------
    adata
        Annotated data matrix.
    obs_key 
        Key for accessing fields of `.obs`.
    groupby
        The key of the observation grouping to consider.
    order
        Order in which to show the groups.
    order_method
        If set to `value`, then order the groups by either a specific 
        value in `adata.obs[obs_key]` as specified by the `orderby`
        argument. If set to `dendrogram`, then order by agglomerative 
        clustering. Will only work if `order` is not specified.
    orderby
        Order the groups according to this value in `adata.obs[obs_key]`.
        Will only work if `order` is not specified and `order_method="value"`.
    plot_dendrogram
        Plot the dendrogram. Will only work if `order_method="dendrogram"`.
    palette
        Colors to use for different values of `adata.obs[obs_key]`.
    rotation
        Rotation of xtick labels.
    figsize
        Figure size. Otherwise the rcParam[‘figure.figsize] value is used. Format is (width, height).
    show
        Show the plot, do not return axis.
    save
        If True or a str, save the figure. A string is appended to the default filename. 
        Infer the filetype if ending on {‘.pdf’, ‘.png’, ‘.svg’}.
    ax
        A matplotlib axes object. Only works if `order_method` is not "dendrogram".

    Returns
    -------
    If ``show==False`` and ``plot_dendrogram!=True``, return a :class:`~matplotlib.axes.Axes`.
    If ``show==False``, ``order_method=="dendrogram``, and ``plot_dendrogram==True``, return a list of :class:`~matplotlib.axes.Axes`.
    Else return ``None``.
    """
    if not isinstance(adata.obs[groupby].dtype, CategoricalDtype) \
        or not not isinstance(adata.obs[groupby].dtype, str):
        raise ValueError(
            f"The column `adata.obs[{groupby!r}]` needs to be categorical, "
            f"but is of dtype {adata.obs[groupby].dtype}."
        )
    if not isinstance(adata.obs[obs_key].dtype, CategoricalDtype) \
        or not not isinstance(adata.obs[obs_key].dtype, str):
        raise ValueError(
            f"The column `adata.obs[{obs_key!r}]` needs to be categorical, "
            f"but is of dtype {adata.obs[obs_key].dtype}."
        )

    palette = _utils.default_palette(palette)

    if figsize is None:
        figsize = rcParams["figure.figsize"]
    if order_method == 'dendrogram' and plot_dendrogram:
        fig, axarr = plt.subplots(
            2, 1, figsize=figsize, gridspec_kw={'height_ratios': [1, 4]}
        )
        fig.subplots_adjust(hspace=0)
        ax_dendro = axarr[0]
        ax = axarr[1]
    else:
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Get all groups and categories. Map categories to colors and
    # groups to the number of cells assigned to those groups
    groups = sorted(set(adata.obs[groupby]))
    cats = sorted(set(adata.obs[obs_key]))
    cat_to_color = {
        cat: color['color']
        for cat, color in zip(cats, palette())
    }
    group_to_size = {
        group: adata.obs.loc[adata.obs[groupby] == group].shape[0]
        for group in groups
    }

    if order is not None:
        assert frozenset(order) == frozenset(groups), f"`order` must include all unique values in {adata.obs[groupby]}"
        groups = order
    else:
        if order_method == 'value':
            if orderby is None:
                raise ValueError(
                    f"If `order_method` is set to `value`, then a category in "
                    f"{adata.obs[obs_key]} must be supplied to `orderby`."
                )
            group_to_val = {
                group: adata.obs.loc[
                    (adata.obs[obs_key] == orderby)
                    & (adata.obs[groupby] == group)
                ].shape[0] / group_to_size[group]
                for group in groups
            }
            groups = sorted(groups, key=lambda x: group_to_val[x])
        elif order_method == 'dendrogram':
            fracs = [[] for i in range(len(groups))]
            for i, group in enumerate(groups):
                for cat in cats:
                    frac = _compute_proportion(
                        adata, obs_key, groupby, group, cat, group_to_size
                    )
                    fracs[i].append(frac)
            Z = linkage(fracs, method='ward')
            leaf_order = leaves_list(Z)
            groups = list(np.array(groups)[leaf_order])

            if plot_dendrogram:
                dendro = dendrogram(
                    Z, ax=ax_dendro, orientation='top', color_threshold=0, 
                    link_color_func=lambda k: 'black'
                )
                ax_dendro.axis('off')
                ax_dendro.set_xticks([])
                ax_dendro.set_yticks([])

    # Create barplot 
    bottom = np.zeros(len(groups)) # Track current heights
    for cat in cats:
        fracs = []
        for group in groups:
            fracs.append(
                _compute_proportion(adata, obs_key, groupby, group, cat, group_to_size)
            )
        bars = ax.bar(groups, fracs, bottom=bottom, color=cat_to_color[cat])
        bottom += fracs

    # Setup legend
    patches = [
        mpatches.Patch(color=cat_to_color[cat], label=cat) for cat in cats
    ]
    ax.legend(
        handles=patches, bbox_to_anchor=(1.02, 1.), loc='upper left', borderaxespad=0,
        frameon=False
    )

    if rotation is not None:
        ax.tick_params(axis="x", labelrotation=rotation)
    ax.margins(0.0, 0)
    ax.set_xlabel(groupby)
    ax.set_ylabel(obs_key)

    show = settings.autoshow if show is None else show
    _utils.savefig_or_show("stacked_barplot", show=show, save=save)
    if show:
        return None
    else:
        if order_method == "dendrogram" and plot_dendrogram:
            return ax, ax_dendro 
        else:
            return ax


def _compute_proportion(adata, obs_key, groupby, group, cat, group_to_size):
    n_cells = adata.obs.loc[
        (adata.obs[obs_key] == cat)
        & (adata.obs[groupby] == group)
    ].shape[0]
    frac = n_cells / group_to_size[group]
    return frac

