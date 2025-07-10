"""Plotting functions for AnnData."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Collection, Mapping, Sequence
from itertools import pairwise, product
from types import NoneType
from typing import TYPE_CHECKING, NamedTuple, TypedDict, cast

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import colormaps, gridspec, patheffects, rcParams
from matplotlib import pyplot as plt
from matplotlib.colors import is_color_like
from packaging.version import Version
from pandas.api.types import CategoricalDtype, is_numeric_dtype

from .. import get
from .. import logging as logg
from .._compat import CSBase, old_positionals
from .._settings import settings
from .._utils import (
    _check_use_raw,
    _doc_params,
    _empty,
    get_literal_vals,
    sanitize_anndata,
)
from . import _utils
from ._docs import (
    doc_common_plot_args,
    doc_scatter_basic,
    doc_show_save_ax,
    doc_vboundnorm,
)
from ._utils import (
    _deprecated_scale,
    _dk,
    check_colornorm,
    scatter_base,
    scatter_group,
    setup_axes,
)

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal, Self

    from anndata import AnnData
    from cycler import Cycler
    from matplotlib.axes import Axes
    from matplotlib.colors import Colormap, ListedColormap, Normalize
    from numpy.typing import NDArray
    from seaborn import FacetGrid
    from seaborn.matrix import ClusterGrid

    from .._utils import Empty
    from ._utils import (
        ColorLike,
        DensityNorm,
        _FontSize,
        _FontWeight,
        _LegendLoc,
    )

    # TODO: is that all?
    _Basis = Literal["pca", "tsne", "umap", "diffmap", "draw_graph_fr"]
    _VarNames = str | Sequence[str]


class VarGroups(NamedTuple):
    labels: Sequence[str]
    """Var labels."""
    positions: Sequence[tuple[int, int]]
    """Var positions.

    Each item in the list should contain the start and end position that the bracket should cover.
    Eg. `[(0, 4), (5, 8)]` means that there are two brackets,
    one for the var_names (eg genes) in positions 0-4 and other for positions 5-8
    """

    @classmethod
    def validate(
        cls, labels: Sequence[str] | None, positions: Sequence[tuple[int, int]] | None
    ) -> Self | None:
        if labels is None and positions is None:
            return None
        if labels is None or positions is None:
            msg = (
                "If var_group_labels or var_group_positions are given, "
                "both have to be given."
            )
            raise ValueError(msg)
        if len(labels) != len(positions):
            msg = (
                "var_group_labels and var_group_positions must have the same length. "
                f"Got {len(labels)=} and {len(positions)=}."
            )
            raise ValueError(msg)
        return None if len(labels) == 0 else cls(labels, positions)


@old_positionals(
    "color",
    "use_raw",
    "layers",
    "sort_order",
    "alpha",
    "basis",
    "groups",
    "components",
    "projection",
    "legend_loc",
    "legend_fontsize",
    "legend_fontweight",
    "legend_fontoutline",
    "color_map",
    # 17 positionals are enough for backwards compatibility
)
@_doc_params(scatter_temp=doc_scatter_basic, show_save_ax=doc_show_save_ax)
def scatter(  # noqa: PLR0913
    adata: AnnData,
    x: str | None = None,
    y: str | None = None,
    *,
    color: str | ColorLike | Collection[str | ColorLike] | None = None,
    use_raw: bool | None = None,
    layers: str | Collection[str] | None = None,
    sort_order: bool = True,
    alpha: float | None = None,
    basis: _Basis | None = None,
    groups: str | Iterable[str] | None = None,
    components: str | Collection[str] | None = None,
    projection: Literal["2d", "3d"] = "2d",
    legend_loc: _LegendLoc | None = "right margin",
    legend_fontsize: float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight | None = None,
    legend_fontoutline: float | None = None,
    color_map: str | Colormap | None = None,
    palette: Cycler | ListedColormap | ColorLike | Sequence[ColorLike] | None = None,
    frameon: bool | None = None,
    right_margin: float | None = None,
    left_margin: float | None = None,
    size: float | None = None,
    marker: str | Sequence[str] = ".",
    title: str | Collection[str] | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: Axes | None = None,
) -> Axes | list[Axes] | None:
    """Scatter plot along observations or variables axes.

    Color the plot using annotations of observations (`.obs`), variables
    (`.var`) or expression of genes (`.var_names`).

    Parameters
    ----------
    adata
        Annotated data matrix.
    x
        x coordinate.
    y
        y coordinate.
    color
        Keys for annotations of observations/cells or variables/genes,
        or a hex color specification, e.g.,
        `'ann1'`, `'#fe57a1'`, or `['ann1', 'ann2']`.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    layers
        Use the `layers` attribute of `adata` if present: specify the layer for
        `x`, `y` and `color`. If `layers` is a string, then it is expanded to
        `(layers, layers, layers)`.
    basis
        String that denotes a plotting tool that computed coordinates.
    {scatter_temp}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    """
    # color can be a obs column name or a matplotlib color specification (or a collection thereof)
    if color is not None:
        color = cast(
            "Collection[str | ColorLike]",
            [color] if isinstance(color, str) or is_color_like(color) else color,
        )
    args = locals()

    if basis is not None:
        return _scatter_obs(**args)
    if x is None or y is None:
        msg = "Either provide a `basis` or `x` and `y`."
        raise ValueError(msg)
    if _check_if_annotations(adata, "obs", x=x, y=y, colors=color, use_raw=use_raw):
        return _scatter_obs(**args)
    if _check_if_annotations(adata, "var", x=x, y=y, colors=color, use_raw=use_raw):
        args_t = {**args, "adata": adata.T}
        axs = _scatter_obs(**args_t)
        # store .uns annotations that were added to the new adata object
        adata.uns = args_t["adata"].uns
        return axs
    msg = (
        "`x`, `y`, and potential `color` inputs must all "
        "come from either `.obs` or `.var`"
    )
    raise ValueError(msg)


def _check_if_annotations(
    adata: AnnData,
    axis_name: Literal["obs", "var"],
    *,
    x: str | None = None,
    y: str | None = None,
    colors: Collection[str | ColorLike] | None = None,
    use_raw: bool | None = None,
) -> bool:
    """Check if `x`, `y`, and `colors` are annotations of `adata`.

    In the case of `colors`, valid matplotlib colors are also accepted.

    If `axis_name` is `obs`, checks in `adata.obs.columns` and `adata.var_names`,
    if `axis_name` is `var`, checks in `adata.var.columns` and `adata.obs_names`.
    """
    annotations: pd.Index[str] = getattr(adata, axis_name).columns
    other_ax_obj = (
        adata.raw if _check_use_raw(adata, use_raw) and axis_name == "obs" else adata
    )
    names: pd.Index[str] = getattr(
        other_ax_obj, "var" if axis_name == "obs" else "obs"
    ).index

    def is_annotation(needle: pd.Index) -> NDArray[np.bool_]:
        return needle.isin({None}) | needle.isin(annotations) | needle.isin(names)

    if not is_annotation(pd.Index([x, y])).all():
        return False

    color_idx = pd.Index(colors if colors is not None else [])
    # Colors are valid
    color_valid: NDArray[np.bool_] = np.fromiter(
        map(is_color_like, color_idx), dtype=np.bool_, count=len(color_idx)
    )
    # Annotation names are valid too
    color_valid[~color_valid] = is_annotation(color_idx[~color_valid])
    return bool(color_valid.all())


def _scatter_obs(  # noqa: PLR0912, PLR0913, PLR0915
    *,
    adata: AnnData,
    x: str | None = None,
    y: str | None = None,
    color: Collection[str | ColorLike] | None = None,
    use_raw: bool | None = None,
    layers: str | Collection[str] | None = None,
    sort_order: bool = True,
    alpha: float | None = None,
    basis: _Basis | None = None,
    groups: str | Iterable[str] | None = None,
    components: str | Collection[str] | None = None,
    projection: Literal["2d", "3d"] = "2d",
    legend_loc: _LegendLoc | None = "right margin",
    legend_fontsize: float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight | None = None,
    legend_fontoutline: float | None = None,
    color_map: str | Colormap | None = None,
    palette: Cycler | ListedColormap | ColorLike | Sequence[ColorLike] | None = None,
    frameon: bool | None = None,
    right_margin: float | None = None,
    left_margin: float | None = None,
    size: float | None = None,
    marker: str | Sequence[str] = ".",
    title: str | Collection[str] | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: Axes | None = None,
) -> Axes | list[Axes] | None:
    """See docstring of scatter."""
    sanitize_anndata(adata)

    use_raw = _check_use_raw(adata, use_raw)

    # Process layers
    if layers in ["X", None] or (isinstance(layers, str) and layers in adata.layers):
        layers = (layers, layers, layers)
    elif isinstance(layers, Collection) and len(layers) == 3:
        layers = tuple(layers)
        for layer in layers:
            if layer not in adata.layers and layer not in ["X", None]:
                msg = (
                    "`layers` should have elements that are "
                    "either None or in adata.layers.keys()."
                )
                raise ValueError(msg)
    else:
        msg = (
            "`layers` should be a string or a collection of strings "
            f"with length 3, had value '{layers}'"
        )
        raise ValueError(msg)
    if use_raw and layers not in [("X", "X", "X"), (None, None, None)]:
        msg = "`use_raw` must be `False` if layers are used."
        raise ValueError(msg)

    if legend_loc not in (valid_legend_locs := get_literal_vals(_utils._LegendLoc)):
        msg = f"Invalid `legend_loc`, need to be one of: {valid_legend_locs}."
        raise ValueError(msg)
    if components is None:
        components = "1,2" if "2d" in projection else "1,2,3"
    if isinstance(components, str):
        components = components.split(",")
    components = np.array(components).astype(int) - 1
    keys = ["grey"] if color is None else color
    if title is not None and isinstance(title, str):
        title = [title]
    highlights = adata.uns.get("highlights", [])
    if basis is not None:
        try:
            # ignore the '0th' diffusion component
            if basis == "diffmap":
                components += 1
            Y = adata.obsm["X_" + basis][:, components]
            # correct the component vector for use in labeling etc.
            if basis == "diffmap":
                components -= 1
        except KeyError:
            msg = f"compute coordinates using visualization tool {basis} first"
            raise KeyError(msg) from None
    elif x is not None and y is not None:
        if use_raw:
            if x in adata.obs.columns:
                x_arr = adata.obs_vector(x)
            else:
                x_arr = adata.raw.obs_vector(x)
            if y in adata.obs.columns:
                y_arr = adata.obs_vector(y)
            else:
                y_arr = adata.raw.obs_vector(y)
        else:
            x_arr = adata.obs_vector(x, layer=layers[0])
            y_arr = adata.obs_vector(y, layer=layers[1])

        Y = np.c_[x_arr, y_arr]
    else:
        msg = "Either provide a `basis` or `x` and `y`."
        raise ValueError(msg)

    if size is None:
        n = Y.shape[0]
        size = 120000 / n

    if legend_fontsize is None:
        legend_fontsize = rcParams["legend.fontsize"]

    if isinstance(palette, Sequence) and not isinstance(palette, str):
        palettes = palette if not is_color_like(palette[0]) else [palette]
    else:
        palettes = [palette for _ in range(len(keys))]
    palettes = [_utils.default_palette(palette) for palette in palettes]

    if basis is not None:
        component_name = (
            "DC"
            if basis == "diffmap"
            else "tSNE"
            if basis == "tsne"
            else "UMAP"
            if basis == "umap"
            else "PC"
            if basis == "pca"
            else "TriMap"
            if basis == "trimap"
            else basis.replace("draw_graph_", "").upper()
            if "draw_graph" in basis
            else basis
        )
    else:
        component_name = None
    axis_labels = (x, y) if component_name is None else None
    show_ticks = component_name is None

    # generate the colors
    color_ids: list[np.ndarray | ColorLike] = []
    categoricals = []
    colorbars = []
    for ikey, key in enumerate(keys):
        c = "white"
        categorical = False  # by default, assume continuous or flat color
        colorbar = None
        # test whether we have categorial or continuous annotation
        if key in adata.obs_keys():
            if isinstance(adata.obs[key].dtype, CategoricalDtype):
                categorical = True
            else:
                c = adata.obs[key].to_numpy()
        # coloring according to gene expression
        elif use_raw and adata.raw is not None and key in adata.raw.var_names:
            c = adata.raw.obs_vector(key)
        elif key in adata.var_names:
            c = adata.obs_vector(key, layer=layers[2])
        elif is_color_like(key):  # a flat color
            c = key
            colorbar = False
        else:
            msg = (
                f"key {key!r} is invalid! pass valid observation annotation, "
                f"one of {adata.obs_keys()} or a gene name {adata.var_names}"
            )
            raise ValueError(msg)
        if colorbar is None:
            colorbar = not categorical
        colorbars.append(colorbar)
        if categorical:
            categoricals.append(ikey)
        color_ids.append(c)

    if right_margin is None and len(categoricals) > 0 and legend_loc == "right margin":
        right_margin = 0.5
    if title is None and keys[0] is not None:
        title = [
            key.replace("_", " ") if not is_color_like(key) else "" for key in keys
        ]

    axs: list[Axes] = scatter_base(
        Y,
        title=title,
        alpha=alpha,
        component_name=component_name,
        axis_labels=axis_labels,
        component_indexnames=components + 1,
        projection=projection,
        colors=color_ids,
        highlights=highlights,
        colorbars=colorbars,
        right_margin=right_margin,
        left_margin=left_margin,
        sizes=[size for _ in keys],
        markers=marker,
        color_map=color_map,
        show_ticks=show_ticks,
        ax=ax,
    )

    def add_centroid(centroids, name, Y, mask):
        Y_mask = Y[mask]
        if Y_mask.shape[0] == 0:
            return
        median = np.median(Y_mask, axis=0)
        i = np.argmin(np.sum(np.abs(Y_mask - median), axis=1))
        centroids[name] = Y_mask[i]

    # loop over all categorical annotation and plot it
    for ikey, pal in zip(categoricals, palettes, strict=False):
        key = keys[ikey]
        _utils.add_colors_for_categorical_sample_annotation(
            adata, key, palette=pal, force_update_colors=palette is not None
        )
        # actually plot the groups
        mask_remaining = np.ones(Y.shape[0], dtype=bool)
        centroids = {}
        if groups is None:
            for iname, name in enumerate(adata.obs[key].cat.categories):
                if name not in settings.categories_to_ignore:
                    mask = scatter_group(
                        axs[ikey],
                        key,
                        iname,
                        adata,
                        Y,
                        projection=projection,
                        size=size,
                        alpha=alpha,
                        marker=marker,
                    )
                    mask_remaining[mask] = False
                    if legend_loc.startswith("on data"):
                        add_centroid(centroids, name, Y, mask)
        else:
            groups = [groups] if isinstance(groups, str) else groups
            for name in groups:
                if name not in set(adata.obs[key].cat.categories):
                    msg = (
                        f"{name!r} is invalid! specify valid name, "
                        f"one of {adata.obs[key].cat.categories}"
                    )
                    raise ValueError(msg)
                else:
                    iname = np.flatnonzero(
                        adata.obs[key].cat.categories.values == name
                    )[0]
                    mask = scatter_group(
                        axs[ikey],
                        key,
                        iname,
                        adata,
                        Y,
                        projection=projection,
                        size=size,
                        alpha=alpha,
                        marker=marker,
                    )
                    if legend_loc.startswith("on data"):
                        add_centroid(centroids, name, Y, mask)
                    mask_remaining[mask] = False
        if mask_remaining.sum() > 0:
            data = [Y[mask_remaining, 0], Y[mask_remaining, 1]]
            if projection == "3d":
                data.append(Y[mask_remaining, 2])
            axs[ikey].scatter(
                *data,
                marker=marker,
                c="lightgrey",
                s=size,
                edgecolors="none",
                zorder=-1,
            )
        legend = None
        if legend_loc.startswith("on data"):
            if legend_fontweight is None:
                legend_fontweight = "bold"
            if legend_fontoutline is not None:
                path_effect = [
                    patheffects.withStroke(linewidth=legend_fontoutline, foreground="w")
                ]
            else:
                path_effect = None
            for name, pos in centroids.items():
                axs[ikey].text(
                    pos[0],
                    pos[1],
                    name,
                    weight=legend_fontweight,
                    verticalalignment="center",
                    horizontalalignment="center",
                    fontsize=legend_fontsize,
                    path_effects=path_effect,
                )

            all_pos = np.zeros((len(adata.obs[key].cat.categories), 2))
            for iname, name in enumerate(adata.obs[key].cat.categories):
                all_pos[iname] = centroids.get(name, [np.nan, np.nan])
            if legend_loc == "on data export":
                filename = settings.writedir / "pos.csv"
                logg.warning(f"exporting label positions to {filename}")
                settings.writedir.mkdir(parents=True, exist_ok=True)
                np.savetxt(filename, all_pos, delimiter=",")
        elif legend_loc == "right margin":
            legend = axs[ikey].legend(
                frameon=False,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                ncol=(
                    1
                    if len(adata.obs[key].cat.categories) <= 14
                    else 2
                    if len(adata.obs[key].cat.categories) <= 30
                    else 3
                ),
                fontsize=legend_fontsize,
            )
        elif legend_loc != "none":
            legend = axs[ikey].legend(
                frameon=False, loc=legend_loc, fontsize=legend_fontsize
            )
        if legend is not None:
            if Version(mpl.__version__) < Version("3.7"):
                _attr = "legendHandles"
            else:
                _attr = "legend_handles"
            for handle in getattr(legend, _attr):
                handle.set_sizes([300.0])

    # draw a frame around the scatter
    frameon = settings._frameon if frameon is None else frameon
    if not frameon and x is None and y is None:
        for ax_ in axs:
            ax_.set_xlabel("")
            ax_.set_ylabel("")
            ax_.set_frame_on(False)

    show = settings.autoshow if show is None else show
    _utils.savefig_or_show("scatter" if basis is None else basis, show=show, save=save)
    if show:
        return None
    if len(keys) > 1:
        return axs
    return axs[0]


@old_positionals(
    "dictionary",
    "indices",
    "labels",
    "color",
    "n_points",
    "log",
    "include_lowest",
    "show",
)
def ranking(  # noqa: PLR0912, PLR0913
    adata: AnnData,
    attr: Literal["var", "obs", "uns", "varm", "obsm"],
    keys: str | Sequence[str],
    *,
    dictionary: str | None = None,
    indices: Sequence[int] | None = None,
    labels: str | Sequence[str] | None = None,
    color: ColorLike = "black",
    n_points: int = 30,
    log: bool = False,
    include_lowest: bool = False,
    show: bool | None = None,
) -> gridspec.GridSpec | None:
    """Plot rankings.

    See, for example, how this is used in pl.pca_loadings.

    Parameters
    ----------
    adata
        The data.
    attr
        The attribute of AnnData that contains the score.
    keys
        The scores to look up an array from the attribute of adata.

    Returns
    -------
    Returns matplotlib gridspec with access to the axes.

    """
    if isinstance(keys, str) and indices is not None:
        scores = getattr(adata, attr)[keys][:, indices]
        keys = [f"{keys[:-1]}{i + 1}" for i in indices]
    elif dictionary is None:
        scores = getattr(adata, attr)[keys]
    else:
        scores = getattr(adata, attr)[dictionary][keys]
    n_panels = len(keys) if isinstance(keys, list) else 1
    if n_panels == 1:
        scores, keys = scores[:, None], [keys]
    if log:
        scores = np.log(scores)
    if labels is None:
        labels = (
            adata.var_names
            if attr in {"var", "varm"}
            else np.arange(scores.shape[0]).astype(str)
        )
    if isinstance(labels, str):
        labels = [labels + str(i + 1) for i in range(scores.shape[0])]
    if n_panels <= 5:
        n_rows, n_cols = 1, n_panels
    else:
        n_rows, n_cols = 2, int(n_panels / 2 + 0.5)
    _ = plt.figure(
        figsize=(
            n_cols * rcParams["figure.figsize"][0],
            n_rows * rcParams["figure.figsize"][1],
        )
    )
    left, bottom = 0.2 / n_cols, 0.13 / n_rows
    gs = gridspec.GridSpec(
        wspace=0.2,
        nrows=n_rows,
        ncols=n_cols,
        left=left,
        bottom=bottom,
        right=1 - (n_cols - 1) * left - 0.01 / n_cols,
        top=1 - (n_rows - 1) * bottom - 0.1 / n_rows,
    )
    for iscore, score in enumerate(scores.T):
        plt.subplot(gs[iscore])
        order_scores = np.argsort(score)[::-1]
        if not include_lowest:
            indices = order_scores[: n_points + 1]
        else:
            indices = order_scores[: n_points // 2]
            neg_indices = order_scores[-(n_points - (n_points // 2)) :]
        txt_args = dict(
            color=color,
            rotation="vertical",
            verticalalignment="bottom",
            horizontalalignment="center",
            fontsize=8,
        )
        for ig, g in enumerate(indices):
            plt.text(ig, score[g], labels[g], **txt_args)
        if include_lowest:
            score_mid = (score[g] + score[neg_indices[0]]) / 2
            if (len(indices) + len(neg_indices)) < len(order_scores):
                plt.text(len(indices), score_mid, "⋮", **txt_args)
                for ig, g in enumerate(neg_indices):
                    plt.text(ig + len(indices) + 2, score[g], labels[g], **txt_args)
            else:
                for ig, g in enumerate(neg_indices):
                    plt.text(ig + len(indices), score[g], labels[g], **txt_args)
            plt.xticks([])
        plt.title(keys[iscore].replace("_", " "))
        if n_panels <= 5 or iscore > n_cols:
            plt.xlabel("ranking")
        plt.xlim(-0.9, n_points + 0.9 + (1 if include_lowest else 0))
        score_min, score_max = (
            np.min(score[neg_indices if include_lowest else indices]),
            np.max(score[indices]),
        )
        plt.ylim(
            (0.95 if score_min > 0 else 1.05) * score_min,
            (1.05 if score_max > 0 else 0.95) * score_max,
        )
    show = settings.autoshow if show is None else show
    if show:
        return None
    return gs


@old_positionals(
    "log",
    "use_raw",
    "stripplot",
    "jitter",
    "size",
    "layer",
    "scale",
    "order",
    "multi_panel",
    "xlabel",
    "ylabel",
    "rotation",
    "show",
    "save",
    "ax",
)
@_doc_params(show_save_ax=doc_show_save_ax)
def violin(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    keys: str | Sequence[str],
    groupby: str | None = None,
    *,
    log: bool = False,
    use_raw: bool | None = None,
    stripplot: bool = True,
    jitter: float | bool = True,
    size: int = 1,
    layer: str | None = None,
    density_norm: DensityNorm = "width",
    order: Sequence[str] | None = None,
    multi_panel: bool | None = None,
    xlabel: str = "",
    ylabel: str | Sequence[str] | None = None,
    rotation: float | None = None,
    show: bool | None = None,
    save: bool | str | None = None,
    ax: Axes | None = None,
    # deprecatd
    scale: DensityNorm | Empty = _empty,
    **kwds,
) -> Axes | FacetGrid | None:
    """Violin plot.

    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    keys
        Keys for accessing variables of `.var_names` or fields of `.obs`.
    groupby
        The key of the observation grouping to consider.
    log
        Plot on logarithmic axis.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    layer
        Name of the AnnData object layer that wants to be plotted. By
        default adata.raw.X is plotted. If `use_raw=False` is set,
        then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    density_norm
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    order
        Order in which to show the categories.
    multi_panel
        Display keys in multiple panels also when `groupby is not None`.
    xlabel
        Label of the x axis. Defaults to `groupby` if `rotation` is `None`,
        otherwise, no label is shown.
    ylabel
        Label of the y axis. If `None` and `groupby` is `None`, defaults
        to `'value'`. If `None` and `groubpy` is not `None`, defaults to `keys`.
    rotation
        Rotation of xtick labels.
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.violinplot`.

    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object if `ax` is `None` else `None`.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.violin(adata, keys='S_score')

    Plot by category. Rotate x-axis labels so that they do not overlap.

    .. plot::
        :context: close-figs

        sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90)

    Set order of categories to be plotted or select specific categories to be plotted.

    .. plot::
        :context: close-figs

        groupby_order = ['CD34+', 'CD19+ B']
        sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90,
            order=groupby_order)

    Plot multiple keys.

    .. plot::
        :context: close-figs

        sc.pl.violin(adata, keys=['S_score', 'G2M_score'], groupby='bulk_labels',
            rotation=90)

    For large datasets consider omitting the overlaid scatter plot.

    .. plot::
        :context: close-figs

        sc.pl.violin(adata, keys='S_score', stripplot=False)

    .. currentmodule:: scanpy

    See Also
    --------
    pl.stacked_violin

    """
    import seaborn as sns  # Slow import, only import if called

    sanitize_anndata(adata)
    use_raw = _check_use_raw(adata, use_raw)
    if isinstance(keys, str):
        keys = [keys]
    keys = list(OrderedDict.fromkeys(keys))  # remove duplicates, preserving the order
    density_norm = _deprecated_scale(density_norm, scale, default="width")
    del scale

    if isinstance(ylabel, str | NoneType):
        ylabel = [ylabel] * (1 if groupby is None else len(keys))
    if groupby is None:
        if len(ylabel) != 1:
            msg = f"Expected number of y-labels to be `1`, found `{len(ylabel)}`."
            raise ValueError(msg)
    elif len(ylabel) != len(keys):
        msg = f"Expected number of y-labels to be `{len(keys)}`, found `{len(ylabel)}`."
        raise ValueError(msg)

    if groupby is not None:
        obs_df = get.obs_df(adata, keys=[groupby, *keys], layer=layer, use_raw=use_raw)
        if kwds.get("palette") is None:
            if not isinstance(adata.obs[groupby].dtype, CategoricalDtype):
                msg = (
                    f"The column `adata.obs[{groupby!r}]` needs to be categorical, "
                    f"but is of dtype {adata.obs[groupby].dtype}."
                )
                raise ValueError(msg)
            _utils.add_colors_for_categorical_sample_annotation(adata, groupby)
            kwds["hue"] = groupby
            kwds["palette"] = dict(
                zip(
                    obs_df[groupby].cat.categories,
                    adata.uns[f"{groupby}_colors"],
                    strict=True,
                )
            )
    else:
        obs_df = get.obs_df(adata, keys=keys, layer=layer, use_raw=use_raw)
    if groupby is None:
        obs_tidy = pd.melt(obs_df, value_vars=keys)
        x = "variable"
        ys = ["value"]
    else:
        obs_tidy = obs_df
        x = groupby
        ys = keys

    if multi_panel and groupby is None and len(ys) == 1:
        # This is a quick and dirty way for adapting scales across several
        # keys if groupby is None.
        y = ys[0]

        g: sns.axisgrid.FacetGrid = sns.catplot(
            y=y,
            data=obs_tidy,
            kind="violin",
            density_norm=density_norm,
            col=x,
            col_order=keys,
            sharey=False,
            cut=0,
            inner=None,
            **kwds,
        )

        if stripplot:
            grouped_df = obs_tidy.groupby(x, observed=True)
            for ax_id, key in zip(range(g.axes.shape[1]), keys, strict=True):
                sns.stripplot(
                    y=y,
                    data=grouped_df.get_group(key),
                    jitter=jitter,
                    size=size,
                    color="black",
                    ax=g.axes[0, ax_id],
                )
        if log:
            g.set(yscale="log")
        g.set_titles(col_template="{col_name}").set_xlabels("")
        if rotation is not None:
            for ax_base in g.axes[0]:
                ax_base.tick_params(axis="x", labelrotation=rotation)
    else:
        # set by default the violin plot cut=0 to limit the extend
        # of the violin plot (see stacked_violin code) for more info.
        kwds.setdefault("cut", 0)
        kwds.setdefault("inner")

        if ax is None:
            axs, _, _, _ = setup_axes(
                ax,
                panels=["x"] if groupby is None else keys,
                show_ticks=True,
                right_margin=0.3,
            )
        else:
            axs = [ax]
        for ax_base, y, ylab in zip(axs, ys, ylabel, strict=True):
            sns.violinplot(
                x=x,
                y=y,
                data=obs_tidy,
                order=order,
                orient="vertical",
                density_norm=density_norm,
                ax=ax_base,
                **kwds,
            )
            if stripplot:
                sns.stripplot(
                    x=x,
                    y=y,
                    data=obs_tidy,
                    order=order,
                    jitter=jitter,
                    color="black",
                    size=size,
                    ax=ax_base,
                )
            if xlabel == "" and groupby is not None and rotation is None:
                xlabel = groupby.replace("_", " ")
            ax_base.set_xlabel(xlabel)
            if ylab is not None:
                ax_base.set_ylabel(ylab)
            if log:
                ax_base.set_yscale("log")
            if rotation is not None:
                ax_base.tick_params(axis="x", labelrotation=rotation)
    show = settings.autoshow if show is None else show
    _utils.savefig_or_show("violin", show=show, save=save)
    if show:
        return None
    if multi_panel and groupby is None and len(ys) == 1:
        return g
    if len(axs) == 1:
        return axs[0]
    return axs


@old_positionals("use_raw", "show", "save")
@_doc_params(show_save_ax=doc_show_save_ax)
def clustermap(
    adata: AnnData,
    obs_keys: str | None = None,
    *,
    use_raw: bool | None = None,
    show: bool | None = None,
    save: bool | str | None = None,
    **kwds,
) -> ClusterGrid | None:
    """Hierarchically-clustered heatmap.

    Wraps :func:`seaborn.clustermap` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    obs_keys
        Categorical annotation to plot with a different color map.
        Currently, only a single key is supported.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    {show_save_ax}
    **kwds
        Keyword arguments passed to :func:`~seaborn.clustermap`.

    Returns
    -------
    If `show` is `False`, a :class:`~seaborn.matrix.ClusterGrid` object
    (see :func:`~seaborn.clustermap`).

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.krumsiek11()
        sc.pl.clustermap(adata)

    .. plot::
        :context: close-figs

        sc.pl.clustermap(adata, obs_keys='cell_type')

    """
    import seaborn as sns  # Slow import, only import if called

    if not isinstance(obs_keys, str | NoneType):
        msg = "Currently, only a single key is supported."
        raise ValueError(msg)
    sanitize_anndata(adata)
    use_raw = _check_use_raw(adata, use_raw)
    X = adata.raw.X if use_raw else adata.X
    if isinstance(X, CSBase):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    if obs_keys is not None:
        row_colors = adata.obs[obs_keys]
        _utils.add_colors_for_categorical_sample_annotation(adata, obs_keys)
        # do this more efficiently... just a quick solution
        lut = dict(
            zip(row_colors.cat.categories, adata.uns[obs_keys + "_colors"], strict=True)
        )
        row_colors = adata.obs[obs_keys].map(lut)
        g = sns.clustermap(df, row_colors=row_colors.values, **kwds)
    else:
        g = sns.clustermap(df, **kwds)
    show = settings.autoshow if show is None else show
    _utils.savefig_or_show("clustermap", show=show, save=save)
    if show:
        plt.show()
        return None
    return g


@old_positionals(
    "use_raw",
    "log",
    "num_categories",
    "dendrogram",
    "gene_symbols",
    "var_group_positions",
    "var_group_labels",
    "var_group_rotation",
    "layer",
    "standard_scale",
    "swap_axes",
    "show_gene_labels",
    "show",
    "save",
    "figsize",
    "vmin",
    "vmax",
    "vcenter",
    "norm",
)
@_doc_params(
    vminmax=doc_vboundnorm,
    show_save_ax=doc_show_save_ax,
    common_plot_args=doc_common_plot_args,
)
def heatmap(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str | Sequence[str],
    *,
    use_raw: bool | None = None,
    log: bool = False,
    num_categories: int = 7,
    dendrogram: bool | str = False,
    gene_symbols: str | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
    var_group_labels: Sequence[str] | None = None,
    var_group_rotation: float | None = None,
    layer: str | None = None,
    standard_scale: Literal["var", "obs"] | None = None,
    swap_axes: bool = False,
    show_gene_labels: bool | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    figsize: tuple[float, float] | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    vcenter: float | None = None,
    norm: Normalize | None = None,
    **kwds,
) -> dict[str, Axes] | None:
    """Heatmap of the expression values of genes.

    If `groupby` is given, the heatmap is ordered by the respective group. For
    example, a list of marker genes can be plotted, ordered by clustering. If
    the `groupby` observation annotation is not categorical the observation
    annotation is turned into a categorical by binning the data into the number
    specified in `num_categories`.

    Parameters
    ----------
    {common_plot_args}
    standard_scale
        Whether or not to standardize that dimension between 0 and 1, meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby`
         categories (if any). By setting `swap_axes` then x are the `groupby` categories and y the `var_names`.
    show_gene_labels
         By default gene labels are shown when there are 50 or less genes. Otherwise the labels are removed.
    {show_save_ax}
    {vminmax}
    **kwds
        Are passed to :func:`matplotlib.pyplot.imshow`.

    Returns
    -------
    Dict of :class:`~matplotlib.axes.Axes`

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.heatmap(adata, markers, groupby='bulk_labels', swap_axes=True)

    .. currentmodule:: scanpy

    See Also
    --------
    pl.rank_genes_groups_heatmap
    tl.rank_genes_groups

    """
    var_names, var_groups = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )
    del var_group_labels, var_group_positions

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw=use_raw,
        log=log,
        num_categories=num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    # check if var_group_labels are a subset of categories:
    var_groups_subset_of_groupby = var_groups is not None and set(
        var_groups.labels
    ).issubset(categories)

    if standard_scale == "obs":
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == "var":
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning("Unknown type for standard_scale, ignored")

    if groupby is None or len(categories) <= 1:
        categorical = False
        # dendrogram can only be computed  between groupby categories
        dendrogram = False
    else:
        categorical = True
        # get categories colors
        if isinstance(groupby, str) and isinstance(
            adata.obs[groupby].dtype, CategoricalDtype
        ):
            # saved category colors only work when groupby is valid adata.obs
            # categorical column. When groupby is a numerical column
            # or when groupby is a list of columns the colors are assigned on the fly,
            # which may create inconsistencies in multiple runs that require sorting
            # of the categories (eg. when dendrogram is plotted).
            if groupby + "_colors" not in adata.uns:
                # if colors are not found, assign a new palette
                # and save it using the same code for embeddings
                from ._tools.scatterplots import _get_palette

                _get_palette(adata, groupby)
            groupby_colors = adata.uns[groupby + "_colors"]
        else:
            # this case happen when adata.obs[groupby] is numeric
            # the values are converted into a category on the fly
            groupby_colors = None

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram_key=_dk(dendrogram),
            var_names=var_names,
            var_groups=var_groups,
            categories=categories,
        )

        var_groups = dendro_data["var_groups"]

        # reorder obs_tidy
        if dendro_data["var_names_idx_ordered"] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data["var_names_idx_ordered"]]
            var_names = [var_names[x] for x in dendro_data["var_names_idx_ordered"]]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data["categories_idx_ordered"]],
            ordered=True,
        )

        # reorder groupby colors
        if groupby_colors is not None:
            groupby_colors = [
                groupby_colors[x] for x in dendro_data["categories_idx_ordered"]
            ]

    if show_gene_labels is None:
        if len(var_names) <= 50:
            show_gene_labels = True
        else:
            show_gene_labels = False
            logg.warning(
                "Gene labels are not shown when more than 50 genes are visualized. "
                "To show gene labels set `show_gene_labels=True`"
            )
    if categorical:
        obs_tidy = obs_tidy.sort_index()

    colorbar_width = 0.2
    norm = check_colornorm(vmin, vmax, vcenter, norm)

    if not swap_axes:
        # define a layout of 2 rows x 4 columns
        # first row is for 'brackets' (if no brackets needed, the height of this row
        # is zero) second row is for main content. This second row is divided into
        # three axes:
        #   first ax is for the categories defined by `groupby`
        #   second ax is for the heatmap
        #   third ax is for the dendrogram
        #   fourth ax is for colorbar

        dendro_width = 1 if dendrogram else 0
        groupby_width = 0.2 if categorical else 0
        if figsize is None:
            height = 6
            heatmap_width = len(var_names) * 0.3 if show_gene_labels else 8
            width = heatmap_width + dendro_width + groupby_width
        else:
            width, height = figsize
            heatmap_width = width - (dendro_width + groupby_width)

        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = (0 if var_groups is None else 0.15, height)
        width_ratios = (
            groupby_width,
            heatmap_width,
            dendro_width,
            colorbar_width,
        )

        fig = plt.figure(figsize=(width, height))
        axs = gridspec.GridSpec(
            nrows=2,
            ncols=4,
            width_ratios=width_ratios,
            wspace=0.15 / width,
            hspace=0.13 / height,
            height_ratios=height_ratios,
        )

        heatmap_ax = fig.add_subplot(axs[1, 1])
        kwds.setdefault("interpolation", "nearest")
        im = heatmap_ax.imshow(obs_tidy.values, aspect="auto", norm=norm, **kwds)

        heatmap_ax.set_ylim(obs_tidy.shape[0] - 0.5, -0.5)
        heatmap_ax.set_xlim(-0.5, obs_tidy.shape[1] - 0.5)
        heatmap_ax.tick_params(axis="y", left=False, labelleft=False)
        heatmap_ax.set_ylabel("")
        heatmap_ax.grid(visible=False)

        if show_gene_labels:
            heatmap_ax.tick_params(axis="x", labelsize="small")
            heatmap_ax.set_xticks(np.arange(len(var_names)))
            heatmap_ax.set_xticklabels(var_names, rotation=90)
        else:
            heatmap_ax.tick_params(axis="x", labelbottom=False, bottom=False)
        # plot colorbar
        _plot_colorbar(im, fig, axs[1, 3])

        if categorical:
            groupby_ax = fig.add_subplot(axs[1, 0])
            (
                label2code,
                ticks,
                labels,
                groupby_cmap,
                norm,
            ) = _plot_categories_as_colorblocks(
                groupby_ax, obs_tidy, colors=groupby_colors, orientation="left"
            )

            # add lines to main heatmap
            line_positions = (
                np.cumsum(obs_tidy.index.value_counts(sort=False))[:-1] - 0.5
            )
            heatmap_ax.hlines(
                line_positions,
                -0.5,
                len(var_names) - 0.5,
                lw=1,
                color="black",
                zorder=10,
                clip_on=False,
            )

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[1, 2], sharey=heatmap_ax)
            _plot_dendrogram(
                dendro_ax, adata, groupby, dendrogram_key=_dk(dendrogram), ticks=ticks
            )

        # plot group legends on top of heatmap_ax (if given)
        if var_groups is not None:
            gene_groups_ax = fig.add_subplot(axs[0, 1], sharex=heatmap_ax)
            _plot_var_groups_brackets(
                gene_groups_ax,
                var_groups=var_groups,
                rotation=var_group_rotation,
                left_adjustment=-0.3,
                right_adjustment=0.3,
            )

    # swap axes case
    else:
        # define a layout of 3 rows x 3 columns
        # The first row is for the dendrogram (if not dendrogram height is zero)
        # second row is for main content. This col is divided into three axes:
        #   first ax is for the heatmap
        #   second ax is for 'brackets' if any (othwerise width is zero)
        #   third ax is for colorbar

        dendro_height = 0.8 if dendrogram else 0
        groupby_height = 0.13 if categorical else 0
        if figsize is None:
            heatmap_height = len(var_names) * 0.18 if show_gene_labels else 4
            width = 10
            height = heatmap_height + dendro_height + groupby_height
        else:
            width, height = figsize
            heatmap_height = height - (dendro_height + groupby_height)

        height_ratios = (dendro_height, heatmap_height, groupby_height)
        # add some space in case 'brackets' want to be plotted on top of the image
        width_ratios = (width, 0 if var_groups is None else 0.14, colorbar_width)

        fig = plt.figure(figsize=(width, height))
        axs = gridspec.GridSpec(
            nrows=3,
            ncols=3,
            wspace=0.25 / width,
            hspace=0.3 / height,
            width_ratios=width_ratios,
            height_ratios=height_ratios,
        )

        # plot heatmap
        heatmap_ax = fig.add_subplot(axs[1, 0])

        kwds.setdefault("interpolation", "nearest")
        im = heatmap_ax.imshow(obs_tidy.T.values, aspect="auto", norm=norm, **kwds)
        heatmap_ax.set_xlim(0 - 0.5, obs_tidy.shape[0] - 0.5)
        heatmap_ax.set_ylim(obs_tidy.shape[1] - 0.5, -0.5)
        heatmap_ax.tick_params(axis="x", bottom=False, labelbottom=False)
        heatmap_ax.set_xlabel("")
        heatmap_ax.grid(visible=False)
        if show_gene_labels:
            heatmap_ax.tick_params(axis="y", labelsize="small", length=1)
            heatmap_ax.set_yticks(np.arange(len(var_names)))
            heatmap_ax.set_yticklabels(var_names, rotation=0)
        else:
            heatmap_ax.tick_params(axis="y", labelleft=False, left=False)

        if categorical:
            groupby_ax = fig.add_subplot(axs[2, 0])
            (
                label2code,
                ticks,
                labels,
                groupby_cmap,
                norm,
            ) = _plot_categories_as_colorblocks(
                groupby_ax, obs_tidy, colors=groupby_colors, orientation="bottom"
            )
            # add lines to main heatmap
            line_positions = (
                np.cumsum(obs_tidy.index.value_counts(sort=False))[:-1] - 0.5
            )
            heatmap_ax.vlines(
                line_positions,
                -0.5,
                len(var_names) - 0.5,
                lw=1,
                color="black",
                zorder=10,
                clip_on=False,
            )

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[0, 0], sharex=heatmap_ax)
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                dendrogram_key=_dk(dendrogram),
                ticks=ticks,
                orientation="top",
            )

        # plot group legends next to the heatmap_ax (if given)
        if var_groups is not None:
            gene_groups_ax = fig.add_subplot(axs[1, 1])
            arr = []
            for idx, (label, pos) in enumerate(zip(*var_groups, strict=True)):
                label_code = label2code[label] if var_groups_subset_of_groupby else idx
                arr += [label_code] * (pos[1] + 1 - pos[0])
            gene_groups_ax.imshow(
                np.array([arr]).T, aspect="auto", cmap=groupby_cmap, norm=norm
            )
            gene_groups_ax.axis("off")

        # plot colorbar
        _plot_colorbar(im, fig, axs[1, 2])

    return_ax_dict = {"heatmap_ax": heatmap_ax}
    if categorical:
        return_ax_dict["groupby_ax"] = groupby_ax
    if dendrogram:
        return_ax_dict["dendrogram_ax"] = dendro_ax
    if var_groups is not None:
        return_ax_dict["gene_groups_ax"] = gene_groups_ax

    _utils.savefig_or_show("heatmap", show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return return_ax_dict


@old_positionals(
    "use_raw",
    "log",
    "dendrogram",
    "gene_symbols",
    "var_group_positions",
    "var_group_labels",
    "layer",
    "show",
    "save",
    "figsize",
)
@_doc_params(show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args)
def tracksplot(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str,
    *,
    use_raw: bool | None = None,
    log: bool = False,
    dendrogram: bool | str = False,
    gene_symbols: str | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
    var_group_labels: Sequence[str] | None = None,
    layer: str | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    figsize: tuple[float, float] | None = None,
    **kwds,
) -> dict[str, Axes] | None:
    """Compact plot of expression of a list of genes.

    In this type of plot each var_name is plotted as a filled line plot where the
    y values correspond to the var_name values and x is each of the cells. Best results
    are obtained when using raw counts that are not log.

    `groupby` is required to sort and order the values using the respective group
    and should be a categorical value.

    Parameters
    ----------
    {common_plot_args}
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.heatmap`.

    Returns
    -------
    A list of :class:`~matplotlib.axes.Axes`.

    Examples
    --------
    Using var_names as list:

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.tracksplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    .. plot::
        :context: close-figs

        markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
        sc.pl.tracksplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    .. currentmodule:: scanpy

    See Also
    --------
    pl.rank_genes_groups_tracksplot: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.

    """
    if groupby not in adata.obs_keys() or adata.obs[groupby].dtype.name != "category":
        msg = (
            "groupby has to be a valid categorical observation. "
            f"Given value: {groupby}, valid categorical observations: "
            f"{[x for x in adata.obs_keys() if adata.obs[x].dtype.name == 'category']}"
        )
        raise ValueError(msg)

    var_names, var_groups = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )
    del var_group_labels, var_group_positions

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw=use_raw,
        log=log,
        num_categories=None,  # TODO: fix this line
        gene_symbols=gene_symbols,
        layer=layer,
    )

    # get categories colors:
    if f"{groupby}_colors" not in adata.uns:
        from ._utils import _set_default_colors_for_categorical_obs

        _set_default_colors_for_categorical_obs(adata, groupby)
    groupby_colors = adata.uns[groupby + "_colors"]

    if dendrogram:
        # compute dendrogram if needed and reorder
        # rows and columns to match leaves order.
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram_key=_dk(dendrogram),
            var_names=var_names,
            var_groups=var_groups,
            categories=categories,
        )
        # reorder obs_tidy
        if dendro_data["var_names_idx_ordered"] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data["var_names_idx_ordered"]]
            var_names = [var_names[x] for x in dendro_data["var_names_idx_ordered"]]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data["categories_idx_ordered"]],
            ordered=True,
        )
        categories = [categories[x] for x in dendro_data["categories_idx_ordered"]]

        groupby_colors = [
            groupby_colors[x] for x in dendro_data["categories_idx_ordered"]
        ]

    obs_tidy = obs_tidy.sort_index()

    # obtain the start and end of each category and make
    # a list of ranges that will be used to plot a different
    # color
    cumsum = [0, *np.cumsum(obs_tidy.index.value_counts(sort=False)).tolist()]
    x_values = list(pairwise(cumsum))

    dendro_height = 1 if dendrogram else 0

    groupby_height = 0.24
    # +2 because of dendrogram on top and categories at bottom
    num_rows = len(var_names) + 2
    if figsize is None:
        width = 12
        track_height = 0.25
    else:
        width, height = figsize
        track_height = (height - (dendro_height + groupby_height)) / len(var_names)

    height_ratios = [dendro_height] + [track_height] * len(var_names) + [groupby_height]
    height = sum(height_ratios)

    obs_tidy = obs_tidy.T

    fig = plt.figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        ncols=2,
        nrows=num_rows,
        wspace=1.0 / width,
        hspace=0,
        height_ratios=height_ratios,
        width_ratios=[width, 0.14],
    )
    axs_list = []
    first_ax = None
    for idx, var in enumerate(var_names):
        ax_idx = idx + 1  # this is because of the dendrogram
        if first_ax is None:
            ax = fig.add_subplot(axs[ax_idx, 0])
            first_ax = ax
        else:
            ax = fig.add_subplot(axs[ax_idx, 0], sharex=first_ax)
        axs_list.append(ax)
        for cat_idx, _category in enumerate(categories):
            x_start, x_end = x_values[cat_idx]
            ax.fill_between(
                range(x_start, x_end),
                0,
                obs_tidy.iloc[idx, x_start:x_end],
                lw=0.1,
                color=groupby_colors[cat_idx],
            )

        # remove the xticks labels except for the last processed plot.
        # Because the plots share the x axis it is redundant and less compact
        # to plot the axis for each plot
        if idx < len(var_names) - 1:
            ax.tick_params(labelbottom=False, labeltop=False, bottom=False, top=False)
            ax.set_xlabel("")
        if log:
            ax.set_yscale("log")
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.grid(visible=False)
        ymin, ymax = ax.get_ylim()
        ymax = int(ymax)
        ax.set_yticks([ymax])
        ax.set_yticklabels([str(ymax)], ha="left", va="top")
        ax.spines["right"].set_position(("axes", 1.01))
        ax.tick_params(
            axis="y",
            labelsize="x-small",
            right=True,
            left=False,
            length=2,
            which="both",
            labelright=True,
            labelleft=False,
            direction="in",
        )
        ax.set_ylabel(var, rotation=0, fontsize="small", ha="right", va="bottom")
        ax.yaxis.set_label_coords(-0.005, 0.1)
    ax.set_xlim(0, x_end)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)

    # the ax to plot the groupby categories is split to add a small space
    # between the rest of the plot and the categories
    axs2 = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=axs[num_rows - 1, 0], height_ratios=[1, 1]
    )

    groupby_ax = fig.add_subplot(axs2[1])

    label2code, ticks, labels, groupby_cmap, norm = _plot_categories_as_colorblocks(
        groupby_ax, obs_tidy.T, colors=groupby_colors, orientation="bottom"
    )
    # add lines to plot
    overlay_ax = fig.add_subplot(axs[1:-1, 0], sharex=first_ax)
    line_positions = np.cumsum(obs_tidy.T.index.value_counts(sort=False))[:-1]
    overlay_ax.vlines(line_positions, 0, 1, lw=0.5, linestyle="--")
    overlay_ax.axis("off")
    overlay_ax.set_ylim(0, 1)

    if dendrogram:
        dendro_ax = fig.add_subplot(axs[0], sharex=first_ax)
        _plot_dendrogram(
            dendro_ax,
            adata,
            groupby,
            dendrogram_key=_dk(dendrogram),
            orientation="top",
            ticks=ticks,
        )

    if var_groups is not None:
        gene_groups_ax = fig.add_subplot(axs[1:-1, 1])
        arr = []
        for idx, pos in enumerate(var_groups.positions):
            arr += [idx] * (pos[1] + 1 - pos[0])

        gene_groups_ax.imshow(
            np.array([arr]).T, aspect="auto", cmap=groupby_cmap, norm=norm
        )
        gene_groups_ax.axis("off")

    return_ax_dict = {"track_axes": axs_list, "groupby_ax": groupby_ax}
    if dendrogram:
        return_ax_dict["dendrogram_ax"] = dendro_ax
    if var_groups is not None:
        return_ax_dict["gene_groups_ax"] = gene_groups_ax

    _utils.savefig_or_show("tracksplot", show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return return_ax_dict


@_doc_params(show_save_ax=doc_show_save_ax)
def dendrogram(
    adata: AnnData,
    groupby: str,
    *,
    dendrogram_key: str | None = None,
    orientation: Literal["top", "bottom", "left", "right"] = "top",
    remove_labels: bool = False,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot a dendrogram of the categories defined in `groupby`.

    See :func:`~scanpy.tl.dendrogram`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        Categorical data column used to create the dendrogram
    dendrogram_key
        Key under with the dendrogram information was stored.
        By default the dendrogram information is stored under
        `.uns[f'dendrogram_{{groupby}}']`.
    orientation
        Origin of the tree. Will grow into the opposite direction.
    remove_labels
        Don’t draw labels. Used e.g. by :func:`scanpy.pl.matrixplot`
        to annotate matrix columns/rows.
    {show_save_ax}

    Returns
    -------
    :class:`matplotlib.axes.Axes`

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.dendrogram(adata, 'bulk_labels')
        sc.pl.dendrogram(adata, 'bulk_labels')

    .. currentmodule:: scanpy

    """
    if ax is None:
        _, ax = plt.subplots()
    _plot_dendrogram(
        ax,
        adata,
        groupby,
        dendrogram_key=dendrogram_key,
        remove_labels=remove_labels,
        orientation=orientation,
    )
    _utils.savefig_or_show("dendrogram", show=show, save=save)
    return ax


@old_positionals(
    "show_correlation_numbers",
    "dendrogram",
    "figsize",
    "show",
    "save",
    "ax",
    "vmin",
    "vmax",
    "vcenter",
    "norm",
)
@_doc_params(show_save_ax=doc_show_save_ax, vminmax=doc_vboundnorm)
def correlation_matrix(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    groupby: str,
    *,
    show_correlation_numbers: bool = False,
    dendrogram: bool | str | None = None,
    figsize: tuple[float, float] | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: Axes | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    vcenter: float | None = None,
    norm: Normalize | None = None,
    **kwds,
) -> list[Axes] | None:
    """Plot the correlation matrix computed as part of `sc.tl.dendrogram`.

    Parameters
    ----------
    adata
    groupby
        Categorical data column used to create the dendrogram
    show_correlation_numbers
        If `show_correlation=True`, plot the correlation on top of each cell.
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the
        hierarchical clustering between the `groupby` categories is added.
        The dendrogram is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously,
        the function is called with default parameters.
    figsize
        By default a figure size that aims to produce a squared correlation
        matrix plot is used. Format is (width, height)
    {show_save_ax}
    {vminmax}
    **kwds
        Only if `show_correlation` is True:
        Are passed to :func:`matplotlib.pyplot.pcolormesh` when plotting the
        correlation heatmap. `cmap` can be used to change the color palette.

    Returns
    -------
    If `show=False`, returns a list of :class:`matplotlib.axes.Axes` objects.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, "bulk_labels")
    >>> sc.pl.correlation_matrix(adata, "bulk_labels")

    """
    dendrogram_key = _get_dendrogram_key(adata, _dk(dendrogram), groupby)

    index = adata.uns[dendrogram_key]["categories_idx_ordered"]
    corr_matrix = adata.uns[dendrogram_key]["correlation_matrix"]
    # reorder matrix columns according to the dendrogram
    if dendrogram is None:
        dendrogram = ax is None
    if dendrogram:
        if ax is not None:
            msg = "Can only plot dendrogram when not plotting to an axis"
            raise ValueError(msg)
        assert (len(index)) == corr_matrix.shape[0]
        corr_matrix = corr_matrix[index, :]
        corr_matrix = corr_matrix[:, index]
        labels = list(adata.obs[groupby].cat.categories)
        labels = np.array(labels).astype("str")[index]
    else:
        labels = adata.obs[groupby].cat.categories
    num_rows = corr_matrix.shape[0]
    colorbar_height = 0.2
    dendrogram_width = 1.8 if dendrogram else 0
    if figsize is None:
        corr_matrix_height = num_rows * 0.6
        height = corr_matrix_height + colorbar_height
        width = corr_matrix_height + dendrogram_width
    else:
        width, height = figsize
        corr_matrix_height = height - colorbar_height

    fig = plt.figure(figsize=(width, height)) if ax is None else None
    # layout with 2 rows and 2  columns:
    # row 1: dendrogram + correlation matrix
    # row 2: nothing + colormap bar (horizontal)
    gs = gridspec.GridSpec(
        nrows=2,
        ncols=2,
        width_ratios=[dendrogram_width, corr_matrix_height],
        height_ratios=[corr_matrix_height, colorbar_height],
        wspace=0.01,
        hspace=0.05,
    )

    axs = []
    corr_matrix_ax = fig.add_subplot(gs[1]) if ax is None else ax
    if dendrogram:
        dendro_ax = fig.add_subplot(gs[0], sharey=corr_matrix_ax)
        _plot_dendrogram(
            dendro_ax,
            adata,
            groupby,
            dendrogram_key=dendrogram_key,
            remove_labels=True,
            orientation="left",
            ticks=np.arange(corr_matrix.shape[0]) + 0.5,
        )
        axs.append(dendro_ax)
    # define some default pcolormesh parameters
    if "edgecolors" not in kwds:
        if corr_matrix.shape[0] > 30:
            # when there are too many rows it is better to remove
            # the black lines surrounding the boxes in the heatmap
            kwds["edgecolors"] = "none"
        else:
            kwds["edgecolors"] = "black"
            kwds.setdefault("linewidth", 0.01)
    if vmax is None and vmin is None and norm is None:
        vmax = 1
        vmin = -1
    norm = check_colornorm(vmin, vmax, vcenter, norm)
    if "cmap" not in kwds:
        # by default use a divergent color map
        kwds["cmap"] = "bwr"

    img_mat = corr_matrix_ax.pcolormesh(corr_matrix, norm=norm, **kwds)
    corr_matrix_ax.set_xlim(0, num_rows)
    corr_matrix_ax.set_ylim(0, num_rows)

    corr_matrix_ax.yaxis.tick_right()
    corr_matrix_ax.set_yticks(np.arange(corr_matrix.shape[0]) + 0.5)
    corr_matrix_ax.set_yticklabels(labels)

    corr_matrix_ax.xaxis.set_tick_params(labeltop=True)
    corr_matrix_ax.xaxis.set_tick_params(labelbottom=False)
    corr_matrix_ax.set_xticks(np.arange(corr_matrix.shape[0]) + 0.5)
    corr_matrix_ax.set_xticklabels(labels, rotation=45, ha="left")

    for ax_name in "xy":
        corr_matrix_ax.tick_params(axis=ax_name, which="both", bottom=False, top=False)

    if show_correlation_numbers:
        for row, col in product(range(num_rows), repeat=2):
            corr_matrix_ax.text(
                row + 0.5,
                col + 0.5,
                f"{corr_matrix[row, col]:.2f}",
                ha="center",
                va="center",
            )

    axs.append(corr_matrix_ax)

    if ax is None:  # Plot colorbar
        colormap_ax = fig.add_subplot(gs[3])
        cobar = plt.colorbar(img_mat, cax=colormap_ax, orientation="horizontal")
        cobar.solids.set_edgecolor("face")
        axs.append(colormap_ax)

    show = settings.autoshow if show is None else show
    _utils.savefig_or_show("correlation_matrix", show=show, save=save)
    if ax is not None or show:
        return None
    return axs


def _prepare_dataframe(  # noqa: PLR0912
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str | Sequence[str] | None = None,
    *,
    use_raw: bool | None = None,
    log: bool = False,
    num_categories: int = 7,
    layer: str | None = None,
    gene_symbols: str | None = None,
) -> tuple[Sequence[str], pd.DataFrame]:
    """Prepare a data frame of categories (`groupby`) × `var_names`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of  `adata.var_names`.
    groupby
        The key of the observation grouping to consider. It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories`.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    log
        Use the log of the values.
    layer
        AnnData layer to use. Takes precedence over `use_raw`
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    gene_symbols
        Key for field in .var that stores gene symbols.

    Returns
    -------
    Tuple of `pandas.DataFrame` and list of categories.

    """
    sanitize_anndata(adata)
    use_raw = _check_use_raw(adata, use_raw, layer=layer)
    if isinstance(var_names, str):
        var_names = [var_names]

    groupby_index = None
    if groupby is not None:
        if isinstance(groupby, str):
            # if not a list, turn into a list
            groupby = [groupby]
        for group in groupby:
            if group not in [*adata.obs_keys(), adata.obs.index.name]:
                if adata.obs.index.name is not None:
                    msg = f' or index name "{adata.obs.index.name}"'
                else:
                    msg = ""
                raise ValueError(
                    "groupby has to be a valid observation. "
                    f"Given {group}, is not in observations: {adata.obs_keys()}" + msg
                )
            if group in adata.obs.columns and group == adata.obs.index.name:
                msg = (
                    f"Given group {group} is both and index and a column level, "
                    "which is ambiguous."
                )
                raise ValueError(msg)
            if group == adata.obs.index.name:
                groupby_index = group
    if groupby_index is not None:
        # obs_tidy contains adata.obs.index
        # and does not need to be given
        groupby = groupby.copy()  # copy to not modify user passed parameter
        groupby.remove(groupby_index)
    keys = list(groupby) + list(np.unique(var_names))
    obs_tidy = get.obs_df(
        adata, keys=keys, layer=layer, use_raw=use_raw, gene_symbols=gene_symbols
    )
    assert np.all(np.array(keys) == np.array(obs_tidy.columns))

    if groupby_index is not None:
        # reset index to treat all columns the same way.
        obs_tidy.reset_index(inplace=True)
        groupby.append(groupby_index)

    if groupby is None:
        categorical = pd.Series(np.repeat("", len(obs_tidy))).astype("category")
    elif len(groupby) == 1 and is_numeric_dtype(obs_tidy[groupby[0]]):
        # if the groupby column is not categorical, turn it into one
        # by subdividing into  `num_categories` categories
        categorical = pd.cut(obs_tidy[groupby[0]], num_categories)
    elif len(groupby) == 1:
        categorical = obs_tidy[groupby[0]].astype("category")
        categorical.name = groupby[0]
    else:
        # join the groupby values  using "_" to make a new 'category'
        categorical = obs_tidy[groupby].apply("_".join, axis=1).astype("category")
        categorical.name = "_".join(groupby)

        # preserve category order
        from itertools import product

        order = {
            "_".join(k): idx
            for idx, k in enumerate(
                product(*(obs_tidy[g].cat.categories for g in groupby))
            )
        }
        categorical = categorical.cat.reorder_categories(
            sorted(categorical.cat.categories, key=lambda x: order[x])
        )
    obs_tidy = obs_tidy[var_names].set_index(categorical)
    categories = obs_tidy.index.categories

    if log:
        obs_tidy = np.log1p(obs_tidy)

    return categories, obs_tidy


def _plot_var_groups_brackets(
    var_groups_ax: Axes,
    *,
    var_groups: VarGroups,
    left_adjustment: float = -0.3,
    right_adjustment: float = 0.3,
    rotation: float | None = None,
    orientation: Literal["top", "right"] = "top",
    wide: bool = False,
) -> None:
    """Draw brackets that represent groups of genes on the give axis.

    For best results, this axis is located on top of an image whose
    x axis contains gene names.

    The `var_groups_ax` should share the x axis with the main ax.

    E.g: `var_groups_ax=fig.add_subplot(axs[0, 0], sharex=dot_ax)`

    This function is used by dotplot, heatmap etc.

    Parameters
    ----------
    var_groups_ax
        In this axis the gene marks are drawn
    var_groups
        Group labels and positions
    left_adjustment
        adjustment to plot the bracket start slightly before or after the first gene position.
        If the value is negative the start is moved before.
    right_adjustment
        adjustment to plot the bracket end slightly before or after the last gene position
        If the value is negative the start is moved before.
    rotation
        rotation degrees for the labels. If not given, small labels (<4 characters) are not
        rotated, otherwise, they are rotated 90 degrees
    orientation
        location of the brackets. Either `top` or `right`

    Returns
    -------
    None

    """
    from matplotlib import patches
    from matplotlib.path import Path

    # get the 'brackets' coordinates as lists of start and end positions
    left = [x[0] + left_adjustment for x in var_groups.positions]
    right = [x[1] + right_adjustment for x in var_groups.positions]

    # verts and codes are used by PathPatch to make the brackets
    verts = []
    codes = []
    if orientation == "top":
        # rotate labels if any of them is longer than 4 characters
        if rotation is None:
            rotation = 90 if max([len(x) for x in var_groups.labels]) > 4 else 0
        for idx in range(len(left)):
            verts.append((left[idx], 0))  # lower-left
            verts.append((left[idx], 0.6))  # upper-left
            verts.append((right[idx], 0.6))  # upper-right
            verts.append((right[idx], 0))  # lower-right

            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            group_x_center = left[idx] + float(right[idx] - left[idx]) / 2
            var_groups_ax.text(
                group_x_center,
                1.1,
                var_groups.labels[idx],
                ha="center",
                va="bottom",
                rotation=rotation,
            )
    else:
        top = left
        bottom = right
        for idx in range(len(top)):
            verts.append((0, top[idx]))  # upper-left
            verts.append((0.4 if wide else 0.15, top[idx]))  # upper-right
            verts.append((0.4 if wide else 0.15, bottom[idx]))  # lower-right
            verts.append((0, bottom[idx]))  # lower-left

            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            diff = bottom[idx] - top[idx]
            group_y_center = top[idx] + float(diff) / 2
            # cut label to fit available space
            label = (
                var_groups.labels[idx][: int(diff * 2)] + "."
                if diff * 2 < len(var_groups.labels[idx])
                else var_groups.labels[idx]
            )
            var_groups_ax.text(
                1.1 if wide else 0.6,
                group_y_center,
                label,
                ha="right",
                va="center",
                rotation=270,
                fontsize="small",
            )

    path = Path(verts, codes)

    patch = patches.PathPatch(path, facecolor="none", lw=1.5)

    var_groups_ax.add_patch(patch)
    var_groups_ax.grid(visible=False)
    var_groups_ax.axis("off")
    # remove y ticks
    var_groups_ax.tick_params(axis="y", left=False, labelleft=False)
    # remove x ticks and labels
    var_groups_ax.tick_params(axis="x", bottom=False, labelbottom=False, labeltop=False)


class _ReorderCats(TypedDict):
    categories_idx_ordered: Sequence[int]
    categories_ordered: Sequence[str]
    var_names_idx_ordered: Sequence[int] | None
    var_names_ordered: Sequence[str] | None
    var_groups: VarGroups | None


def _reorder_categories_after_dendrogram(
    adata: AnnData,
    groupby: str | Sequence[str],
    *,
    dendrogram_key: str | None,
    var_names: Sequence[str],
    var_groups: VarGroups | None,
    categories: Sequence[str],
) -> _ReorderCats:
    """Reorder the the groupby observations based on the dendrogram results.

    The function checks if a dendrogram has already been precomputed.
    If not, `sc.tl.dendrogram` is run with default parameters.

    The results found in `.uns[dendrogram_key]` are used to reorder `var_groups`.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    dendro_info = adata.uns[
        _get_dendrogram_key(adata, dendrogram_key, groupby, validate_groupby=True)
    ]

    if categories is None:
        categories = adata.obs[dendro_info["groupby"]].cat.categories

    # order of groupby categories
    categories_idx_ordered = dendro_info["categories_idx_ordered"]
    categories_ordered = dendro_info["categories_ordered"]

    if len(categories) != len(categories_idx_ordered):
        msg = (
            "Incompatible observations. Dendrogram data has "
            f"{len(categories_idx_ordered)} categories but current groupby "
            f"observation {groupby!r} contains {len(categories)} categories. "
            "Most likely the underlying groupby observation changed after the "
            "initial computation of `sc.tl.dendrogram`. "
            "Please run `sc.tl.dendrogram` again.'"
        )
        raise ValueError(msg)

    # reorder var_groups (if any)
    if var_groups is None:
        var_names_idx_ordered = None
    elif set(var_groups.labels) == set(categories):
        positions_ordered = []
        labels_ordered = []
        position_start = 0
        var_names_idx_ordered = []
        for cat_name in categories_ordered:
            idx = var_groups.labels.index(cat_name)
            position = var_groups.positions[idx]
            _var_names = var_names[position[0] : position[1] + 1]
            var_names_idx_ordered.extend(range(position[0], position[1] + 1))
            positions_ordered.append(
                (position_start, position_start + len(_var_names) - 1)
            )
            position_start += len(_var_names)
            labels_ordered.append(var_groups.labels[idx])
        var_groups = VarGroups(labels_ordered, positions_ordered)
    else:
        logg.warning(
            "Groups are not reordered because the `groupby` categories "
            "and the `var_group_labels` are different.\n"
            f"categories: {_format_first_three_categories(categories)}\n"
            f"var_group_labels: {_format_first_three_categories(var_groups.labels)}"
        )
        var_names_idx_ordered = list(range(len(var_names)))

    if var_names_idx_ordered is not None:
        var_names_ordered = [var_names[x] for x in var_names_idx_ordered]
    else:
        var_names_ordered = None

    return _ReorderCats(
        categories_idx_ordered=categories_idx_ordered,
        categories_ordered=dendro_info["categories_ordered"],
        var_names_idx_ordered=var_names_idx_ordered,
        var_names_ordered=var_names_ordered,
        var_groups=var_groups,
    )


def _format_first_three_categories(categories):
    """Clean up warning message."""
    categories = list(categories)
    if len(categories) > 3:
        categories = [*categories[:3], "etc."]
    return ", ".join(categories)


def _get_dendrogram_key(
    adata: AnnData,
    dendrogram_key: str | None,
    groupby: str | Sequence[str],
    *,
    validate_groupby: bool = False,
) -> str:
    # the `dendrogram_key` can be a bool an NoneType or the name of the
    # dendrogram key. By default the name of the dendrogram key is 'dendrogram'
    if dendrogram_key is None:
        if isinstance(groupby, str):
            dendrogram_key = f"dendrogram_{groupby}"
        elif isinstance(groupby, Sequence):
            dendrogram_key = f"dendrogram_{'_'.join(groupby)}"
        else:
            msg = f"groupby has wrong type: {type(groupby).__name__}."
            raise AssertionError(msg)

    if dendrogram_key not in adata.uns:
        from ..tools._dendrogram import dendrogram

        logg.warning(
            f"dendrogram data not found (using key={dendrogram_key}). "
            "Running `sc.tl.dendrogram` with default parameters. For fine "
            "tuning it is recommended to run `sc.tl.dendrogram` independently."
        )
        dendrogram(adata, groupby, key_added=dendrogram_key)

    if "dendrogram_info" not in adata.uns[dendrogram_key]:
        msg = (
            f"The given dendrogram key ({dendrogram_key!r}) does not contain "
            "valid dendrogram information."
        )
        raise ValueError(msg)

    if validate_groupby:
        existing_groupby = adata.uns[dendrogram_key]["groupby"]
        if groupby != existing_groupby:
            msg = (
                "Incompatible observations. The precomputed dendrogram contains "
                f"information for the observation: {groupby!r} while the plot is "
                f"made for the observation: {existing_groupby!r}. "
                "Please run `sc.tl.dendrogram` using the right observation.'"
            )
            raise ValueError(msg)

    return dendrogram_key


def _plot_dendrogram(  # noqa: PLR0915
    dendro_ax: Axes,
    adata: AnnData,
    groupby: str | Sequence[str],
    *,
    dendrogram_key: str | None = None,
    orientation: Literal["top", "bottom", "left", "right"] = "right",
    remove_labels: bool = True,
    ticks: Collection[float] | None = None,
):
    """Plot a dendrogram on the given ax.

    Uses the precomputed dendrogram information stored in `.uns[dendrogram_key]`.
    """
    dendrogram_key = _get_dendrogram_key(adata, dendrogram_key, groupby)

    def translate_pos(pos_list, new_ticks, old_ticks):
        """Transform the dendrogram coordinates to a given new position.

        The xlabel_pos and orig_ticks should be of the same length.

        This is mostly done for the heatmap case, where the position of the
        dendrogram leaves needs to be adjusted depending on the category size.

        Parameters
        ----------
        pos_list
            list of dendrogram positions that should be translated
        new_ticks
            sorted list of goal tick positions (e.g. [0,1,2,3] )
        old_ticks
            sorted list of original tick positions (e.g. [5, 15, 25, 35]),
            This list is usually the default position used by
            `scipy.cluster.hierarchy.dendrogram`.

        Returns
        -------
        translated list of positions

        Examples
        --------
        >>> translate_pos(
        ...     [5, 15, 20, 21],
        ...     [0, 1, 2, 3],
        ...     [5, 15, 25, 35],
        ... )
        [0, 1, 1.5, 1.6]

        """
        # of given coordinates.

        if not isinstance(old_ticks, list):
            # assume that the list is a numpy array
            old_ticks = old_ticks.tolist()
        new_xs = []
        for x_val in pos_list:
            if x_val in old_ticks:
                new_x_val = new_ticks[old_ticks.index(x_val)]
            else:
                # find smaller and bigger indices
                idx_next = np.searchsorted(old_ticks, x_val, side="left")
                idx_prev = idx_next - 1
                old_min = old_ticks[idx_prev]
                old_max = old_ticks[idx_next]
                new_min = new_ticks[idx_prev]
                new_max = new_ticks[idx_next]
                new_x_val = ((x_val - old_min) / (old_max - old_min)) * (
                    new_max - new_min
                ) + new_min
            new_xs.append(new_x_val)
        return new_xs

    dendro_info = adata.uns[dendrogram_key]["dendrogram_info"]
    leaves = dendro_info["ivl"]
    icoord = np.array(dendro_info["icoord"])
    dcoord = np.array(dendro_info["dcoord"])

    orig_ticks = np.arange(5, len(leaves) * 10 + 5, 10).astype(float)
    # check that ticks has the same length as orig_ticks
    if ticks is not None and len(orig_ticks) != len(ticks):
        logg.warning(
            "ticks argument does not have the same size as orig_ticks. "
            "The argument will be ignored"
        )
        ticks = None

    for xs, ys in zip(icoord, dcoord, strict=True):
        if ticks is not None:
            xs = translate_pos(xs, ticks, orig_ticks)  # noqa: PLW2901
        if orientation in ["right", "left"]:
            xs, ys = ys, xs  # noqa: PLW2901
        dendro_ax.plot(xs, ys, color="#555555")

    dendro_ax.tick_params(bottom=False, top=False, left=False, right=False)
    ticks = ticks if ticks is not None else orig_ticks
    if orientation in ["right", "left"]:
        dendro_ax.set_yticks(ticks)
        dendro_ax.set_yticklabels(leaves, fontsize="small", rotation=0)
        dendro_ax.tick_params(labelbottom=False, labeltop=False)
        if orientation == "left":
            xmin, xmax = dendro_ax.get_xlim()
            dendro_ax.set_xlim(xmax, xmin)
            dendro_ax.tick_params(labelleft=False, labelright=True)
    else:
        dendro_ax.set_xticks(ticks)
        dendro_ax.set_xticklabels(leaves, fontsize="small", rotation=90)
        dendro_ax.tick_params(labelleft=False, labelright=False)
        if orientation == "bottom":
            ymin, ymax = dendro_ax.get_ylim()
            dendro_ax.set_ylim(ymax, ymin)
            dendro_ax.tick_params(labeltop=True, labelbottom=False)

    if remove_labels:
        dendro_ax.tick_params(
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )

    dendro_ax.grid(visible=False)

    dendro_ax.spines["right"].set_visible(False)
    dendro_ax.spines["top"].set_visible(False)
    dendro_ax.spines["left"].set_visible(False)
    dendro_ax.spines["bottom"].set_visible(False)


def _plot_categories_as_colorblocks(
    groupby_ax: Axes,
    obs_tidy: pd.DataFrame,
    colors=None,
    orientation: Literal["top", "bottom", "left", "right"] = "left",
    cmap_name: str = "tab20",
):
    """Plot categories as colored blocks.

    If orientation is 'left', the categories are plotted vertically,
    otherwise they are plotted horizontally.

    Parameters
    ----------
    groupby_ax
    obs_tidy
    colors
        Sequence of valid color names to use for each category.
    orientation
    cmap_name
        Name of colormap to use, in case colors is None

    Returns
    -------
    ticks position, labels, colormap

    """
    groupby = obs_tidy.index.name
    from matplotlib.colors import BoundaryNorm, ListedColormap

    if colors is None:
        groupby_cmap = colormaps.get_cmap(cmap_name)
    else:
        groupby_cmap = ListedColormap(colors, groupby + "_cmap")
    norm = BoundaryNorm(np.arange(groupby_cmap.N + 1) - 0.5, groupby_cmap.N)

    # determine groupby label positions such that they appear
    # centered next/below to the color code rectangle assigned to the category
    value_sum = 0
    ticks = []  # list of centered position of the labels
    labels = []
    label2code = {}  # dictionary of numerical values asigned to each label
    for code, (label, value) in enumerate(
        obs_tidy.index.value_counts(sort=False).items()
    ):
        ticks.append(value_sum + (value / 2))
        labels.append(label)
        value_sum += value
        label2code[label] = code

    groupby_ax.grid(visible=False)

    if orientation == "left":
        groupby_ax.imshow(
            np.array([[label2code[lab] for lab in obs_tidy.index]]).T,
            aspect="auto",
            cmap=groupby_cmap,
            norm=norm,
        )
        if len(labels) > 1:
            groupby_ax.set_yticks(ticks)
            groupby_ax.set_yticklabels(labels)

        # remove y ticks
        groupby_ax.tick_params(axis="y", left=False, labelsize="small")
        # remove x ticks and labels
        groupby_ax.tick_params(axis="x", bottom=False, labelbottom=False)

        # remove surrounding lines
        groupby_ax.spines["right"].set_visible(False)
        groupby_ax.spines["top"].set_visible(False)
        groupby_ax.spines["left"].set_visible(False)
        groupby_ax.spines["bottom"].set_visible(False)

        groupby_ax.set_ylabel(groupby)
    else:
        groupby_ax.imshow(
            np.array([[label2code[lab] for lab in obs_tidy.index]]),
            aspect="auto",
            cmap=groupby_cmap,
            norm=norm,
        )
        if len(labels) > 1:
            groupby_ax.set_xticks(ticks)
            # if the labels are small do not rotate them
            rotation = 0 if max(len(str(x)) for x in labels) < 3 else 90
            groupby_ax.set_xticklabels(labels, rotation=rotation)

        # remove x ticks
        groupby_ax.tick_params(axis="x", bottom=False, labelsize="small")
        # remove y ticks and labels
        groupby_ax.tick_params(axis="y", left=False, labelleft=False)

        # remove surrounding lines
        groupby_ax.spines["right"].set_visible(False)
        groupby_ax.spines["top"].set_visible(False)
        groupby_ax.spines["left"].set_visible(False)
        groupby_ax.spines["bottom"].set_visible(False)

        groupby_ax.set_xlabel(groupby)

    return label2code, ticks, labels, groupby_cmap, norm


def _plot_colorbar(mappable, fig, subplot_spec, max_cbar_height: float = 4.0):
    """Plot a vertical color bar based on mappable.

    The height of the colorbar is min(figure-height, max_cmap_height).

    Parameters
    ----------
    mappable
        The image to which the colorbar applies.
    fig
        The figure object
    subplot_spec
        The gridspec subplot. Eg. axs[1,2]
    max_cbar_height
        The maximum colorbar height

    Returns
    -------
    color bar ax

    """
    width, height = fig.get_size_inches()
    if height > max_cbar_height:
        # to make the colorbar shorter, the
        # ax is split and the lower portion is used.
        axs2 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=subplot_spec,
            height_ratios=[height - max_cbar_height, max_cbar_height],
        )
        heatmap_cbar_ax = fig.add_subplot(axs2[1])
    else:
        heatmap_cbar_ax = fig.add_subplot(subplot_spec)
    plt.colorbar(mappable, cax=heatmap_cbar_ax)
    return heatmap_cbar_ax


def _check_var_names_type(
    var_names: _VarNames | Mapping[str, _VarNames],
    var_group_labels: Sequence[str] | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
) -> tuple[Sequence[str], VarGroups | None]:
    """If var_names is a dict, set the `var_group_labels` and `var_group_positions`.

    Returns
    -------
    var_names, var_groups

    """
    from ._baseplot_class import _var_groups

    if isinstance(var_names, Mapping):
        return _var_groups(var_names)

    if isinstance(var_names, str):
        var_names = [var_names]
    return var_names, VarGroups.validate(var_group_labels, var_group_positions)
