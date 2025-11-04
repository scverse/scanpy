from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping, Sequence
from itertools import cycle, islice
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
from cycler import Cycler, cycler
from matplotlib import axes, colormaps, gridspec, rcParams, ticker
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like
from matplotlib.figure import SubplotParams
from matplotlib.patches import Circle

from .. import logging as logg
from .._compat import old_positionals, warn
from .._settings import settings
from .._utils import NeighborsView, _empty
from . import palettes

if TYPE_CHECKING:
    from collections.abc import Collection

    from anndata import AnnData
    from matplotlib.colors import Colormap
    from matplotlib.figure import Figure
    from matplotlib.typing import MarkerType
    from numpy.typing import ArrayLike
    from PIL.Image import Image

    from .._utils import Empty

__all__ = [
    "ColorLike",
    "DensityNorm",
    "VBound",
    "_AxesSubplot",
    "_FontSize",
    "_FontWeight",
    "_LegendLoc",
    "_deprecated_scale",
    "_dk",
    "add_colors_for_categorical_sample_annotation",
    "check_colornorm",
    "check_projection",
    "circles",
    "default_palette",
    "fix_kwds",
    "make_grid_spec",
    "matrix",
    "plot_arrows",
    "plot_edges",
    "savefig_or_show",
    "scatter_base",
    "scatter_group",
    "set_colors_for_categorical_obs",
    "set_default_colors_for_categorical_obs",
    "setup_axes",
    "ticks_formatter",
    "timeseries_as_heatmap",
    "timeseries_subplot",
    "validate_palette",
]

# TODO: more
type DensityNorm = Literal["area", "count", "width"]

# These are needed by _wraps_plot_scatter
type VBound = str | float | Callable[[Sequence[float]], float]
type _FontWeight = Literal[
    "light", "normal", "medium", "semibold", "bold", "heavy", "black"
]
type _FontSize = Literal[
    "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"
]
type _LegendLoc = Literal[
    "none",
    "right margin",
    "on data",
    "on data export",
    "best",
    "upper right",
    "upper left",
    "lower left",
    "lower right",
    "right",
    "center left",
    "center right",
    "lower center",
    "upper center",
    "center",
]
type ColorLike = str | tuple[float, ...]


class _AxesSubplot(Axes, axes.SubplotBase):
    """Intersection between Axes and SubplotBase: Has methods of both."""


# -------------------------------------------------------------------------------
# Simple plotting functions
# -------------------------------------------------------------------------------


@old_positionals(
    "xlabel",
    "ylabel",
    "xticks",
    "yticks",
    "title",
    "colorbar_shrink",
    "color_map",
    "show",
    "save",
    "ax",
)
def matrix(  # noqa: PLR0913
    matrix: ArrayLike | Image,
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    xticks: Collection[str] | None = None,
    yticks: Collection[str] | None = None,
    title: str | None = None,
    colorbar_shrink: float = 0.5,
    color_map: str | Colormap | None = None,
    show: bool | None = None,
    ax: Axes | None = None,
    # deprecated
    save: bool | str | None = None,
) -> None:
    """Plot a matrix."""
    if ax is None:
        ax = plt.gca()
    img = ax.imshow(matrix, cmap=color_map)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if xticks is not None:
        ax.set_xticks(range(len(xticks)), xticks, rotation="vertical")
    if yticks is not None:
        ax.set_yticks(range(len(yticks)), yticks)
    plt.colorbar(
        img, shrink=colorbar_shrink, ax=ax
    )  # need a figure instance for colorbar
    savefig_or_show("matrix", show=show, save=save)


def timeseries(X, **kwargs):  # noqa: N803
    """Plot X. See timeseries_subplot."""
    plt.figure(
        figsize=tuple(2 * s for s in rcParams["figure.figsize"]),
        subplotpars=SubplotParams(left=0.12, right=0.98, bottom=0.13),
    )
    timeseries_subplot(X, **kwargs)


def timeseries_subplot(  # noqa: PLR0912, PLR0913
    X: np.ndarray,  # noqa: N803
    *,
    time=None,
    color=None,
    var_names=(),
    highlights_x=(),
    xlabel="",
    ylabel="gene expression",
    yticks=None,
    xlim=None,
    legend=True,
    palette: Sequence[str] | Cycler | None = None,
    color_map="viridis",
    ax: Axes | None = None,
    marker: str | Sequence[str] = ".",
):
    """Plot X.

    Parameters
    ----------
    X
        Call this with:
        X with one column, color categorical.
        X with one column, color continuous.
        X with n columns, color is of length n.

    """
    use_color_map = color is not None and isinstance(color[0], float | np.floating)
    palette = default_palette(palette)
    x_range = np.arange(X.shape[0]) if time is None else time
    if X.ndim == 1:
        X = X[:, None]  # noqa: N806
    if X.shape[1] > 1:
        colors = islice(cycle(palette.by_key()["color"]), X.shape[1])
        subsets = [(x_range, X[:, i]) for i in range(X.shape[1])]
    elif use_color_map:
        colors = [color]
        subsets = [(x_range, X[:, 0])]
    else:
        levels, _ = np.unique(color, return_inverse=True)
        colors = islice(cycle(palette.by_key()["color"]), len(levels))
        subsets = [(x_range[color == level], X[color == level, :]) for level in levels]

    if isinstance(marker, str):
        marker = [marker]
    if len(marker) != len(subsets) and len(marker) == 1:
        marker = [marker[0]] * len(subsets)
    if not (has_var_names := (len(var_names) > 0)):
        var_names = [""] * len(subsets)

    if ax is None:
        ax = plt.subplot()
    for (x, y), m, c, var_name in zip(subsets, marker, colors, var_names, strict=True):
        ax.scatter(
            x,
            y,
            marker=m,
            edgecolor="face",
            s=rcParams["lines.markersize"],
            c=c,
            label=var_name,
            rasterized=settings._vector_friendly,
            **(dict(cmap=color_map) if use_color_map else {}),
        )
    ylim = ax.get_ylim()
    for h in highlights_x:
        ax.plot([h, h], [ylim[0], ylim[1]], "--", color="black")
    ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if yticks is not None:
        ax.set_yticks(yticks)
    if has_var_names and legend:
        ax.legend(frameon=False)


def timeseries_as_heatmap(
    X: np.ndarray,  # noqa: N803
    *,
    var_names: Collection[str] = (),
    highlights_x=(),
    color_map=None,
):
    """Plot timeseries as heatmap.

    Parameters
    ----------
    X
        Data array.
    var_names
        Array of strings naming variables stored in columns of X.

    """
    if len(var_names) == 0:
        var_names = np.arange(X.shape[1])
    if var_names.ndim == 2:
        var_names = var_names[:, 0]

    X = X.T  # noqa: N806

    _, ax = plt.subplots(figsize=(1.5 * 4, 2 * 4))
    img = ax.imshow(
        np.array(X, dtype=np.float64),
        aspect="auto",
        interpolation="nearest",
        cmap=color_map,
    )
    plt.colorbar(img, shrink=0.5)
    plt.yticks(range(X.shape[0]), var_names)
    for h in highlights_x:
        plt.plot([h, h], [0, X.shape[0]], "--", color="black")
    plt.xlim([0, X.shape[1] - 1])
    plt.ylim([0, X.shape[0] - 1])


# -------------------------------------------------------------------------------
# Colors in addition to matplotlib's colors
# -------------------------------------------------------------------------------


_ADDITIONAL_COLORS = {
    "gold2": "#eec900",
    "firebrick3": "#cd2626",
    "khaki2": "#eee685",
    "slategray3": "#9fb6cd",
    "palegreen3": "#7ccd7c",
    "tomato2": "#ee5c42",
    "grey80": "#cccccc",
    "grey90": "#e5e5e5",
    "wheat4": "#8b7e66",
    "grey65": "#a6a6a6",
    "grey10": "#1a1a1a",
    "grey20": "#333333",
    "grey50": "#7f7f7f",
    "grey30": "#4d4d4d",
    "grey40": "#666666",
    "antiquewhite2": "#eedfcc",
    "grey77": "#c4c4c4",
    "snow4": "#8b8989",
    "chartreuse3": "#66cd00",
    "yellow4": "#8b8b00",
    "darkolivegreen2": "#bcee68",
    "olivedrab3": "#9acd32",
    "azure3": "#c1cdcd",
    "violetred": "#d02090",
    "mediumpurple3": "#8968cd",
    "purple4": "#551a8b",
    "seagreen4": "#2e8b57",
    "lightblue3": "#9ac0cd",
    "orchid3": "#b452cd",
    "indianred 3": "#cd5555",
    "grey60": "#999999",
    "mediumorchid1": "#e066ff",
    "plum3": "#cd96cd",
    "palevioletred3": "#cd6889",
}

# -------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------


def _savefig(writekey, dpi=None, ext=None):
    """Save current figure to file.

    The `filename` is generated as follows:

        filename = settings.figdir / f"{writekey}{settings.plot_suffix}.{settings.file_format_figs}"
    """
    if dpi is None:
        # we need this as in notebooks, the internal figures are also influenced by 'savefig.dpi' this...
        if (
            not isinstance(rcParams["savefig.dpi"], str)
            and rcParams["savefig.dpi"] < 150
        ):
            if settings._low_resolution_warning:
                logg.warning(
                    "You are using a low resolution (dpi<150) for saving figures.\n"
                    "Consider running `set_figure_params(dpi_save=...)`, which will "
                    "adjust `matplotlib.rcParams['savefig.dpi']`"
                )
                settings._low_resolution_warning = False
        else:
            dpi = rcParams["savefig.dpi"]
    settings.figdir.mkdir(parents=True, exist_ok=True)
    if ext is None:
        ext = settings.file_format_figs
    filename = settings.figdir / f"{writekey}{settings.plot_suffix}.{ext}"
    # output the following msg at warning level; it's really important for the user
    logg.warning(f"saving figure to file {filename}")
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")


def savefig_or_show(
    writekey: str,
    *,
    show: bool | None = None,
    dpi: int | None = None,
    ext: str | None = None,
    save: bool | str | None = None,
):
    if isinstance(save, str):
        # check whether `save` contains a figure extension
        if ext is None:
            for try_ext in [".svg", ".pdf", ".png"]:
                if save.endswith(try_ext):
                    ext = try_ext[1:]
                    save = save.replace(try_ext, "")
                    break
        # append it
        writekey += save
        save = True
    if do_save := settings.autosave if save is None else save:
        if save:  # `save=True | "some-str"` argument has been used
            msg = (
                "Argument `save` is deprecated and will be removed in a future version. "
                "Use `sc.pl.plot(show=False).figure.savefig()` instead."
            )
            warn(msg, FutureWarning)
        _savefig(writekey, dpi=dpi, ext=ext)
    if settings.autoshow if show is None else show:
        plt.show()
    if do_save:
        plt.close()  # clear figure


def default_palette(
    palette: str | Sequence[str] | Cycler | None = None,
) -> str | Cycler:
    if palette is None:
        return rcParams["axes.prop_cycle"]
    elif not isinstance(palette, str | Cycler):
        return cycler(color=palette)
    else:
        return palette


def validate_palette(adata: AnnData, key: str) -> None:
    """Validate and update the list of colors in `adata.uns[f'{key}_colors']`.

    Not only valid matplotlib colors are checked but also if the color name
    is a valid R color name, in which case it will be translated to a valid name
    """
    color_key = f"{key}_colors"
    raw_palette = adata.uns[color_key]
    try:
        # check if the color is a valid R color and translate it
        # to a valid hex color value
        palette = [
            color if is_color_like(color) else _ADDITIONAL_COLORS[color]
            for color in raw_palette
        ]
    except KeyError as e:
        logg.warning(
            f"The following color value found in adata.uns['{key}_colors'] "
            f"is not valid: {e.args[0]!r}. Default colors will be used instead."
        )
        set_default_colors_for_categorical_obs(adata, key)
        palette = None
    # Don’t modify if nothing changed
    if palette is None or np.array_equal(palette, adata.uns[color_key]):
        return
    adata.uns[color_key] = palette


def set_colors_for_categorical_obs(
    adata, value_to_plot, palette: str | Sequence[str] | Cycler
):
    """Set `adata.uns[f'{value_to_plot}_colors']` according to the given palette.

    Parameters
    ----------
    adata
        annData object
    value_to_plot
        name of a valid categorical observation
    palette
        Palette should be either a valid :func:`~matplotlib.pyplot.colormaps` string,
        a sequence of colors (in a format that can be understood by matplotlib,
        eg. RGB, RGBS, hex, or a cycler object with key='color'

    Returns
    -------
    None

    """
    from matplotlib.colors import to_hex

    if adata.obs[value_to_plot].dtype == bool:
        categories = (
            adata.obs[value_to_plot].astype(str).astype("category").cat.categories
        )
    else:
        categories = adata.obs[value_to_plot].cat.categories
    # check is palette is a valid matplotlib colormap
    if isinstance(palette, str) and palette in colormaps:
        # this creates a palette from a colormap. E.g. 'Accent, Dark2, tab20'
        cmap = colormaps.get_cmap(palette)
        colors_list = [to_hex(x) for x in cmap(np.linspace(0, 1, len(categories)))]
    elif isinstance(palette, Mapping):
        colors_list = [to_hex(palette[k], keep_alpha=True) for k in categories]
    else:
        # check if palette is a list and convert it to a cycler, thus
        # it doesnt matter if the list is shorter than the categories length:
        if isinstance(palette, Sequence):
            if len(palette) < len(categories):
                logg.warning(
                    "Length of palette colors is smaller than the number of "
                    f"categories (palette length: {len(palette)}, "
                    f"categories length: {len(categories)}. "
                    "Some categories will have the same color."
                )
            try:  # check that colors are valid
                _color_list = [
                    color if is_color_like(color) else _ADDITIONAL_COLORS[color]
                    for color in palette
                ]
            except KeyError as e:
                msg = (
                    "The following color value of the given palette "
                    f"is not valid: {e.args[0]!r}"
                )
                raise ValueError(msg) from None

            palette = cycler(color=_color_list)
        if not isinstance(palette, Cycler):
            msg = (
                "Please check that the value of 'palette' is a valid "
                "matplotlib colormap string (eg. Set2), a  list of color names "
                "or a cycler with a 'color' key."
            )
            raise ValueError(msg)
        if "color" not in palette.keys:
            msg = "Please set the palette key 'color'."
            raise ValueError(msg)

        cc = palette()
        colors_list = [to_hex(next(cc)["color"]) for x in range(len(categories))]

    adata.uns[f"{value_to_plot}_colors"] = colors_list


def set_default_colors_for_categorical_obs(adata, value_to_plot):
    """Set `adata.uns[f'{value_to_plot}_colors']` using default color palettes.

    Parameters
    ----------
    adata
        AnnData object
    value_to_plot
        Name of a valid categorical observation

    Returns
    -------
    None

    """
    if adata.obs[value_to_plot].dtype == bool:
        categories = (
            adata.obs[value_to_plot].astype(str).astype("category").cat.categories
        )
    else:
        categories = adata.obs[value_to_plot].cat.categories

    length = len(categories)

    # check if default matplotlib palette has enough colors
    if len(rcParams["axes.prop_cycle"].by_key()["color"]) >= length:
        cc = rcParams["axes.prop_cycle"]()
        palette = [next(cc)["color"] for _ in range(length)]

    elif length <= 20:
        palette = palettes.default_20
    elif length <= 28:
        palette = palettes.default_28
    elif length <= len(palettes.default_102):  # 103 colors
        palette = palettes.default_102
    else:
        palette = ["grey" for _ in range(length)]
        logg.info(
            f"the obs value {value_to_plot!r} has more than 103 categories. Uniform "
            "'grey' color will be used for all categories."
        )

    set_colors_for_categorical_obs(adata, value_to_plot, palette[:length])


def add_colors_for_categorical_sample_annotation(
    adata, key, *, palette=None, force_update_colors=False
):
    color_key = f"{key}_colors"
    colors_needed = len(adata.obs[key].cat.categories)
    if palette and force_update_colors:
        set_colors_for_categorical_obs(adata, key, palette)
    elif color_key in adata.uns and len(adata.uns[color_key]) <= colors_needed:
        validate_palette(adata, key)
    else:
        set_default_colors_for_categorical_obs(adata, key)


def plot_edges(axs, adata, basis, edges_width, edges_color, *, neighbors_key=None):
    import networkx as nx

    if not isinstance(axs, Sequence):
        axs = [axs]

    if neighbors_key is None:
        neighbors_key = "neighbors"
    if neighbors_key not in adata.uns:
        msg = "`edges=True` requires `pp.neighbors` to be run before."
        raise ValueError(msg)
    neighbors = NeighborsView(adata, neighbors_key)
    g = nx.Graph(neighbors["connectivities"])
    basis_key = _get_basis(adata, basis)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ax in axs:
            edge_collection = nx.draw_networkx_edges(
                g,
                adata.obsm[basis_key],
                ax=ax,
                width=edges_width,
                edge_color=edges_color,
            )
            edge_collection.set_zorder(-2)
            edge_collection.set_rasterized(settings._vector_friendly)


def plot_arrows(axs, adata, basis, arrows_kwds=None):
    if not isinstance(axs, Sequence):
        axs = [axs]
    v_prefix = next(
        (p for p in ["velocity", "Delta"] if f"{p}_{basis}" in adata.obsm), None
    )
    if v_prefix is None:
        msg = (
            "`arrows=True` requires "
            f"`'velocity_{basis}'` from scvelo or "
            f"`'Delta_{basis}'` from velocyto."
        )
        raise ValueError(msg)
    if v_prefix == "velocity":
        logg.warning(
            "The module `scvelo` has improved plotting facilities. "
            "Prefer using `scv.pl.velocity_embedding` to `arrows=True`."
        )

    basis_key = _get_basis(adata, basis)
    x = adata.obsm[basis_key]
    v = adata.obsm[f"{v_prefix}_{basis}"]
    for ax in axs:
        quiver_kwds = arrows_kwds if arrows_kwds is not None else {}
        ax.quiver(
            x[:, 0],
            x[:, 1],
            v[:, 0],
            v[:, 1],
            **quiver_kwds,
            rasterized=settings._vector_friendly,
        )


def scatter_group(
    ax: Axes,
    key: str,
    cat_code: int,
    adata: AnnData,
    y: np.ndarray,
    *,
    projection: Literal["2d", "3d"] = "2d",
    size: int = 3,
    alpha: float | None = None,
    marker: MarkerType = ".",
):
    """Scatter of group using representation of data Y."""
    mask_obs = adata.obs[key].cat.categories[cat_code] == adata.obs[key].values
    color = adata.uns[f"{key}_colors"][cat_code]
    if not isinstance(color[0], str):
        from matplotlib.colors import rgb2hex

        color = rgb2hex(adata.uns[f"{key}_colors"][cat_code])
    if not is_color_like(color):
        msg = f"{color!r} is not a valid matplotlib color."
        raise ValueError(msg)
    data = [y[mask_obs, 0], y[mask_obs, 1]]
    if projection == "3d":
        data.append(y[mask_obs, 2])
    ax.scatter(
        *data,
        marker=marker,
        alpha=alpha,
        c=color,
        edgecolors="none",
        s=size,
        label=adata.obs[key].cat.categories[cat_code],
        rasterized=settings._vector_friendly,
    )
    return mask_obs


def setup_axes(  # noqa: PLR0912
    ax: Axes | Sequence[Axes] | None = None,
    *,
    panels="blue",
    colorbars=(False,),
    right_margin=None,
    left_margin=None,
    projection: Literal["2d", "3d"] = "2d",
    show_ticks=False,
):
    """Grid of axes for plotting, legends and colorbars."""
    check_projection(projection)
    if left_margin is not None:
        msg = "We currently don’t support `left_margin`."
        raise NotImplementedError(msg)
    if np.any(colorbars) and right_margin is None:
        right_margin = 1 - rcParams["figure.subplot.right"] + 0.21  # 0.25
    elif right_margin is None:
        right_margin = 1 - rcParams["figure.subplot.right"] + 0.06  # 0.10
    # make a list of right margins for each panel
    if not isinstance(right_margin, list):
        right_margin_list = [right_margin for i in range(len(panels))]
    else:
        right_margin_list = right_margin

    # make a figure with len(panels) panels in a row side by side
    top_offset = 1 - rcParams["figure.subplot.top"]
    bottom_offset = 0.15 if show_ticks else 0.08
    left_offset = 1 if show_ticks else 0.3  # in units of base_height
    base_height = rcParams["figure.figsize"][1]
    height = base_height
    base_width = rcParams["figure.figsize"][0]
    if show_ticks:
        base_width *= 1.1

    draw_region_width = (
        base_width - left_offset - top_offset - 0.5
    )  # this is kept constant throughout

    right_margin_factor = sum([1 + right_margin for right_margin in right_margin_list])
    width_without_offsets = (
        right_margin_factor * draw_region_width
    )  # this is the total width that keeps draw_region_width

    right_offset = (len(panels) - 1) * left_offset
    figure_width = width_without_offsets + left_offset + right_offset
    draw_region_width_frac = draw_region_width / figure_width
    left_offset_frac = left_offset / figure_width
    right_offset_frac = (  # noqa: F841  # TODO Does this need fixing?
        1 - (len(panels) - 1) * left_offset_frac
    )

    if ax is None:
        plt.figure(
            figsize=(figure_width, height),
            subplotpars=SubplotParams(left=0, right=1, bottom=bottom_offset),
        )
    left_positions = [left_offset_frac, left_offset_frac + draw_region_width_frac]
    for i in range(1, len(panels)):
        right_margin = right_margin_list[i - 1]
        left_positions.append(
            left_positions[-1] + right_margin * draw_region_width_frac
        )
        left_positions.append(left_positions[-1] + draw_region_width_frac)
    panel_pos = [[bottom_offset], [1 - top_offset], left_positions]

    axs = []
    if ax is None:
        for icolor, _color in enumerate(panels):
            left = panel_pos[2][2 * icolor]
            bottom = panel_pos[0][0]
            width = draw_region_width / figure_width
            height = panel_pos[1][0] - bottom
            if projection == "2d":
                ax = plt.axes([left, bottom, width, height])
            elif projection == "3d":
                ax = plt.axes([left, bottom, width, height], projection="3d")
            axs.append(ax)
    else:
        axs = ax if isinstance(ax, Sequence) else [ax]

    return axs, panel_pos, draw_region_width, figure_width


def scatter_base(  # noqa: PLR0912, PLR0913, PLR0915
    y: np.ndarray,
    /,
    *,
    colors: str | Sequence[ColorLike | np.ndarray] = "blue",
    sort_order=True,
    alpha=None,
    highlights=(),
    right_margin=None,
    left_margin=None,
    projection: Literal["2d", "3d"] = "2d",
    title=None,
    component_name="DC",
    component_indexnames=(1, 2, 3),
    axis_labels=None,
    colorbars=(False,),
    sizes=(1,),
    markers=".",
    color_map="viridis",
    show_ticks=True,
    ax=None,
) -> Axes | list[Axes]:
    """Plot scatter plot of data.

    Parameters
    ----------
    y
        Data array.
    projection

    Returns
    -------
    Depending on whether supplying a single array or a list of arrays,
    return a single axis or a list of axes.

    """
    if isinstance(highlights, Mapping):
        highlights_indices = sorted(highlights)
        highlights_labels = [highlights[i] for i in highlights_indices]
    else:
        highlights_indices = map(int, highlights)
        highlights_labels = map(str, highlights)
    # if we have a single array, transform it into a list with a single array
    if isinstance(colors, str):
        colors = [colors]
    if isinstance(markers, str):
        markers = [markers]
    if len(sizes) != len(colors) and len(sizes) == 1:
        sizes = [sizes[0] for _ in range(len(colors))]
    if len(markers) != len(colors) and len(markers) == 1:
        markers = [markers[0] for _ in range(len(colors))]
    axs, panel_pos, draw_region_width, _figure_width = setup_axes(
        ax,
        panels=colors,
        colorbars=colorbars,
        projection=projection,
        right_margin=right_margin,
        left_margin=left_margin,
        show_ticks=show_ticks,
    )
    for icolor, color_spec in enumerate(colors):
        ax = axs[icolor]
        marker = markers[icolor]
        bottom = panel_pos[0][0]
        height = panel_pos[1][0] - bottom
        y_sort = y
        if not is_color_like(color_spec) and sort_order:
            sort = np.argsort(color_spec)
            color = color_spec[sort]
            y_sort = y[sort]
        else:
            color = color_spec
        if projection == "2d":
            data = y_sort[:, 0], y_sort[:, 1]
        elif projection == "3d":
            data = y_sort[:, 0], y_sort[:, 1], y_sort[:, 2]
        else:
            msg = f"Unknown projection {projection!r} not in '2d', '3d'"
            raise ValueError(msg)
        if not isinstance(color, str) or color != "white":
            sct = ax.scatter(
                *data,
                marker=marker,
                c=color,
                alpha=alpha,
                edgecolors="none",  # 'face',
                s=sizes[icolor],
                cmap=color_map,
                rasterized=settings._vector_friendly,
            )
        if colorbars[icolor]:
            width = 0.006 * draw_region_width / len(colors)
            left = (
                panel_pos[2][2 * icolor + 1]
                + (1.2 if projection == "3d" else 0.2) * width
            )
            rectangle = [left, bottom, width, height]
            fig = plt.gcf()
            ax_cb = fig.add_axes(rectangle)
            _ = plt.colorbar(
                sct, format=ticker.FuncFormatter(ticks_formatter), cax=ax_cb
            )
        # set the title
        if title is not None:
            ax.set_title(title[icolor])
        # output highlighted data points
        for ihighlight, highlight_text in zip(
            highlights_indices, highlights_labels, strict=True
        ):
            data = [y[ihighlight, 0]], [y[ihighlight, 1]]
            if "3d" in projection:
                data = [y[ihighlight, 0]], [y[ihighlight, 1]], [y[ihighlight, 2]]
            ax.scatter(
                *data,
                c="black",
                facecolors="black",
                edgecolors="black",
                marker="x",
                s=10,
                zorder=20,
            )
            # the following is a Python 2 compatibility hack
            ax.text(
                *([d[0] for d in data] + [highlight_text]),
                zorder=20,
                fontsize=10,
                color="black",
            )
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
            if "3d" in projection:
                ax.set_zticks([])
    # set default axis_labels
    if axis_labels is None:
        axis_labels = [
            [component_name + str(i) for i in component_indexnames]
            for _ in range(len(axs))
        ]
    else:
        axis_labels = [axis_labels for _ in range(len(axs))]
    for iax, ax in enumerate(axs):
        ax.set_xlabel(axis_labels[iax][0])
        ax.set_ylabel(axis_labels[iax][1])
        if "3d" in projection:
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[iax][2], labelpad=-7)
    for ax in axs:
        # scale limits to match data
        ax.autoscale_view()
    return axs


def ticks_formatter(x, pos) -> str:
    return f"{x:.3f}".rstrip("0").rstrip(".")


def check_projection(projection):
    """Validate projection argument."""
    if projection not in {"2d", "3d"}:
        msg = f"Projection must be '2d' or '3d', was '{projection}'."
        raise ValueError(msg)


def circles(
    x, y, *, s, ax, marker=None, c="b", vmin=None, vmax=None, scale_factor=1.0, **kwargs
) -> PatchCollection:
    """Make a scatter plot of circles.

    Similar to pl.scatter, but the size of circles are in data scale.

    Taken from here: <https://gist.github.com/syrte/592a062c562cd2a98a83>

    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    pl.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (https://opensource.org/license/bsd-3-clause/)

    """
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.
    if scale_factor != 1.0:
        x = x * scale_factor
        y = y * scale_factor
    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_) for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if isinstance(c, np.ndarray) and np.issubdtype(c.dtype, np.number):
        collection.set_array(np.ma.masked_invalid(c))
        collection.set_clim(vmin, vmax)
    else:
        collection.set_facecolor(c)

    ax.add_collection(collection)

    return collection


def make_grid_spec(
    ax_or_figsize: tuple[int, int] | _AxesSubplot,
    *,
    nrows: int,
    ncols: int,
    wspace: float | None = None,
    hspace: float | None = None,
    width_ratios: Sequence[float] | None = None,
    height_ratios: Sequence[float] | None = None,
) -> tuple[Figure, gridspec.GridSpecBase]:
    kw = dict(
        wspace=wspace,
        hspace=hspace,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )
    if isinstance(ax_or_figsize, tuple):
        fig = plt.figure(figsize=ax_or_figsize)
        return fig, gridspec.GridSpec(nrows, ncols, **kw)
    else:
        ax = ax_or_figsize
        ax.axis("off")
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        return ax.figure, ax.get_subplotspec().subgridspec(nrows, ncols, **kw)


def fix_kwds(kwds_dict, **kwargs):
    """Merge the parameters into a single consolidated dictionary.

    Given a dictionary of plot parameters (`kwds_dict`) and a dict of `kwds`,
    this function prevents argument duplication errors.

    If `kwds_dict` an kwargs have the same key, only the value in `kwds_dict` is kept.

    Parameters
    ----------
    kwds_dict
        kwds dictionary
    kwargs

    Returns
    -------
    `kwds_dict` merged with `kwargs`

    Examples
    --------
    >>> def _example(**kwds):
    ...     return fix_kwds(kwds, key1="value1", key2="value2")
    >>> _example(key1="value10", key3="value3")
    {'key1': 'value10', 'key2': 'value2', 'key3': 'value3'}

    """
    kwargs.update(kwds_dict)

    return kwargs


def _get_basis(adata: AnnData, basis: str):
    if basis in adata.obsm:
        basis_key = basis

    elif f"X_{basis}" in adata.obsm:
        basis_key = f"X_{basis}"

    return basis_key


def check_colornorm(vmin=None, vmax=None, vcenter=None, norm=None):
    from matplotlib.colors import Normalize

    try:
        from matplotlib.colors import TwoSlopeNorm as DivNorm
    except ImportError:
        # matplotlib<3.2
        from matplotlib.colors import DivergingNorm as DivNorm

    if norm is not None:
        if (vmin is not None) or (vmax is not None) or (vcenter is not None):
            msg = "Passing both norm and vmin/vmax/vcenter is not allowed."
            raise ValueError(msg)
    elif vcenter is not None:
        norm = DivNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    return norm


@overload
def _deprecated_scale(
    density_norm: DensityNorm,
    scale: DensityNorm | Empty,
    *,
    default: DensityNorm,
) -> DensityNorm: ...


@overload
def _deprecated_scale(
    density_norm: DensityNorm | Empty,
    scale: DensityNorm | Empty,
    *,
    default: DensityNorm | Empty = _empty,
) -> DensityNorm | Empty: ...


def _deprecated_scale(
    density_norm: DensityNorm | Empty,
    scale: DensityNorm | Empty,
    *,
    default: DensityNorm | Empty = _empty,
) -> DensityNorm | Empty:
    if scale is _empty:
        return density_norm
    if density_norm != default:
        msg = "can’t specify both `scale` and `density_norm`"
        raise ValueError(msg)
    msg = "`scale` is deprecated, use `density_norm` instead"
    warn(msg, FutureWarning)
    return scale


def _dk(dendrogram: bool | str | None) -> str | None:  # noqa: FBT001
    """Convert the `dendrogram` parameter to a `dendrogram_key` parameter."""
    return None if isinstance(dendrogram, bool) else dendrogram
