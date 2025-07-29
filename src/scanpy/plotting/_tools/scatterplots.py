from __future__ import annotations

import inspect
from collections.abc import Mapping, Sequence  # noqa: TC003
from copy import copy
from functools import partial
from itertools import combinations, product
from numbers import Integral
from typing import (
    TYPE_CHECKING,
    Any,  # noqa: TC003
    Literal,  # noqa: TC003
)

import numpy as np
import pandas as pd
from anndata import AnnData  # noqa: TC002
from cycler import Cycler  # noqa: TC002
from matplotlib import colormaps, colors, patheffects, rcParams
from matplotlib import pyplot as plt
from matplotlib.axes import Axes  # noqa: TC002
from matplotlib.colors import (
    Colormap,  # noqa: TC002
    Normalize,
)
from matplotlib.figure import Figure  # noqa: TC002
from numpy.typing import NDArray  # noqa: TC002
from packaging.version import Version

from ... import logging as logg
from ..._compat import deprecated
from ..._settings import settings
from ..._utils import (
    Empty,  # noqa: TC001
    _doc_params,
    _empty,
    sanitize_anndata,
)
from ...get import _check_mask
from ...tools._draw_graph import _Layout  # noqa: TC001
from .. import _utils
from .._docs import (
    doc_adata_color_etc,
    doc_edges_arrows,
    doc_scatter_embedding,
    doc_scatter_spatial,
    doc_show_save_ax,
)
from .._utils import (
    ColorLike,  # noqa: TC001
    VBound,  # noqa: TC001
    _FontSize,  # noqa: TC001
    _FontWeight,  # noqa: TC001
    _LegendLoc,  # noqa: TC001
    check_colornorm,
    check_projection,
    circles,
)

if TYPE_CHECKING:
    from collections.abc import Collection


@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def embedding(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    basis: str,
    *,
    color: str | Sequence[str] | None = None,
    mask_obs: NDArray[np.bool_] | str | None = None,
    gene_symbols: str | None = None,
    use_raw: bool | None = None,
    sort_order: bool = True,
    edges: bool = False,
    edges_width: float = 0.1,
    edges_color: str | Sequence[float] | Sequence[str] = "grey",
    neighbors_key: str | None = None,
    arrows: bool = False,
    arrows_kwds: Mapping[str, Any] | None = None,
    groups: str | Sequence[str] | None = None,
    components: str | Sequence[str] | None = None,
    dimensions: tuple[int, int] | Sequence[tuple[int, int]] | None = None,
    layer: str | None = None,
    projection: Literal["2d", "3d"] = "2d",
    scale_factor: float | None = None,
    color_map: Colormap | str | None = None,
    cmap: Colormap | str | None = None,
    palette: str | Sequence[str] | Cycler | None = None,
    na_color: ColorLike = "lightgray",
    na_in_legend: bool = True,
    size: float | Sequence[float] | None = None,
    frameon: bool | None = None,
    legend_fontsize: float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight = "bold",
    legend_loc: _LegendLoc | None = "right margin",
    legend_fontoutline: int | None = None,
    colorbar_loc: str | None = "right",
    vmax: VBound | Sequence[VBound] | None = None,
    vmin: VBound | Sequence[VBound] | None = None,
    vcenter: VBound | Sequence[VBound] | None = None,
    norm: Normalize | Sequence[Normalize] | None = None,
    add_outline: bool | None = False,
    outline_width: tuple[float, float] = (0.3, 0.05),
    outline_color: tuple[str, str] = ("black", "white"),
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: float | None = None,
    title: str | Sequence[str] | None = None,
    show: bool | None = None,
    save: bool | str | None = None,
    ax: Axes | None = None,
    return_fig: bool | None = None,
    marker: str | Sequence[str] = ".",
    **kwargs,
) -> Figure | Axes | list[Axes] | None:
    """Scatter plot for user specified embedding basis (e.g. umap, pca, etc).

    Parameters
    ----------
    basis
        Name of the `obsm` basis to use.
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    """
    #####################
    # Argument handling #
    #####################

    check_projection(projection)
    sanitize_anndata(adata)

    basis_values = _get_basis(adata, basis)
    dimensions = _components_to_dimensions(
        components, dimensions, projection=projection, total_dims=basis_values.shape[1]
    )
    args_3d = dict(projection="3d") if projection == "3d" else {}

    # Checking the mask format and if used together with groups
    if groups is not None and mask_obs is not None:
        msg = "Groups and mask arguments are incompatible."
        raise ValueError(msg)
    mask_obs = _check_mask(adata, mask_obs, "obs")

    # Figure out if we're using raw
    if use_raw is None:
        # check if adata.raw is set
        use_raw = layer is None and adata.raw is not None
    if use_raw and layer is not None:
        msg = (
            "Cannot use both a layer and the raw representation. "
            f"Was passed: {use_raw=!r}, {layer=!r}."
        )
        raise ValueError(msg)
    if use_raw and adata.raw is None:
        msg = (
            "`use_raw` is set to True but AnnData object does not have raw. "
            "Please check."
        )
        raise ValueError(msg)

    if isinstance(groups, str):
        groups = [groups]

    # Color map
    if color_map is not None:
        if cmap is not None:
            msg = "Cannot specify both `color_map` and `cmap`."
            raise ValueError(msg)
        else:
            cmap = color_map
    cmap = copy(colormaps.get_cmap(cmap))
    cmap.set_bad(na_color)
    # Prevents warnings during legend creation
    na_color = colors.to_hex(na_color, keep_alpha=True)

    # by default turn off edge color. Otherwise, for
    # very small sizes the edge will not reduce its size
    # (https://github.com/scverse/scanpy/issues/293)
    kwargs.setdefault("edgecolor", "none")

    # Vectorized arguments

    # turn color into a python list
    color = [color] if isinstance(color, str) or color is None else list(color)

    # turn marker into a python list
    marker = [marker] if isinstance(marker, str) else list(marker)

    if title is not None:
        # turn title into a python list if not None
        title = [title] if isinstance(title, str) else list(title)

    # turn vmax and vmin into a sequence
    if isinstance(vmax, str) or not isinstance(vmax, Sequence):
        vmax = [vmax]
    if isinstance(vmin, str) or not isinstance(vmin, Sequence):
        vmin = [vmin]
    if isinstance(vcenter, str) or not isinstance(vcenter, Sequence):
        vcenter = [vcenter]
    if isinstance(norm, Normalize) or not isinstance(norm, Sequence):
        norm = [norm]

    # Size
    if "s" in kwargs and size is None:
        size = kwargs.pop("s")
    if size is not None:
        # check if size is any type of sequence, and if so
        # set as ndarray
        if (
            size is not None
            and isinstance(size, Sequence | pd.Series | np.ndarray)
            and len(size) == adata.shape[0]
        ):
            size = np.array(size, dtype=float)
    else:
        size = 120000 / adata.shape[0]

    ##########
    # Layout #
    ##########
    # Most of the code is for the case when multiple plots are required

    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams["figure.figsize"][0] + 0.02

    if components is not None:
        color, dimensions = list(zip(*product(color, dimensions), strict=True))

    color, dimensions, marker = _broadcast_args(color, dimensions, marker)

    # 'color' is a list of names that want to be plotted.
    # Eg. ['Gene1', 'louvain', 'Gene2'].
    # component_list is a list of components [[0,1], [1,2]]
    if (
        not isinstance(color, str) and isinstance(color, Sequence) and len(color) > 1
    ) or len(dimensions) > 1:
        if ax is not None:
            msg = (
                "Cannot specify `ax` when plotting multiple panels "
                "(each for a given value of 'color')."
            )
            raise ValueError(msg)

        # each plot needs to be its own panel
        fig, grid = _panel_grid(hspace, wspace, ncols, len(color))
    else:
        grid = None
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, **args_3d)

    ############
    # Plotting #
    ############
    axs = []

    # use itertools.product to make a plot for each color and for each component
    # For example if color=[gene1, gene2] and components=['1,2, '2,3'].
    # The plots are: [
    #     color=gene1, components=[1,2], color=gene1, components=[2,3],
    #     color=gene2, components = [1, 2], color=gene2, components=[2,3],
    # ]
    for count, (value_to_plot, dims) in enumerate(zip(color, dimensions, strict=True)):
        kwargs_scatter = kwargs.copy()  # is potentially mutated for each plot
        color_source_vector = _get_color_source_vector(
            adata,
            value_to_plot,
            layer=layer,
            mask_obs=mask_obs,
            use_raw=use_raw,
            gene_symbols=gene_symbols,
            groups=groups,
        )
        color_vector, color_type = _color_vector(
            adata,
            value_to_plot,
            values=color_source_vector,
            palette=palette,
            na_color=na_color,
        )

        # Order points
        order = slice(None)
        if sort_order and value_to_plot is not None and color_type == "cont":
            # Higher values plotted on top, null values on bottom
            order = np.argsort(-color_vector, kind="stable")[::-1]
        elif sort_order and color_type == "cat":
            # Null points go on bottom
            order = np.argsort(~pd.isnull(color_source_vector), kind="stable")
        # Set orders
        if isinstance(size, np.ndarray):
            size = np.array(size)[order]
        color_source_vector = color_source_vector[order]
        color_vector = color_vector[order]
        coords = basis_values[:, dims][order, :]

        # if plotting multiple panels, get the ax from the grid spec
        # else use the ax value (either user given or created previously)
        if grid:
            ax = plt.subplot(grid[count], **args_3d)
            axs.append(ax)
        if not (settings._frameon if frameon is None else frameon):
            ax.axis("off")
        if title is None:
            if value_to_plot is not None:
                ax.set_title(value_to_plot)
            else:
                ax.set_title("")
        else:
            try:
                ax.set_title(title[count])
            except IndexError:
                logg.warning(
                    "The title list is shorter than the number of panels. "
                    "Using 'color' value instead for some plots."
                )
                ax.set_title(value_to_plot)

        if color_type == "cont":
            vmin_float, vmax_float, vcenter_float, norm_obj = _get_vboundnorm(
                vmin, vmax, vcenter, norm=norm, index=count, colors=color_vector
            )
            kwargs_scatter["norm"] = check_colornorm(
                vmin_float,
                vmax_float,
                vcenter_float,
                norm_obj,
            )
            kwargs_scatter["cmap"] = cmap

        # make the scatter plot
        if projection == "3d":
            cax = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                coords[:, 2],
                c=color_vector,
                rasterized=settings._vector_friendly,
                marker=marker[count],
                **kwargs_scatter,
            )
        else:
            scatter = (
                partial(ax.scatter, s=size, plotnonfinite=True)
                if scale_factor is None
                else partial(
                    circles, s=size, ax=ax, scale_factor=scale_factor
                )  # size in circles is radius
            )

            if add_outline:
                # the default outline is a black edge followed by a
                # thin white edged added around connected clusters.
                # To add an outline
                # three overlapping scatter plots are drawn:
                # First black dots with slightly larger size,
                # then, white dots a bit smaller, but still larger
                # than the final dots. Then the final dots are drawn
                # with some transparency.

                bg_width, gap_width = outline_width
                point = np.sqrt(size)
                gap_size = (point + (point * gap_width) * 2) ** 2
                bg_size = (np.sqrt(gap_size) + (point * bg_width) * 2) ** 2
                # the default black and white colors can be changes using
                # the contour_config parameter
                bg_color, gap_color = outline_color

                # remove edge from kwargs if present
                # because edge needs to be set to None
                kwargs_scatter["edgecolor"] = "none"
                # For points, if user did not set alpha, set alpha to 0.7
                kwargs_scatter.setdefault("alpha", 0.7)

                # remove alpha and color mapping for outline
                kwargs_outline = {
                    k: v
                    for k, v in kwargs.items()
                    if k not in {"alpha", "cmap", "norm"}
                }

                for s, c in [(bg_size, bg_color), (gap_size, gap_color)]:
                    ax.scatter(
                        coords[:, 0],
                        coords[:, 1],
                        s=s,
                        c=c,
                        rasterized=settings._vector_friendly,
                        marker=marker[count],
                        **kwargs_outline,
                    )

            cax = scatter(
                coords[:, 0],
                coords[:, 1],
                c=color_vector,
                rasterized=settings._vector_friendly,
                marker=marker[count],
                **kwargs_scatter,
            )

        # remove y and x ticks
        ax.set_yticks([])
        ax.set_xticks([])
        if projection == "3d":
            ax.set_zticks([])

        # set default axis_labels
        name = _basis2name(basis)
        axis_labels = [name + str(d + 1) for d in dims]

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        if projection == "3d":
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[2], labelpad=-7)
        ax.autoscale_view()

        if edges:
            _utils.plot_edges(
                ax, adata, basis, edges_width, edges_color, neighbors_key=neighbors_key
            )
        if arrows:
            _utils.plot_arrows(ax, adata, basis, arrows_kwds)

        if value_to_plot is None:
            # if only dots were plotted without an associated value
            # there is not need to plot a legend or a colorbar
            continue

        if legend_fontoutline is not None:
            path_effect = [
                patheffects.withStroke(linewidth=legend_fontoutline, foreground="w")
            ]
        else:
            path_effect = None

        # Adding legends
        if color_type == "cat":
            _add_categorical_legend(
                ax,
                color_source_vector,
                palette=_get_palette(adata, value_to_plot),
                scatter_array=coords,
                legend_loc=legend_loc,
                legend_fontweight=legend_fontweight,
                legend_fontsize=legend_fontsize,
                legend_fontoutline=path_effect,
                na_color=na_color,
                na_in_legend=na_in_legend,
                multi_panel=bool(grid),
            )
        elif colorbar_loc is not None:
            plt.colorbar(
                cax, ax=ax, pad=0.01, fraction=0.08, aspect=30, location=colorbar_loc
            )

    if return_fig is True:
        return fig
    axs = axs if grid else ax
    _utils.savefig_or_show(basis, show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return axs


def _panel_grid(hspace, wspace, ncols, num_panels):
    from matplotlib import gridspec

    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    # each panel will have the size of rcParams['figure.figsize']
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams["figure.figsize"][0] * (1 + wspace),
            n_panels_y * rcParams["figure.figsize"][1],
        ),
    )
    left = 0.2 / n_panels_x
    bottom = 0.13 / n_panels_y
    gs = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        left=left,
        right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
        bottom=bottom,
        top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
        hspace=hspace,
        wspace=wspace,
    )
    return fig, gs


def _get_vboundnorm(
    vmin: Sequence[VBound],
    vmax: Sequence[VBound],
    vcenter: Sequence[VBound],
    *,
    norm: Sequence[Normalize],
    index: int,
    colors: Sequence[float],
) -> tuple[float | None, float | None]:
    """Evaluate the value of `vmin`, `vmax` and `vcenter`.

    Each could be a str in which case is interpreted as a percentile and should
    be specified in the form `pN` where `N` is the percentile.
    Eg. for a percentile of 85 the format would be `p85`.
    Floats are accepted as `p99.9`.

    Alternatively, `vmin`/`vmax` could be a function that is applied to
    the list of color values (`colors`). E.g.

    >>> def my_vmax(colors):
    ...     return np.percentile(colors, p=80)

    Parameters
    ----------
    index
        This index of the plot
    colors
        Values for the plot

    Returns
    -------
    (vmin, vmax, vcenter, norm) containing None or float values for
    vmin, vmax, vcenter and matplotlib.colors.Normalize  or None for norm.

    """
    out = []
    for v_name, v in [("vmin", vmin), ("vmax", vmax), ("vcenter", vcenter)]:
        if len(v) == 1:
            # this case usually happens when the user sets eg vmax=0.9, which
            # is internally converted into list of len=1, but is expected that this
            # value applies to all plots.
            v_value = v[0]
        else:
            try:
                v_value = v[index]
            except IndexError:
                logg.error(
                    f"The parameter {v_name} is not valid. If setting multiple {v_name} values,"
                    f"check that the length of the {v_name} list is equal to the number "
                    "of plots. "
                )
                v_value = None

        if v_value is not None:
            if isinstance(v_value, str) and v_value.startswith("p"):
                try:
                    float(v_value[1:])
                except ValueError:
                    logg.error(
                        f"The parameter {v_name}={v_value} for plot number {index + 1} is not valid. "
                        f"Please check the correct format for percentiles."
                    )
                # interpret value of vmin/vmax as quantile with the following syntax 'p99.9'
                v_value = np.nanpercentile(colors, q=float(v_value[1:]))
            elif callable(v_value):
                # interpret vmin/vmax as function
                v_value = v_value(colors)
                if not isinstance(v_value, float):
                    logg.error(
                        f"The return of the function given for {v_name} is not valid. "
                        "Please check that the function returns a number."
                    )
                    v_value = None
            else:
                try:
                    float(v_value)
                except ValueError:
                    logg.error(
                        f"The given {v_name}={v_value} for plot number {index + 1} is not valid. "
                        f"Please check that the value given is a valid number, a string "
                        f"starting with 'p' for percentiles or a valid function."
                    )
                    v_value = None
        out.append(v_value)
    out.append(norm[0] if len(norm) == 1 else norm[index])
    return tuple(out)


def _wraps_plot_scatter(wrapper):
    """Update the wrapper function to use the correct signature."""
    params = inspect.signature(embedding, eval_str=True).parameters.copy()
    wrapper_sig = inspect.signature(wrapper, eval_str=True)
    wrapper_params = wrapper_sig.parameters.copy()

    params.pop("basis")
    params.pop("kwargs")
    wrapper_params.pop("adata")

    params.update(wrapper_params)
    annotations = {
        k: v.annotation
        for k, v in params.items()
        if v.annotation != inspect.Parameter.empty
    }
    if wrapper_sig.return_annotation is not inspect.Signature.empty:
        annotations["return"] = wrapper_sig.return_annotation

    wrapper.__signature__ = inspect.Signature(
        list(params.values()), return_annotation=wrapper_sig.return_annotation
    )
    wrapper.__annotations__ = annotations

    return wrapper


# API


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def umap(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in UMAP basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.umap(adata)

    Colour points by discrete variable (Louvain clusters).

    .. plot::
        :context: close-figs

        sc.pl.umap(adata, color="louvain")

    Colour points by gene expression.

    .. plot::
        :context: close-figs

        sc.pl.umap(adata, color="HES4")

    Plot muliple umaps for different gene expressions.

    .. plot::
        :context: close-figs

        sc.pl.umap(adata, color=["HES4", "TNFRSF4"])

    .. currentmodule:: scanpy

    See Also
    --------
    tl.umap

    """
    return embedding(adata, "umap", **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def tsne(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in tSNE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.tsne(adata)
        sc.pl.tsne(adata, color='bulk_labels')

    .. currentmodule:: scanpy

    See Also
    --------
    tl.tsne

    """
    return embedding(adata, "tsne", **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def diffmap(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.diffmap(adata)
        sc.pl.diffmap(adata, color='bulk_labels')

    .. currentmodule:: scanpy

    See Also
    --------
    tl.diffmap

    """
    return embedding(adata, "diffmap", **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def draw_graph(
    adata: AnnData, *, layout: _Layout | None = None, **kwargs
) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in graph-drawing basis.

    Parameters
    ----------
    {adata_color_etc}
    layout
        One of the :func:`~scanpy.tl.draw_graph` layouts.
        By default, the last computed layout is used.
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.draw_graph(adata)
        sc.pl.draw_graph(adata, color=['phase', 'bulk_labels'])

    .. currentmodule:: scanpy

    See Also
    --------
    tl.draw_graph

    """
    if layout is None:
        layout = str(adata.uns["draw_graph"]["params"]["layout"])
    basis = f"draw_graph_{layout}"
    if f"X_{basis}" not in adata.obsm_keys():
        msg = f"Did not find {basis} in adata.obs. Did you compute layout {layout}?"
        raise ValueError(msg)

    return embedding(adata, basis, **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def pca(
    adata: AnnData,
    *,
    annotate_var_explained: bool = False,
    show: bool | None = None,
    return_fig: bool | None = None,
    save: bool | str | None = None,
    **kwargs,
) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in PCA coordinates.

    Use the parameter `annotate_var_explained` to annotate the explained variance.

    Parameters
    ----------
    {adata_color_etc}
    annotate_var_explained
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.pl.pca(adata)

    Colour points by discrete variable (Louvain clusters).

    .. plot::
        :context: close-figs

        sc.pl.pca(adata, color="louvain")

    Colour points by gene expression.

    .. plot::
        :context: close-figs

        sc.pl.pca(adata, color="CST3")

    .. currentmodule:: scanpy

    See Also
    --------
    pp.pca

    """
    if not annotate_var_explained:
        return embedding(
            adata, "pca", show=show, return_fig=return_fig, save=save, **kwargs
        )
    if "pca" not in adata.obsm and "X_pca" not in adata.obsm:
        msg = (
            f"Could not find entry in `obsm` for 'pca'.\n"
            f"Available keys are: {list(adata.obsm.keys())}."
        )
        raise KeyError(msg)

    label_dict = {
        f"PC{i + 1}": f"PC{i + 1} ({round(v * 100, 2)}%)"
        for i, v in enumerate(adata.uns["pca"]["variance_ratio"])
    }

    if return_fig is True:
        # edit axis labels in returned figure
        fig = embedding(adata, "pca", return_fig=return_fig, **kwargs)
        for ax in fig.axes:
            if xlabel := label_dict.get(ax.xaxis.get_label().get_text()):
                ax.set_xlabel(xlabel)
            if ylabel := label_dict.get(ax.yaxis.get_label().get_text()):
                ax.set_ylabel(ylabel)
        return fig

    # get the axs, edit the labels and apply show and save from user
    axs = embedding(adata, "pca", show=False, save=False, **kwargs)
    if isinstance(axs, list):
        for ax in axs:
            ax.set_xlabel(label_dict[ax.xaxis.get_label().get_text()])
            ax.set_ylabel(label_dict[ax.yaxis.get_label().get_text()])
    else:
        axs.set_xlabel(label_dict[axs.xaxis.get_label().get_text()])
        axs.set_ylabel(label_dict[axs.yaxis.get_label().get_text()])
    _utils.savefig_or_show("pca", show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return axs


@deprecated("Use `squidpy.pl.spatial_scatter` instead.")
@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_spatial=doc_scatter_spatial,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def spatial(  # noqa: PLR0913
    adata: AnnData,
    *,
    basis: str = "spatial",
    img: np.ndarray | None = None,
    img_key: str | None | Empty = _empty,
    library_id: str | None | Empty = _empty,
    crop_coord: tuple[int, int, int, int] | None = None,
    alpha_img: float = 1.0,
    bw: bool | None = False,
    size: float = 1.0,
    scale_factor: float | None = None,
    spot_size: float | None = None,
    na_color: ColorLike | None = None,
    show: bool | None = None,
    return_fig: bool | None = None,
    save: bool | str | None = None,
    **kwargs,
) -> Figure | Axes | list[Axes] | None:
    """Scatter plot in spatial coordinates.

    .. deprecated:: 1.11.0
       Use :func:`squidpy.pl.spatial_scatter` instead.

    This function allows overlaying data on top of images.
    Use the parameter `img_key` to see the image in the background
    And the parameter `library_id` to select the image.
    By default, `'hires'` and `'lowres'` are attempted.

    Use `crop_coord`, `alpha_img`, and `bw` to control how it is displayed.
    Use `size` to scale the size of the Visium spots plotted on top.

    As this function is designed to for imaging data, there are two key assumptions
    about how coordinates are handled:

    1. The origin (e.g `(0, 0)`) is at the top left – as is common convention
    with image data.

    2. Coordinates are in the pixel space of the source image, so an equal
    aspect ratio is assumed.

    If your anndata object has a `"spatial"` entry in `.uns`, the `img_key`
    and `library_id` parameters to find values for `img`, `scale_factor`,
    and `spot_size` arguments. Alternatively, these values be passed directly.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_spatial}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    This function behaves very similarly to other embedding plots like
    :func:`~scanpy.pl.umap`

    >>> import scanpy as sc
    >>> adata = sc.datasets.visium_sge("Targeted_Visium_Human_Glioblastoma_Pan_Cancer")
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.pl.spatial(adata, color="log1p_n_genes_by_counts")

    See Also
    --------
    :func:`scanpy.datasets.visium_sge`
        Example visium data.

    """
    # get default image params if available
    library_id, spatial_data = _check_spatial_data(adata.uns, library_id)
    img, img_key = _check_img(spatial_data, img, img_key, bw=bw)
    spot_size = _check_spot_size(spatial_data, spot_size)
    scale_factor = _check_scale_factor(
        spatial_data, img_key=img_key, scale_factor=scale_factor
    )
    crop_coord = _check_crop_coord(crop_coord, scale_factor)
    na_color = _check_na_color(na_color, img=img)

    cmap_img = "gray" if bw else None
    circle_radius = size * scale_factor * spot_size * 0.5

    axs = embedding(
        adata,
        basis=basis,
        scale_factor=scale_factor,
        size=circle_radius,
        na_color=na_color,
        show=False,
        save=False,
        **kwargs,
    )
    if not isinstance(axs, list):
        axs = [axs]
    for ax in axs:
        cur_coords = np.concatenate([ax.get_xlim(), ax.get_ylim()])
        if img is not None:
            ax.imshow(img, cmap=cmap_img, alpha=alpha_img)
        else:
            ax.set_aspect("equal")
            ax.invert_yaxis()
        if crop_coord is not None:
            ax.set_xlim(crop_coord[0], crop_coord[1])
            ax.set_ylim(crop_coord[3], crop_coord[2])
        else:
            ax.set_xlim(cur_coords[0], cur_coords[1])
            ax.set_ylim(cur_coords[3], cur_coords[2])
    _utils.savefig_or_show("show", show=show, save=save)
    if return_fig:
        return axs[0].figure
    show = settings.autoshow if show is None else show
    if show:
        return None
    return axs


# Helpers
def _components_to_dimensions(
    components: str | Collection[str] | None,
    dimensions: Collection[int] | Collection[Collection[int]] | None,
    *,
    projection: Literal["2d", "3d"] = "2d",
    total_dims: int,
) -> list[Collection[int]]:
    """Normalize components/ dimensions args for embedding plots."""
    # TODO: Deprecate components kwarg
    ndims = {"2d": 2, "3d": 3}[projection]
    if components is None and dimensions is None:
        dimensions = [tuple(i for i in range(ndims))]
    elif components is not None and dimensions is not None:
        msg = "Cannot provide both dimensions and components"
        raise ValueError(msg)

    # TODO: Consider deprecating this
    # If components is not None, parse them and set dimensions
    if components == "all":
        dimensions = list(combinations(range(total_dims), ndims))
    elif components is not None:
        if isinstance(components, str):
            components = [components]
        # Components use 1 based indexing
        dimensions = [[int(dim) - 1 for dim in c.split(",")] for c in components]

    if all(isinstance(el, Integral) for el in dimensions):
        dimensions = [dimensions]
    # if all(isinstance(el, Collection) for el in dimensions):
    for dims in dimensions:
        if len(dims) != ndims or not all(isinstance(d, Integral) for d in dims):
            raise ValueError()

    return dimensions


def _add_categorical_legend(  # noqa: PLR0913
    ax: Axes,
    color_source_vector,
    *,
    palette: dict,
    legend_loc: _LegendLoc | None,
    legend_fontweight,
    legend_fontsize,
    legend_fontoutline,
    multi_panel,
    na_color,
    na_in_legend: bool,
    scatter_array=None,
):
    """Add a legend to the passed Axes."""
    if na_in_legend and pd.isnull(color_source_vector).any():
        if "NA" in color_source_vector:
            msg = "No fallback for null labels has been defined if NA already in categories."
            raise NotImplementedError(msg)
        color_source_vector = color_source_vector.add_categories("NA").fillna("NA")
        palette = palette.copy()
        palette["NA"] = na_color
    if color_source_vector.dtype == bool:
        cats = pd.Categorical(color_source_vector.astype(str)).categories
    else:
        cats = color_source_vector.categories

    if multi_panel is True:
        # Shrink current axis by 10% to fit legend and match
        # size of plots that are not categorical
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.91, box.height])

    if legend_loc == "on data":
        # identify centroids to put labels

        all_pos = (
            pd.DataFrame(scatter_array, columns=["x", "y"])
            .groupby(color_source_vector, observed=True)
            .median()
            # Have to sort_index since if observed=True and categorical is unordered
            # the order of values in .index is undefined. Related issue:
            # https://github.com/pandas-dev/pandas/issues/25167
            .sort_index()
        )

        for label, x_pos, y_pos in all_pos.itertuples():
            ax.text(
                x_pos,
                y_pos,
                label,
                weight=legend_fontweight,
                verticalalignment="center",
                horizontalalignment="center",
                fontsize=legend_fontsize,
                path_effects=legend_fontoutline,
            )
    elif legend_loc not in {None, "none"}:
        for label in cats:
            ax.scatter([], [], c=palette[label], label=label)
        if legend_loc == "right margin":
            ax.legend(
                frameon=False,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                ncol=(1 if len(cats) <= 14 else 2 if len(cats) <= 30 else 3),
                fontsize=legend_fontsize,
            )
        else:
            ax.legend(loc=legend_loc, fontsize=legend_fontsize)


def _get_basis(adata: AnnData, basis: str) -> np.ndarray:
    """Get array for basis from anndata. Just tries to add 'X_'."""
    if basis in adata.obsm:
        return adata.obsm[basis]
    elif f"X_{basis}" in adata.obsm:
        return adata.obsm[f"X_{basis}"]
    else:
        msg = f"Could not find {basis!r} or 'X_{basis}' in .obsm"
        raise KeyError(msg)


def _get_color_source_vector(
    adata: AnnData,
    value_to_plot: str,
    *,
    mask_obs: NDArray[np.bool_] | None = None,
    use_raw: bool = False,
    gene_symbols: str | None = None,
    layer: str | None = None,
    groups: Sequence[str] | None = None,
) -> np.ndarray | pd.api.extensions.ExtensionArray:
    """Get array from adata that colors will be based on."""
    if value_to_plot is None:
        # Points will be plotted with `na_color`. Ideally this would work
        # with the "bad color" in a color map but that throws a warning. Instead
        # _color_vector handles this.
        # https://github.com/matplotlib/matplotlib/issues/18294
        return np.broadcast_to(np.nan, adata.n_obs)
    if (
        gene_symbols is not None
        and value_to_plot not in adata.obs.columns
        and value_to_plot not in adata.var_names
    ):
        # We should probably just make an index for this, and share it over runs
        # TODO: Throw helpful error if this doesn't work
        value_to_plot = adata.var.index[adata.var[gene_symbols] == value_to_plot][0]
    if use_raw and value_to_plot not in adata.obs.columns:
        values = adata.raw.obs_vector(value_to_plot)
    else:
        values = adata.obs_vector(value_to_plot, layer=layer)
    if mask_obs is not None:
        values = values.copy()
        values[~mask_obs] = np.nan
    if groups and isinstance(values, pd.Categorical):
        values = values.remove_categories(values.categories.difference(groups))
    return values


def _get_palette(adata, values_key: str, palette=None):
    color_key = f"{values_key}_colors"
    if adata.obs[values_key].dtype == bool:
        values = pd.Categorical(adata.obs[values_key].astype(str))
    else:
        values = pd.Categorical(adata.obs[values_key])
    if palette:
        _utils._set_colors_for_categorical_obs(adata, values_key, palette)
    elif color_key not in adata.uns or len(adata.uns[color_key]) < len(
        values.categories
    ):
        #  set a default palette in case that no colors or too few colors are found
        _utils._set_default_colors_for_categorical_obs(adata, values_key)
    else:
        _utils._validate_palette(adata, values_key)
    return dict(
        zip(
            values.categories,
            adata.uns[color_key][: len(values.categories)],
            strict=True,
        )
    )


def _color_vector(
    adata: AnnData,
    values_key: str | None,
    *,
    values: np.ndarray | pd.api.extensions.ExtensionArray,
    palette: str | Sequence[str] | Cycler | None,
    na_color: ColorLike = "lightgray",
) -> tuple[np.ndarray | pd.api.extensions.ExtensionArray, Literal["cat", "na", "cont"]]:
    """Map array of values to array of hex (plus alpha) codes.

    For categorical data, the return value is list of colors taken
    from the category palette or from the given `palette` value.

    For continuous values, the input array is returned (may change in future).
    """
    ###
    # when plotting, the color of the dots is determined for each plot
    # the data is either categorical or continuous and the data could be in
    # 'obs' or in 'var'
    to_hex = partial(colors.to_hex, keep_alpha=True)
    if values_key is None:
        return np.broadcast_to(to_hex(na_color), adata.n_obs), "na"
    if values.dtype == bool:
        values = pd.Categorical(values.astype(str))
    elif not isinstance(values, pd.Categorical):
        return values, "cont"

    color_map = {
        k: to_hex(v)
        for k, v in _get_palette(adata, values_key, palette=palette).items()
    }
    # If color_map does not have unique values, this can be slow as the
    # result is not categorical
    if Version(pd.__version__) < Version("2.1.0"):
        color_vector = pd.Categorical(values.map(color_map))
    else:
        color_vector = pd.Categorical(values.map(color_map, na_action="ignore"))
    # Set color to 'missing color' for all missing values
    if color_vector.isna().any():
        color_vector = color_vector.add_categories([to_hex(na_color)])
        color_vector = color_vector.fillna(to_hex(na_color))
    return color_vector, "cat"


def _basis2name(basis):
    """Convert the 'basis' into the proper name."""
    component_name = (
        "DC"
        if basis == "diffmap"
        else "tSNE"
        if basis == "tsne"
        else "UMAP"
        if basis == "umap"
        else "PC"
        if basis == "pca"
        else basis.replace("draw_graph_", "").upper()
        if "draw_graph" in basis
        else basis
    )
    return component_name


def _check_spot_size(spatial_data: Mapping | None, spot_size: float | None) -> float:
    """Resolve spot_size value.

    This is a required argument for spatial plots.
    """
    if spatial_data is None and spot_size is None:
        msg = (
            "When .uns['spatial'][library_id] does not exist, spot_size must be "
            "provided directly."
        )
        raise ValueError(msg)
    elif spot_size is None:
        return spatial_data["scalefactors"]["spot_diameter_fullres"]
    else:
        return spot_size


def _check_scale_factor(
    spatial_data: Mapping | None,
    img_key: str | None,
    scale_factor: float | None,
) -> float:
    """Resolve scale_factor, defaults to 1."""
    if scale_factor is not None:
        return scale_factor
    elif spatial_data is not None and img_key is not None:
        return spatial_data["scalefactors"][f"tissue_{img_key}_scalef"]
    else:
        return 1.0


def _check_spatial_data(
    uns: Mapping, library_id: str | None | Empty
) -> tuple[str | None, Mapping | None]:
    """Given a mapping, try and extract a library id/ mapping with spatial data.

    Assumes this is `.uns` from how we parse visium data.
    """
    spatial_mapping = uns.get("spatial", {})
    if library_id is _empty:
        if len(spatial_mapping) > 1:
            msg = (
                "Found multiple possible libraries in `.uns['spatial']. Please specify."
                f" Options are:\n\t{list(spatial_mapping.keys())}"
            )
            raise ValueError(msg)
        elif len(spatial_mapping) == 1:
            library_id = next(iter(spatial_mapping.keys()))
        else:
            library_id = None
    spatial_data = spatial_mapping[library_id] if library_id is not None else None
    return library_id, spatial_data


def _check_img(
    spatial_data: Mapping | None,
    img: np.ndarray | None,
    img_key: None | str | Empty,
    *,
    bw: bool = False,
) -> tuple[np.ndarray | None, str | None]:
    """Resolve image for spatial plots."""
    if img is None and spatial_data is not None and img_key is _empty:
        img_key = next(
            (k for k in ["hires", "lowres"] if k in spatial_data["images"]),
        )  # Throws StopIteration Error if keys not present
    if img is None and spatial_data is not None and img_key is not None:
        img = spatial_data["images"][img_key]
    if bw:
        img = np.dot(img[..., :3], [0.2989, 0.5870, 0.1140])
    return img, img_key


def _check_crop_coord(
    crop_coord: tuple | None,
    scale_factor: float,
) -> tuple[float, float, float, float]:
    """Handle cropping with image or basis."""
    if crop_coord is None:
        return None
    if len(crop_coord) != 4:
        msg = "Invalid crop_coord of length {len(crop_coord)}(!=4)"
        raise ValueError(msg)
    crop_coord = tuple(c * scale_factor for c in crop_coord)
    return crop_coord


def _check_na_color(
    na_color: ColorLike | None, *, img: np.ndarray | None = None
) -> ColorLike:
    if na_color is None:
        na_color = (0.0, 0.0, 0.0, 0.0) if img is not None else "lightgray"
    return na_color


def _broadcast_args(*args):
    """Broadcasts arguments to a common length."""
    lens = [len(arg) for arg in args]
    longest = max(lens)
    if not (set(lens) == {1, longest} or set(lens) == {longest}):
        msg = f"Could not broadcast together arguments with shapes: {lens}."
        raise ValueError(msg)
    return [[arg[0] for _ in range(longest)] if len(arg) == 1 else arg for arg in args]
