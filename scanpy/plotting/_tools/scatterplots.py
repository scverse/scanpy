import collections.abc as cabc
from copy import copy
from numbers import Integral
from itertools import combinations, product
from typing import (
    Collection,
    Union,
    Optional,
    Sequence,
    Any,
    Mapping,
    List,
    Tuple,
)
from warnings import warn

import numpy as np
import pandas as pd
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pandas.api.types import is_categorical_dtype
from matplotlib import pyplot as pl, colors
from matplotlib.cm import get_cmap
from matplotlib import rcParams
from matplotlib import patheffects
from matplotlib.colors import Colormap, Normalize
from functools import partial

from .. import _utils
from .._utils import (
    _IGraphLayout,
    _FontWeight,
    _FontSize,
    ColorLike,
    VBound,
    circles,
    check_projection,
    check_colornorm,
)
from .._docs import (
    doc_adata_color_etc,
    doc_edges_arrows,
    doc_scatter_embedding,
    doc_scatter_spatial,
    doc_show_save_ax,
)
from ... import logging as logg
from ..._settings import settings
from ..._utils import sanitize_anndata, _doc_params, Empty, _empty
from ..._compat import Literal


@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def embedding(
    adata: AnnData,
    basis: str,
    *,
    color: Union[str, Sequence[str], None] = None,
    gene_symbols: Optional[str] = None,
    use_raw: Optional[bool] = None,
    sort_order: bool = True,
    edges: bool = False,
    edges_width: float = 0.1,
    edges_color: Union[str, Sequence[float], Sequence[str]] = 'grey',
    neighbors_key: Optional[str] = None,
    arrows: bool = False,
    arrows_kwds: Optional[Mapping[str, Any]] = None,
    groups: Optional[str] = None,
    components: Union[str, Sequence[str]] = None,
    dimensions: Optional[Union[Tuple[int, int], Sequence[Tuple[int, int]]]] = None,
    layer: Optional[str] = None,
    projection: Literal['2d', '3d'] = '2d',
    scale_factor: Optional[float] = None,
    color_map: Union[Colormap, str, None] = None,
    cmap: Union[Colormap, str, None] = None,
    palette: Union[str, Sequence[str], Cycler, None] = None,
    na_color: ColorLike = "lightgray",
    na_in_legend: bool = True,
    size: Union[float, Sequence[float], None] = None,
    frameon: Optional[bool] = None,
    legend_fontsize: Union[int, float, _FontSize, None] = None,
    legend_fontweight: Union[int, _FontWeight] = 'bold',
    legend_loc: str = 'right margin',
    legend_fontoutline: Optional[int] = None,
    colorbar_loc: Optional[str] = "right",
    vmax: Union[VBound, Sequence[VBound], None] = None,
    vmin: Union[VBound, Sequence[VBound], None] = None,
    vcenter: Union[VBound, Sequence[VBound], None] = None,
    norm: Union[Normalize, Sequence[Normalize], None] = None,
    add_outline: Optional[bool] = False,
    outline_width: Tuple[float, float] = (0.3, 0.05),
    outline_color: Tuple[str, str] = ('black', 'white'),
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: Optional[float] = None,
    title: Union[str, Sequence[str], None] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Axes] = None,
    return_fig: Optional[bool] = None,
    **kwargs,
) -> Union[Figure, Axes, None]:
    """\
    Scatter plot for user specified embedding basis (e.g. umap, pca, etc)

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
    args_3d = dict(projection='3d') if projection == '3d' else {}

    # Figure out if we're using raw
    if use_raw is None:
        # check if adata.raw is set
        use_raw = layer is None and adata.raw is not None
    if use_raw and layer is not None:
        raise ValueError(
            "Cannot use both a layer and the raw representation. Was passed:"
            f"use_raw={use_raw}, layer={layer}."
        )
    if use_raw and adata.raw is None:
        raise ValueError(
            "`use_raw` is set to True but AnnData object does not have raw. "
            "Please check."
        )

    if isinstance(groups, str):
        groups = [groups]

    # Color map
    if color_map is not None:
        if cmap is not None:
            raise ValueError("Cannot specify both `color_map` and `cmap`.")
        else:
            cmap = color_map
    cmap = copy(get_cmap(cmap))
    cmap.set_bad(na_color)
    kwargs["cmap"] = cmap
    # Prevents warnings during legend creation
    na_color = colors.to_hex(na_color, keep_alpha=True)

    if 'edgecolor' not in kwargs:
        # by default turn off edge color. Otherwise, for
        # very small sizes the edge will not reduce its size
        # (https://github.com/scverse/scanpy/issues/293)
        kwargs['edgecolor'] = 'none'

    # Vectorized arguments

    # turn color into a python list
    color = [color] if isinstance(color, str) or color is None else list(color)
    if title is not None:
        # turn title into a python list if not None
        title = [title] if isinstance(title, str) else list(title)

    # turn vmax and vmin into a sequence
    if isinstance(vmax, str) or not isinstance(vmax, cabc.Sequence):
        vmax = [vmax]
    if isinstance(vmin, str) or not isinstance(vmin, cabc.Sequence):
        vmin = [vmin]
    if isinstance(vcenter, str) or not isinstance(vcenter, cabc.Sequence):
        vcenter = [vcenter]
    if isinstance(norm, Normalize) or not isinstance(norm, cabc.Sequence):
        norm = [norm]

    # Size
    if 's' in kwargs and size is None:
        size = kwargs.pop('s')
    if size is not None:
        # check if size is any type of sequence, and if so
        # set as ndarray
        if (
            size is not None
            and isinstance(size, (cabc.Sequence, pd.Series, np.ndarray))
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
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02

    if components is not None:
        color, dimensions = list(zip(*product(color, dimensions)))

    color, dimensions = _broadcast_args(color, dimensions)

    # 'color' is a list of names that want to be plotted.
    # Eg. ['Gene1', 'louvain', 'Gene2'].
    # component_list is a list of components [[0,1], [1,2]]
    if (
        not isinstance(color, str)
        and isinstance(color, cabc.Sequence)
        and len(color) > 1
    ) or len(dimensions) > 1:
        if ax is not None:
            raise ValueError(
                "Cannot specify `ax` when plotting multiple panels "
                "(each for a given value of 'color')."
            )

        # each plot needs to be its own panel
        fig, grid = _panel_grid(hspace, wspace, ncols, len(color))
    else:
        grid = None
        if ax is None:
            fig = pl.figure()
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
    for count, (value_to_plot, dims) in enumerate(zip(color, dimensions)):
        color_source_vector = _get_color_source_vector(
            adata,
            value_to_plot,
            layer=layer,
            use_raw=use_raw,
            gene_symbols=gene_symbols,
            groups=groups,
        )
        color_vector, categorical = _color_vector(
            adata,
            value_to_plot,
            color_source_vector,
            palette=palette,
            na_color=na_color,
        )

        # Order points
        order = slice(None)
        if sort_order is True and value_to_plot is not None and categorical is False:
            # Higher values plotted on top, null values on bottom
            order = np.argsort(-color_vector, kind="stable")[::-1]
        elif sort_order and categorical:
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
            ax = pl.subplot(grid[count], **args_3d)
            axs.append(ax)
        if not (settings._frameon if frameon is None else frameon):
            ax.axis('off')
        if title is None:
            if value_to_plot is not None:
                ax.set_title(value_to_plot)
            else:
                ax.set_title('')
        else:
            try:
                ax.set_title(title[count])
            except IndexError:
                logg.warning(
                    "The title list is shorter than the number of panels. "
                    "Using 'color' value instead for some plots."
                )
                ax.set_title(value_to_plot)

        if not categorical:
            vmin_float, vmax_float, vcenter_float, norm_obj = _get_vboundnorm(
                vmin, vmax, vcenter, norm, count, color_vector
            )
            normalize = check_colornorm(
                vmin_float,
                vmax_float,
                vcenter_float,
                norm_obj,
            )
        else:
            normalize = None

        # make the scatter plot
        if projection == '3d':
            cax = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                coords[:, 2],
                marker=".",
                c=color_vector,
                rasterized=settings._vector_friendly,
                norm=normalize,
                **kwargs,
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
                kwargs['edgecolor'] = 'none'

                # remove alpha for outline
                alpha = kwargs.pop('alpha') if 'alpha' in kwargs else None

                ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    s=bg_size,
                    marker=".",
                    c=bg_color,
                    rasterized=settings._vector_friendly,
                    norm=normalize,
                    **kwargs,
                )
                ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    s=gap_size,
                    marker=".",
                    c=gap_color,
                    rasterized=settings._vector_friendly,
                    norm=normalize,
                    **kwargs,
                )
                # if user did not set alpha, set alpha to 0.7
                kwargs['alpha'] = 0.7 if alpha is None else alpha

            cax = scatter(
                coords[:, 0],
                coords[:, 1],
                marker=".",
                c=color_vector,
                rasterized=settings._vector_friendly,
                norm=normalize,
                **kwargs,
            )

        # remove y and x ticks
        ax.set_yticks([])
        ax.set_xticks([])
        if projection == '3d':
            ax.set_zticks([])

        # set default axis_labels
        name = _basis2name(basis)
        axis_labels = [name + str(d + 1) for d in dims]

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        if projection == '3d':
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[2], labelpad=-7)
        ax.autoscale_view()

        if edges:
            _utils.plot_edges(ax, adata, basis, edges_width, edges_color, neighbors_key)
        if arrows:
            _utils.plot_arrows(ax, adata, basis, arrows_kwds)

        if value_to_plot is None:
            # if only dots were plotted without an associated value
            # there is not need to plot a legend or a colorbar
            continue

        if legend_fontoutline is not None:
            path_effect = [
                patheffects.withStroke(linewidth=legend_fontoutline, foreground='w')
            ]
        else:
            path_effect = None

        # Adding legends
        if categorical:
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
            pl.colorbar(
                cax, ax=ax, pad=0.01, fraction=0.08, aspect=30, location=colorbar_loc
            )

    if return_fig is True:
        return fig
    axs = axs if grid else ax
    _utils.savefig_or_show(basis, show=show, save=save)
    if show is False:
        return axs


def _panel_grid(hspace, wspace, ncols, num_panels):
    from matplotlib import gridspec

    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    # each panel will have the size of rcParams['figure.figsize']
    fig = pl.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
            n_panels_y * rcParams['figure.figsize'][1],
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
    norm: Sequence[Normalize],
    index: int,
    color_vector: Sequence[float],
) -> Tuple[Union[float, None], Union[float, None]]:
    """
    Evaluates the value of vmin, vmax and vcenter, which could be a
    str in which case is interpreted as a percentile and should
    be specified in the form 'pN' where N is the percentile.
    Eg. for a percentile of 85 the format would be 'p85'.
    Floats are accepted as p99.9

    Alternatively, vmin/vmax could be a function that is applied to
    the list of color values (`color_vector`).  E.g.

    def my_vmax(color_vector): np.percentile(color_vector, p=80)


    Parameters
    ----------
    index
        This index of the plot
    color_vector
        List or values for the plot

    Returns
    -------

    (vmin, vmax, vcenter, norm) containing None or float values for
    vmin, vmax, vcenter and matplotlib.colors.Normalize  or None for norm.

    """
    out = []
    for v_name, v in [('vmin', vmin), ('vmax', vmax), ('vcenter', vcenter)]:
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
            if isinstance(v_value, str) and v_value.startswith('p'):
                try:
                    float(v_value[1:])
                except ValueError:
                    logg.error(
                        f"The parameter {v_name}={v_value} for plot number {index + 1} is not valid. "
                        f"Please check the correct format for percentiles."
                    )
                # interpret value of vmin/vmax as quantile with the following syntax 'p99.9'
                v_value = np.nanpercentile(color_vector, q=float(v_value[1:]))
            elif callable(v_value):
                # interpret vmin/vmax as function
                v_value = v_value(color_vector)
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
    import inspect

    params = inspect.signature(embedding).parameters.copy()
    wrapper_sig = inspect.signature(wrapper)
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
def umap(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in UMAP basis.

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

    See also
    --------
    tl.umap
    """
    return embedding(adata, 'umap', **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def tsne(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in tSNE basis.

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

    See also
    --------
    tl.tsne
    """
    return embedding(adata, 'tsne', **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def diffmap(adata, **kwargs) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in Diffusion Map basis.

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

    See also
    --------
    tl.diffmap
    """
    return embedding(adata, 'diffmap', **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    edges_arrows=doc_edges_arrows,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def draw_graph(
    adata: AnnData, *, layout: Optional[_IGraphLayout] = None, **kwargs
) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in graph-drawing basis.

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

    See also
    --------
    tl.draw_graph
    """
    if layout is None:
        layout = str(adata.uns['draw_graph']['params']['layout'])
    basis = 'draw_graph_' + layout
    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError(
            'Did not find {} in adata.obs. Did you compute layout {}?'.format(
                'draw_graph_' + layout, layout
            )
        )

    return embedding(adata, basis, **kwargs)


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def pca(
    adata,
    *,
    annotate_var_explained: bool = False,
    show: Optional[bool] = None,
    return_fig: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    **kwargs,
) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in PCA coordinates.

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

    See also
    --------
    tl.pca
    pp.pca
    """
    if not annotate_var_explained:
        return embedding(
            adata, 'pca', show=show, return_fig=return_fig, save=save, **kwargs
        )
    else:

        if 'pca' not in adata.obsm.keys() and 'X_pca' not in adata.obsm.keys():
            raise KeyError(
                f"Could not find entry in `obsm` for 'pca'.\n"
                f"Available keys are: {list(adata.obsm.keys())}."
            )

        label_dict = {
            'PC{}'.format(i + 1): 'PC{} ({}%)'.format(i + 1, round(v * 100, 2))
            for i, v in enumerate(adata.uns['pca']['variance_ratio'])
        }

        if return_fig is True:
            # edit axis labels in returned figure
            fig = embedding(adata, 'pca', return_fig=return_fig, **kwargs)
            for ax in fig.axes:
                ax.set_xlabel(label_dict[ax.xaxis.get_label().get_text()])
                ax.set_ylabel(label_dict[ax.yaxis.get_label().get_text()])
            return fig

        else:
            # get the axs, edit the labels and apply show and save from user
            axs = embedding(adata, 'pca', show=False, save=False, **kwargs)
            if isinstance(axs, list):
                for ax in axs:
                    ax.set_xlabel(label_dict[ax.xaxis.get_label().get_text()])
                    ax.set_ylabel(label_dict[ax.yaxis.get_label().get_text()])
            else:
                axs.set_xlabel(label_dict[axs.xaxis.get_label().get_text()])
                axs.set_ylabel(label_dict[axs.yaxis.get_label().get_text()])
            _utils.savefig_or_show('pca', show=show, save=save)
            if show is False:
                return axs


@_wraps_plot_scatter
@_doc_params(
    adata_color_etc=doc_adata_color_etc,
    scatter_spatial=doc_scatter_spatial,
    scatter_bulk=doc_scatter_embedding,
    show_save_ax=doc_show_save_ax,
)
def spatial(
    adata,
    *,
    basis: str = "spatial",
    img: Union[np.ndarray, None] = None,
    img_key: Union[str, None, Empty] = _empty,
    library_id: Union[str, Empty] = _empty,
    crop_coord: Tuple[int, int, int, int] = None,
    alpha_img: float = 1.0,
    bw: Optional[bool] = False,
    size: float = 1.0,
    scale_factor: Optional[float] = None,
    spot_size: Optional[float] = None,
    na_color: Optional[ColorLike] = None,
    show: Optional[bool] = None,
    return_fig: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    **kwargs,
) -> Union[Axes, List[Axes], None]:
    """\
    Scatter plot in spatial coordinates.

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

    >>> adata = sc.datasets.visium_sge("Targeted_Visium_Human_Glioblastoma_Pan_Cancer")
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.pl.spatial(adata, color="log1p_n_genes_by_counts")

    See Also
    --------
    :func:`scanpy.datasets.visium_sge`
        Example visium data.
    :tutorial:`spatial/basic-analysis`
        Tutorial on spatial analysis.
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

    if bw:
        cmap_img = "gray"
    else:
        cmap_img = None
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
    _utils.savefig_or_show('show', show=show, save=save)
    if show is False or return_fig is True:
        return axs


# Helpers
def _components_to_dimensions(
    components: Optional[Union[str, Collection[str]]],
    dimensions: Optional[Union[Collection[int], Collection[Collection[int]]]],
    *,
    projection: Literal["2d", "3d"] = "2d",
    total_dims: int,
) -> List[Collection[int]]:
    """Normalize components/ dimensions args for embedding plots."""
    # TODO: Deprecate components kwarg
    ndims = {"2d": 2, "3d": 3}[projection]
    if components is None and dimensions is None:
        dimensions = [tuple(i for i in range(ndims))]
    elif components is not None and dimensions is not None:
        raise ValueError("Cannot provide both dimensions and components")

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


def _add_categorical_legend(
    ax,
    color_source_vector,
    palette: dict,
    legend_loc: str,
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
            raise NotImplementedError(
                "No fallback for null labels has been defined if NA already in categories."
            )
        color_source_vector = color_source_vector.add_categories("NA").fillna("NA")
        palette = palette.copy()
        palette["NA"] = na_color
    cats = color_source_vector.categories

    if multi_panel is True:
        # Shrink current axis by 10% to fit legend and match
        # size of plots that are not categorical
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.91, box.height])

    if legend_loc == 'right margin':
        for label in cats:
            ax.scatter([], [], c=palette[label], label=label)
        ax.legend(
            frameon=False,
            loc='center left',
            bbox_to_anchor=(1, 0.5),
            ncol=(1 if len(cats) <= 14 else 2 if len(cats) <= 30 else 3),
            fontsize=legend_fontsize,
        )
    elif legend_loc == 'on data':
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
                verticalalignment='center',
                horizontalalignment='center',
                fontsize=legend_fontsize,
                path_effects=legend_fontoutline,
            )


def _get_basis(adata: AnnData, basis: str) -> np.ndarray:
    """Get array for basis from anndata. Just tries to add 'X_'."""
    if basis in adata.obsm:
        return adata.obsm[basis]
    elif f"X_{basis}" in adata.obsm:
        return adata.obsm[f"X_{basis}"]
    else:
        raise KeyError(f"Could not find '{basis}' or 'X_{basis}' in .obsm")


def _get_color_source_vector(
    adata, value_to_plot, use_raw=False, gene_symbols=None, layer=None, groups=None
):
    """
    Get array from adata that colors will be based on.
    """
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
        value_to_plot = adata.var.index[adata.var[gene_symbols] == value_to_plot][
            0
        ]  # TODO: Throw helpful error if this doesn't work
    if use_raw and value_to_plot not in adata.obs.columns:
        values = adata.raw.obs_vector(value_to_plot)
    else:
        values = adata.obs_vector(value_to_plot, layer=layer)
    if groups and is_categorical_dtype(values):
        values = values.replace(values.categories.difference(groups), np.nan)
    return values


def _get_palette(adata, values_key: str, palette=None):
    color_key = f"{values_key}_colors"
    values = pd.Categorical(adata.obs[values_key])
    if palette:
        _utils._set_colors_for_categorical_obs(adata, values_key, palette)
    elif color_key not in adata.uns or len(adata.uns[color_key]) < len(
        values.categories
    ):
        #  set a default palette in case that no colors or few colors are found
        _utils._set_default_colors_for_categorical_obs(adata, values_key)
    else:
        _utils._validate_palette(adata, values_key)
    return dict(zip(values.categories, adata.uns[color_key]))


def _color_vector(
    adata, values_key: str, values, palette, na_color="lightgray"
) -> Tuple[np.ndarray, bool]:
    """
    Map array of values to array of hex (plus alpha) codes.

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
        return np.broadcast_to(to_hex(na_color), adata.n_obs), False
    if not is_categorical_dtype(values):
        return values, False
    else:  # is_categorical_dtype(values)
        color_map = {
            k: to_hex(v)
            for k, v in _get_palette(adata, values_key, palette=palette).items()
        }
        # If color_map does not have unique values, this can be slow as the
        # result is not categorical
        color_vector = pd.Categorical(values.map(color_map))

        # Set color to 'missing color' for all missing values
        if color_vector.isna().any():
            color_vector = color_vector.add_categories([to_hex(na_color)])
            color_vector = color_vector.fillna(to_hex(na_color))
        return color_vector, True


def _basis2name(basis):
    """
    converts the 'basis' into the proper name.
    """

    component_name = (
        'DC'
        if basis == 'diffmap'
        else 'tSNE'
        if basis == 'tsne'
        else 'UMAP'
        if basis == 'umap'
        else 'PC'
        if basis == 'pca'
        else basis.replace('draw_graph_', '').upper()
        if 'draw_graph' in basis
        else basis
    )
    return component_name


def _check_spot_size(
    spatial_data: Optional[Mapping], spot_size: Optional[float]
) -> float:
    """
    Resolve spot_size value.

    This is a required argument for spatial plots.
    """
    if spatial_data is None and spot_size is None:
        raise ValueError(
            "When .uns['spatial'][library_id] does not exist, spot_size must be "
            "provided directly."
        )
    elif spot_size is None:
        return spatial_data['scalefactors']['spot_diameter_fullres']
    else:
        return spot_size


def _check_scale_factor(
    spatial_data: Optional[Mapping],
    img_key: Optional[str],
    scale_factor: Optional[float],
) -> float:
    """Resolve scale_factor, defaults to 1."""
    if scale_factor is not None:
        return scale_factor
    elif spatial_data is not None and img_key is not None:
        return spatial_data['scalefactors'][f"tissue_{img_key}_scalef"]
    else:
        return 1.0


def _check_spatial_data(
    uns: Mapping, library_id: Union[Empty, None, str]
) -> Tuple[Optional[str], Optional[Mapping]]:
    """
    Given a mapping, try and extract a library id/ mapping with spatial data.

    Assumes this is `.uns` from how we parse visium data.
    """
    spatial_mapping = uns.get("spatial", {})
    if library_id is _empty:
        if len(spatial_mapping) > 1:
            raise ValueError(
                "Found multiple possible libraries in `.uns['spatial']. Please specify."
                f" Options are:\n\t{list(spatial_mapping.keys())}"
            )
        elif len(spatial_mapping) == 1:
            library_id = list(spatial_mapping.keys())[0]
        else:
            library_id = None
    if library_id is not None:
        spatial_data = spatial_mapping[library_id]
    else:
        spatial_data = None
    return library_id, spatial_data


def _check_img(
    spatial_data: Optional[Mapping],
    img: Optional[np.ndarray],
    img_key: Union[None, str, Empty],
    bw: bool = False,
) -> Tuple[Optional[np.ndarray], Optional[str]]:
    """
    Resolve image for spatial plots.
    """
    if img is None and spatial_data is not None and img_key is _empty:
        img_key = next(
            (k for k in ['hires', 'lowres'] if k in spatial_data['images']),
        )  # Throws StopIteration Error if keys not present
    if img is None and spatial_data is not None and img_key is not None:
        img = spatial_data["images"][img_key]
    if bw:
        img = np.dot(img[..., :3], [0.2989, 0.5870, 0.1140])
    return img, img_key


def _check_crop_coord(
    crop_coord: Optional[tuple],
    scale_factor: float,
) -> Tuple[float, float, float, float]:
    """Handle cropping with image or basis."""
    if crop_coord is None:
        return None
    if len(crop_coord) != 4:
        raise ValueError("Invalid crop_coord of length {len(crop_coord)}(!=4)")
    crop_coord = tuple(c * scale_factor for c in crop_coord)
    return crop_coord


def _check_na_color(
    na_color: Optional[ColorLike], *, img: Optional[np.ndarray] = None
) -> ColorLike:
    if na_color is None:
        if img is not None:
            na_color = (0.0, 0.0, 0.0, 0.0)
        else:
            na_color = "lightgray"
    return na_color


def _broadcast_args(*args):
    """Broadcasts arguments to a common length."""
    from itertools import repeat

    lens = [len(arg) for arg in args]
    longest = max(lens)
    if not (set(lens) == {1, longest} or set(lens) == {longest}):
        raise ValueError(f"Could not broadast together arguments with shapes: {lens}.")
    return list(
        [[arg[0] for _ in range(longest)] if len(arg) == 1 else arg for arg in args]
    )
