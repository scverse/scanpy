from typing import Optional, Literal, Union, List  # Special
from typing import Collection, Iterable, Sequence  # ABCs
from typing import Tuple  # classes

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, ListedColormap, is_color_like
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

from ._tools.scatterplots import (
    _get_color_values,
    _get_data_points,
    gene_symbol_column,
    VMinMax,
    _get_vmin_vmax,
    _basis2name,
    _add_legend_or_colorbar,
)
from ._utils import (
    _FontSize,
    _FontWeight,
    ColorLike,
    savefig_or_show,
    make_grid_spec,
    _AxesSubplot,
    make_projection_available,
)
from .. import settings

_Aes = Literal['x', 'y', 'z', 'color', 'alpha', 'size', 'layer', 'basis']
_Basis = Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr']


def scatter(
    adata: AnnData,
    x: Union[str, Collection[str], np.ndarray, None] = None,
    y: Union[str, Collection[str], np.ndarray, None] = None,
    basis: Union[_Basis, Collection[_Basis], None] = None,
    *,
    # additional point aesthetics
    z: Union[str, Collection[str], np.ndarray, None] = None,
    color: Union[str, Collection[str], np.ndarray, None] = None,
    size: Union[str, Collection[str], np.ndarray, None] = None,
    alpha: Union[str, Collection[str], np.ndarray, None] = None,
    # TODO
    groups: Optional[str, Iterable[str]] = None,
    sort_order: bool = True,
    # mapping
    use_raw: Optional[bool] = None,
    layer: Union[str, Collection[str], None] = None,
    color_map: Union[str, Colormap] = None,
    palette: Union[Cycler, ListedColormap, ColorLike, Sequence[ColorLike]] = None,
    projection: Literal['2d', '3d'] = '2d',  # TODO: z
    gene_symbols: Optional[str],
    # vminmax
    vmax: Union[VMinMax, Sequence[VMinMax], None] = None,
    vmin: Union[VMinMax, Sequence[VMinMax], None] = None,
    # legend
    legend_loc: str = 'right margin',
    legend_fontsize: Union[int, float, _FontSize, None] = None,
    legend_fontweight: Union[int, _FontWeight, None] = None,
    legend_fontoutline: float = None,
    # global figure params
    title: Optional[str] = None,
    frameon: Optional[bool] = None,
    figsize: Optional[Tuple[int, int]] = None,
    # gridspec
    ncols: Optional[int] = None,
    nrows: Optional[int] = None,
    hspace: Optional[float] = None,
    wspace: Optional[float] = None,
    # show save ax
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    ax: Union[Axes, _AxesSubplot, None] = None,
    return_fig: bool = False,
) -> Union[Figure, Axes, List[Axes]]:
    params = locals()
    if use_raw is None and adata.raw is not None:
        use_raw = True
    make_projection_available(projection)
    ax_parent = ax
    del ax

    # extract param sets for plots
    multi = {
        k: v
        for k, v in dict(
            x=x, y=y, color=color, alpha=alpha, size=size, layer=layer, basis=basis
        ).items()
        if not (
            v is None
            or isinstance(v, str)
            or (not isinstance(v, np.ndarray) and len(v) > 1)
        )
    }
    for k in multi.keys() | {'ax'}:
        del params[k]
    n_plots = sum(map(len, multi.values())) or 1
    if n_plots > 1:
        if nrows is None:
            nrows = 1 if ncols is None else np.ceil(n_plots / ncols)
        if ncols is None:
            ncols = np.ceil(n_plots / nrows)
        fig, grid_spec = make_grid_spec(
            ax_parent or figsize, nrows, ncols, wspace, hspace, frameon=frameon
        )
        axes = []
        multi_params = [{k: v} for k, vs in multi.items() for v in vs]
        for i, (gs, ps) in enumerate(zip(grid_spec, multi_params)):
            a = fig.add_subplot(gs, projection=projection)
            do_scatter(a, i, multi=True, **params, **ps)
        ax = axes
    else:  # do a single plot
        if ax_parent is None:
            fig = plt.figure(figsize=figsize, frameon=frameon)
            ax = fig.add_subplot(projection=projection)
        else:
            fig = ax_parent.figure
            ax = ax_parent

        do_scatter(ax, 0, multi=False, **params)

    if return_fig:
        return fig
    savefig_or_show(basis or "scatter", show=show, save=save)
    if show is False:
        return ax


def do_scatter(
    ax: Union[Axes, Axes3D],
    index: int,
    *,
    multi: bool,
    adata: AnnData,
    x: Union[str, np.ndarray, None] = None,
    y: Union[str, np.ndarray, None] = None,
    z: Union[str, np.ndarray, None] = None,
    basis: Optional[_Basis] = None,
    color: Union[str, np.ndarray, None] = None,
    size: Union[str, np.ndarray, None] = None,
    alpha: Union[str, np.ndarray, None] = None,
    # TODO
    groups: Union[str, Iterable[str]] = None,
    sort_order: bool = True,
    # mapping
    use_raw: Optional[bool] = None,
    layer: Optional[str] = None,
    color_map: Union[str, Colormap] = None,
    palette: Union[Cycler, ListedColormap, ColorLike, Sequence[ColorLike]] = None,
    projection: Literal['2d', '3d'] = '2d',  # TODO: z
    gene_symbols: Optional[str],
    # vminmax
    vmax: Union[VMinMax, Sequence[VMinMax], None] = None,
    vmin: Union[VMinMax, Sequence[VMinMax], None] = None,
    # legend
    legend_loc: str = 'right margin',
    legend_fontsize: Union[int, float, _FontSize, None] = None,
    legend_fontweight: Union[int, _FontWeight, None] = None,
    legend_fontoutline: float = None,
) -> None:
    dims = [x, y, z] if projection == '3d' else [x, y]
    basis_key = f"X_{basis}"
    name = _basis2name(basis)
    data_points = [None] * len(dims)
    axis_labels = [None] * len(dims)
    # TODO: restructure _get_data_points?
    for i, d in enumerate(dims):
        if isinstance(d, np.ndarray):
            data_points[i] = d
        elif isinstance(d, (type(None), int)) or d in '123':
            if d in '123':
                d = int(d) - 1
            elif d is None:
                d = i
            data_points[i] = adata.obsm[basis_key][:, d]
            axis_labels[i] = f"{name}{i+1}"
        else:  # isinstance(d, str)
            d = gene_symbol_column(adata, gene_symbols, d)
            data_points[i] = adata.obs_vector(d, layer=layer)
            axis_labels[i] = d
    data_points = np.array(data_points)

    if not is_color_like(color):
        if sort_order:
            sort = np.argsort(color)
            color = color[sort]
            data_points = data_points[:, sort]
        color_vector = color
        categorical = pd.is_categorical(color)
    else:
        color_vector, categorical = _get_color_values(
            adata,
            color,
            layer=layer,
            groups=groups,
            palette=palette,
            use_raw=use_raw,
            gene_symbols=gene_symbols,
        )
    if categorical:
        vmin = vmax = None
    else:
        vmin, vmax = _get_vmin_vmax(vmin, vmax, index, color_vector)

    ax.scatter(
        *data_points,
        marker='.',
        c=color_vector,
        alpha=alpha,
        s=size,
        edgecolors='none',
        cmap=color_map,
        rasterized=settings._vector_friendly,
        vmin=vmin,
        vmax=vmax,
    )
    for n, l in zip("xyz", axis_labels):
        getattr(ax, f"{n}axis").set_label_lext(l)
    _add_legend_or_colorbar(
        adata,
        ax,
        cax,
        categorical,
        color,
        legend_loc,
        data_points,
        legend_fontweight,
        legend_fontsize,
        legend_fontoutline,
        groups,
        multi,
    )
