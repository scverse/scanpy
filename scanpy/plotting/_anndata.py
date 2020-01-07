"""Plotting functions for AnnData.
"""
import collections.abc as cabc
from itertools import product
from typing import Optional, Union, Mapping  # Special
from typing import Sequence, Collection, Iterable  # ABCs
from typing import Tuple, List  # Classes

import numpy as np
import pandas as pd
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib import patheffects
from matplotlib.colors import is_color_like, Colormap, ListedColormap

from .. import get
from .. import logging as logg
from .._settings import settings
from .._utils import sanitize_anndata, _doc_params
from .._compat import Literal
from . import _utils
from ._utils import scatter_base, scatter_group, setup_axes
from ._utils import ColorLike, _FontWeight, _FontSize
from ._docs import doc_scatter_basic, doc_show_save_ax, doc_common_plot_args


VALID_LEGENDLOCS = {
    'none',
    'right margin',
    'on data',
    'on data export',
    'best',
    'upper right',
    'upper left',
    'lower left',
    'lower right',
    'right',
    'center left',
    'center right',
    'lower center',
    'upper center',
    'center',
}

# TODO: is that all?
_Basis = Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr']
_VarNames = Union[str, Sequence[str]]


@_doc_params(scatter_temp=doc_scatter_basic, show_save_ax=doc_show_save_ax)
def scatter(
    adata: AnnData,
    x: Optional[str] = None,
    y: Optional[str] = None,
    color: Union[str, Collection[str]] = None,
    use_raw: Optional[bool] = None,
    layers: Union[str, Collection[str]] = None,
    sort_order: bool = True,
    alpha: Optional[float] = None,
    basis: Optional[_Basis] = None,
    groups: Union[str, Iterable[str]] = None,
    components: Union[str, Collection[str]] = None,
    projection: Literal['2d', '3d'] = '2d',
    legend_loc: str = 'right margin',
    legend_fontsize: Union[int, float, _FontSize, None] = None,
    legend_fontweight: Union[int, _FontWeight, None] = None,
    legend_fontoutline: float = None,
    color_map: Union[str, Colormap] = None,
    palette: Union[
        Cycler, ListedColormap, ColorLike, Sequence[ColorLike]
    ] = None,
    frameon: Optional[bool] = None,
    right_margin: Optional[float] = None,
    left_margin: Optional[float] = None,
    size: Union[int, float, None] = None,
    title: Optional[str] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    ax: Optional[Axes] = None,
):
    """\
    Scatter plot along observations or variables axes.

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
        Use `raw` attribute of `adata` if present.
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
    args = locals()
    if basis is not None:
        return _scatter_obs(**args)
    if x is None or y is None:
        raise ValueError('Either provide a `basis` or `x` and `y`.')
    if (
        (x in adata.obs.keys() or x in adata.var.index)
        and (y in adata.obs.keys() or y in adata.var.index)
        and (
            color is None
            or color in adata.obs.keys()
            or color in adata.var.index
        )
    ):
        return _scatter_obs(**args)
    if (
        (x in adata.var.keys() or x in adata.obs.index)
        and (y in adata.var.keys() or y in adata.obs.index)
        and (
            color is None
            or color in adata.var.keys()
            or color in adata.obs.index
        )
    ):
        adata_T = adata.T
        axs = _scatter_obs(
            adata=adata_T,
            **{name: val for name, val in args.items() if name != 'adata'},
        )
        # store .uns annotations that were added to the new adata object
        adata.uns = adata_T.uns
        return axs
    raise ValueError(
        '`x`, `y`, and potential `color` inputs must all '
        'come from either `.obs` or `.var`'
    )


def _scatter_obs(
    adata: AnnData,
    x=None,
    y=None,
    color=None,
    use_raw=None,
    layers=None,
    sort_order=True,
    alpha=None,
    basis=None,
    groups=None,
    components=None,
    projection: Literal['2d', '3d'] = '2d',
    legend_loc='right margin',
    legend_fontsize=None,
    legend_fontweight=None,
    legend_fontoutline=None,
    color_map=None,
    palette=None,
    frameon=None,
    right_margin=None,
    left_margin=None,
    size=None,
    title=None,
    show=None,
    save=None,
    ax=None,
):
    """See docstring of scatter."""
    sanitize_anndata(adata)
    from scipy.sparse import issparse

    if use_raw is None and adata.raw is not None:
        use_raw = True

    # Process layers
    if layers in ['X', None] or (
        isinstance(layers, str) and layers in adata.layers.keys()
    ):
        layers = (layers, layers, layers)
    elif isinstance(layers, cabc.Collection) and len(layers) == 3:
        layers = tuple(layers)
        for layer in layers:
            if layer not in adata.layers.keys() and layer not in ['X', None]:
                raise ValueError(
                    '`layers` should have elements that are '
                    'either None or in adata.layers.keys().'
                )
    else:
        raise ValueError(
            "`layers` should be a string or a collection of strings "
            f"with length 3, had value '{layers}'"
        )
    if use_raw and layers not in [('X', 'X', 'X'), (None, None, None)]:
        ValueError('`use_raw` must be `False` if layers are used.')

    if legend_loc not in VALID_LEGENDLOCS:
        raise ValueError(
            f'Invalid `legend_loc`, need to be one of: {VALID_LEGENDLOCS}.'
        )
    if components is None:
        components = '1,2' if '2d' in projection else '1,2,3'
    if isinstance(components, str):
        components = components.split(',')
    components = np.array(components).astype(int) - 1
    keys = (
        ['grey']
        if color is None
        else [color]
        if isinstance(color, str)
        else color
    )
    if title is not None and isinstance(title, str):
        title = [title]
    highlights = adata.uns['highlights'] if 'highlights' in adata.uns else []
    if basis is not None:
        try:
            # ignore the '0th' diffusion component
            if basis == 'diffmap':
                components += 1
            Y = adata.obsm['X_' + basis][:, components]
            # correct the component vector for use in labeling etc.
            if basis == 'diffmap':
                components -= 1
        except KeyError:
            raise KeyError(
                f'compute coordinates using visualization tool {basis} first'
            )
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
        raise ValueError('Either provide a `basis` or `x` and `y`.')

    if size is None:
        n = Y.shape[0]
        size = 120000 / n

    if legend_loc.startswith('on data') and legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']
    elif legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']

    palette_was_none = False
    if palette is None:
        palette_was_none = True
    if isinstance(palette, cabc.Sequence):
        if not is_color_like(palette[0]):
            palettes = palette
        else:
            palettes = [palette]
    else:
        palettes = [palette for _ in range(len(keys))]
    for i, palette in enumerate(palettes):
        palettes[i] = _utils.default_palette(palette)

    if basis is not None:
        component_name = (
            'DC'
            if basis == 'diffmap'
            else 'tSNE'
            if basis == 'tsne'
            else 'UMAP'
            if basis == 'umap'
            else 'PC'
            if basis == 'pca'
            else 'TriMap'
            if basis == 'trimap'
            else basis.replace('draw_graph_', '').upper()
            if 'draw_graph' in basis
            else basis
        )
    else:
        component_name = None
    axis_labels = (x, y) if component_name is None else None
    show_ticks = True if component_name is None else False

    # generate the colors
    color_ids = []
    categoricals = []
    colorbars = []
    for ikey, key in enumerate(keys):
        c = 'white'
        categorical = False  # by default, assume continuous or flat color
        colorbar = None
        # test whether we have categorial or continuous annotation
        if key in adata.obs_keys():
            if is_categorical_dtype(adata.obs[key]):
                categorical = True
            else:
                c = adata.obs[key]
        # coloring according to gene expression
        elif use_raw and adata.raw is not None and key in adata.raw.var_names:
            c = adata.raw.obs_vector(key)
        elif key in adata.var_names:
            c = adata.obs_vector(key, layer=layers[2])
        elif is_color_like(key):  # a flat color
            c = key
            colorbar = False
        else:
            raise ValueError(
                f'key {key!r} is invalid! pass valid observation annotation, '
                f'one of {adata.obs_keys()} or a gene name {adata.var_names}'
            )
        if colorbar is None:
            colorbar = not categorical
        colorbars.append(colorbar)
        if categorical:
            categoricals.append(ikey)
        color_ids.append(c)

    if right_margin is None and len(categoricals) > 0:
        if legend_loc == 'right margin':
            right_margin = 0.5
    if title is None and keys[0] is not None:
        title = [
            key.replace('_', ' ') if not is_color_like(key) else ''
            for key in keys
        ]

    axs = scatter_base(
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
    for i, ikey in enumerate(categoricals):
        palette = palettes[i]
        key = keys[ikey]
        _utils.add_colors_for_categorical_sample_annotation(
            adata, key, palette, force_update_colors=not palette_was_none
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
                        projection,
                        size=size,
                        alpha=alpha,
                    )
                    mask_remaining[mask] = False
                    if legend_loc.startswith('on data'):
                        add_centroid(centroids, name, Y, mask)
        else:
            groups = [groups] if isinstance(groups, str) else groups
            for name in groups:
                if name not in set(adata.obs[key].cat.categories):
                    raise ValueError(
                        f'{name!r} is invalid! specify valid name, '
                        f'one of {adata.obs[key].cat.categories}'
                    )
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
                        projection,
                        size=size,
                        alpha=alpha,
                    )
                    if legend_loc.startswith('on data'):
                        add_centroid(centroids, name, Y, mask)
                    mask_remaining[mask] = False
        if mask_remaining.sum() > 0:
            data = [Y[mask_remaining, 0], Y[mask_remaining, 1]]
            if projection == '3d':
                data.append(Y[mask_remaining, 2])
            axs[ikey].scatter(
                *data,
                marker='.',
                c='lightgrey',
                s=size,
                edgecolors='none',
                zorder=-1,
            )
        legend = None
        if legend_loc.startswith('on data'):
            if legend_fontweight is None:
                legend_fontweight = 'bold'
            if legend_fontoutline is not None:
                path_effect = [
                    patheffects.withStroke(
                        linewidth=legend_fontoutline, foreground='w'
                    )
                ]
            else:
                path_effect = None
            for name, pos in centroids.items():
                axs[ikey].text(
                    pos[0],
                    pos[1],
                    name,
                    weight=legend_fontweight,
                    verticalalignment='center',
                    horizontalalignment='center',
                    fontsize=legend_fontsize,
                    path_effects=path_effect,
                )

            all_pos = np.zeros((len(adata.obs[key].cat.categories), 2))
            for iname, name in enumerate(adata.obs[key].cat.categories):
                if name in centroids:
                    all_pos[iname] = centroids[name]
                else:
                    all_pos[iname] = [np.nan, np.nan]
            _utils._tmp_cluster_pos = all_pos
            if legend_loc == 'on data export':
                filename = settings.writedir / 'pos.csv'
                logg.warning(f'exporting label positions to {filename}')
                settings.writedir.mkdir(parents=True, exist_ok=True)
                np.savetxt(filename, all_pos, delimiter=',')
        elif legend_loc == 'right margin':
            legend = axs[ikey].legend(
                frameon=False,
                loc='center left',
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
        elif legend_loc != 'none':
            legend = axs[ikey].legend(
                frameon=False, loc=legend_loc, fontsize=legend_fontsize
            )
        if legend is not None:
            for handle in legend.legendHandles:
                handle.set_sizes([300.0])

    # draw a frame around the scatter
    frameon = settings._frameon if frameon is None else frameon
    if not frameon and x is None and y is None:
        for ax in axs:
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_frame_on(False)

    show = settings.autoshow if show is None else show
    _utils.savefig_or_show(
        'scatter' if basis is None else basis, show=show, save=save
    )
    if not show:
        return axs if len(keys) > 1 else axs[0]


def ranking(
    adata: AnnData,
    attr: Literal['var', 'obs', 'uns', 'varm', 'obsm'],
    keys: Union[str, Sequence[str]],
    dictionary=None,
    indices=None,
    labels=None,
    color='black',
    n_points=30,
    log=False,
    include_lowest=False,
    show=None,
):
    """\
    Plot rankings.

    See, for example, how this is used in pl.pca_ranking.

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
        keys = [f'{keys[:-1]}{i + 1}' for i in indices]
    else:
        if dictionary is None:
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
            if attr in {'var', 'varm'}
            else np.arange(scores.shape[0]).astype(str)
        )
    if isinstance(labels, str):
        labels = [labels + str(i + 1) for i in range(scores.shape[0])]
    if n_panels <= 5:
        n_rows, n_cols = 1, n_panels
    else:
        n_rows, n_cols = 2, int(n_panels / 2 + 0.5)
    fig = pl.figure(
        figsize=(
            n_cols * rcParams['figure.figsize'][0],
            n_rows * rcParams['figure.figsize'][1],
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
        pl.subplot(gs[iscore])
        order_scores = np.argsort(score)[::-1]
        if not include_lowest:
            indices = order_scores[: n_points + 1]
        else:
            indices = order_scores[: n_points // 2]
            neg_indices = order_scores[-(n_points - (n_points // 2)) :]
        txt_args = dict(
            color=color,
            rotation='vertical',
            verticalalignment='bottom',
            horizontalalignment='center',
            fontsize=8,
        )
        for ig, g in enumerate(indices):
            pl.text(ig, score[g], labels[g], **txt_args)
        if include_lowest:
            score_mid = (score[g] + score[neg_indices[0]]) / 2
            pl.text(len(indices), score_mid, '⋮', **txt_args)
            for ig, g in enumerate(neg_indices):
                pl.text(ig + len(indices) + 2, score[g], labels[g], **txt_args)
            pl.xticks([])
        pl.title(keys[iscore].replace('_', ' '))
        if n_panels <= 5 or iscore > n_cols:
            pl.xlabel('ranking')
        pl.xlim(-0.9, n_points + 0.9 + (1 if include_lowest else 0))
        score_min, score_max = (
            np.min(score[neg_indices if include_lowest else indices]),
            np.max(score[indices]),
        )
        pl.ylim(
            (0.95 if score_min > 0 else 1.05) * score_min,
            (1.05 if score_max > 0 else 0.95) * score_max,
        )
    show = settings.autoshow if show is None else show
    if not show:
        return gs


@_doc_params(show_save_ax=doc_show_save_ax)
def violin(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    groupby: Optional[str] = None,
    log: bool = False,
    use_raw: Optional[bool] = None,
    stripplot: bool = True,
    jitter: Union[float, bool] = True,
    size: int = 1,
    layer: Optional[str] = None,
    scale: Literal['area', 'count', 'width'] = 'width',
    order: Optional[Sequence[str]] = None,
    multi_panel: Optional[bool] = None,
    xlabel: str = '',
    rotation: Optional[float] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    ax: Optional[Axes] = None,
    **kwds,
):
    """\
    Violin plot.

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
        Use `raw` attribute of `adata` if present.
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
    scale
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
    rotation
        Rotation of xtick labels.
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.violinplot`.

    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object if `ax` is `None` else `None`.
    """
    import seaborn as sns  # Slow import, only import if called

    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None:
        use_raw = True
    if isinstance(keys, str):
        keys = [keys]
    if groupby is not None:
        obs_df = get.obs_df(
            adata, keys=[groupby] + keys, layer=layer, use_raw=use_raw
        )
    else:
        obs_df = get.obs_df(adata, keys=keys, layer=layer, use_raw=use_raw)
    if groupby is None:
        obs_tidy = pd.melt(obs_df, value_vars=keys)
        x = 'variable'
        ys = ['value']
    else:
        obs_tidy = obs_df
        x = groupby
        ys = keys

    # set by default the violin plot cut=0 to limit the extend
    # of the violin plot (see stacked_violin code) for more info.
    if 'cut' not in kwds:
        kwds['cut'] = 0

    if multi_panel and groupby is None and len(ys) == 1:
        # This is a quick and dirty way for adapting scales across several
        # keys if groupby is None.
        y = ys[0]
        g = sns.FacetGrid(obs_tidy, col=x, col_order=keys, sharey=False)
        # don't really know why this gives a warning without passing `order`
        g = g.map(
            sns.violinplot,
            y,
            inner=None,
            orient='vertical',
            scale=scale,
            order=keys,
            **kwds,
        )
        if stripplot:
            g = g.map(
                sns.stripplot,
                y,
                orient='vertical',
                jitter=jitter,
                size=size,
                order=keys,
                color='black',
            )
        if log:
            g.set(yscale='log')
        g.set_titles(col_template='{col_name}').set_xlabels('')
        if rotation is not None:
            for ax in g.axes[0]:
                ax.tick_params(axis='x', labelrotation=rotation)
    else:
        if ax is None:
            axs, _, _, _ = setup_axes(
                ax=ax,
                panels=['x'] if groupby is None else keys,
                show_ticks=True,
                right_margin=0.3,
            )
        else:
            axs = [ax]
        for ax, y in zip(axs, ys):
            ax = sns.violinplot(
                x,
                y=y,
                data=obs_tidy,
                inner=None,
                order=order,
                orient='vertical',
                scale=scale,
                ax=ax,
                **kwds,
            )
            if stripplot:
                ax = sns.stripplot(
                    x,
                    y=y,
                    data=obs_tidy,
                    order=order,
                    jitter=jitter,
                    color='black',
                    size=size,
                    ax=ax,
                )
            if xlabel == '' and groupby is not None and rotation is None:
                xlabel = groupby.replace('_', ' ')
            ax.set_xlabel(xlabel)
            if log:
                ax.set_yscale('log')
            if rotation is not None:
                ax.tick_params(axis='x', labelrotation=rotation)
    show = settings.autoshow if show is None else show
    _utils.savefig_or_show('violin', show=show, save=save)
    if not show:
        if multi_panel and groupby is None and len(ys) == 1:
            return g
        elif len(axs) == 1:
            return axs[0]
        else:
            return axs


@_doc_params(show_save_ax=doc_show_save_ax)
def clustermap(
    adata: AnnData,
    obs_keys: str = None,
    use_raw: Optional[bool] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    **kwds,
):
    """\
    Hierarchically-clustered heatmap.

    Wraps :func:`seaborn.clustermap` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    obs_keys
        Categorical annotation to plot with a different color map.
        Currently, only a single key is supported.
    use_raw
        Use `raw` attribute of `adata` if present.
    {show_save_ax}
    **kwds
        Keyword arguments passed to :func:`~seaborn.clustermap`.

    Returns
    -------
    If `show` is `False`, a :class:`~seaborn.ClusterGrid` object
    (see :func:`~seaborn.clustermap`).

    Examples
    --------
    Soon to come with figures. In the meanwile, see :func:`~seaborn.clustermap`.

    >>> import scanpy as sc
    >>> adata = sc.datasets.krumsiek11()
    >>> sc.pl.clustermap(adata, obs_keys='cell_type')
    """
    import seaborn as sns  # Slow import, only import if called

    if not isinstance(obs_keys, (str, type(None))):
        raise ValueError('Currently, only a single key is supported.')
    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None:
        use_raw = True
    X = adata.raw.X if use_raw else adata.X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    if obs_keys is not None:
        row_colors = adata.obs[obs_keys]
        _utils.add_colors_for_categorical_sample_annotation(adata, obs_keys)
        # do this more efficiently... just a quick solution
        lut = dict(
            zip(row_colors.cat.categories, adata.uns[obs_keys + '_colors'])
        )
        row_colors = adata.obs[obs_keys].map(lut)
        g = sns.clustermap(df, row_colors=row_colors.values, **kwds)
    else:
        g = sns.clustermap(df, **kwds)
    show = settings.autoshow if show is None else show
    _utils.savefig_or_show('clustermap', show=show, save=save)
    if show:
        pl.show()
    else:
        return g


@_doc_params(
    show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args
)
def stacked_violin(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    log: bool = False,
    use_raw: Optional[bool] = None,
    num_categories: int = 7,
    figsize: Optional[Tuple[float, float]] = None,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    standard_scale: Optional[Literal['var', 'obs']] = None,
    var_group_rotation: Optional[float] = None,
    layer: Optional[str] = None,
    stripplot: bool = False,
    jitter: Union[float, bool] = False,
    size: int = 1,
    scale: Literal['area', 'count', 'width'] = 'width',
    order: Optional[Sequence[str]] = None,
    swap_axes: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    row_palette: str = 'muted',
    **kwds,
):
    """\
    Stacked violin plots.

    Makes a compact image composed of individual violin plots
    (from :func:`~seaborn.violinplot`) stacked on top of each other.
    Useful to visualize gene expression per cluster.

    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    {common_plot_args}
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    order
        Order in which to show the categories. Note: if `dendrogram=True`
        the categories order will be given by the dendrogram and `order`
        will be ignored.
    scale
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    row_palette
        The row palette determines the colors to use for the stacked violins.
        The value should be a valid seaborn or matplotlib palette name
        (see :func:`~seaborn.color_palette`).
        Alternatively, a single color name or hex value can be passed,
        e.g. `'red'` or `'#cc33ff'`.
    standard_scale
        Whether or not to standardize a dimension between 0 and 1,
        meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby` categories.
         By setting `swap_axes` then x are the `groupby` categories and y the `var_names`. When swapping
         axes var_group_positions are no longer used
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.violinplot`.

    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`

    Examples
    -------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)

    See also
    --------
    rank_genes_groups_stacked_violin: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """
    import seaborn as sns  # Slow import, only import if called

    if use_raw is None and adata.raw is not None:
        use_raw = True
    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )
    has_var_groups = (
        True
        if var_group_positions is not None and len(var_group_positions) > 0
        else False
    )
    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    if standard_scale == 'obs':
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    if 'color' in kwds:
        row_palette = kwds['color']
        # remove color from kwds in case is set to avoid an error caused by
        # double parameters
        del kwds['color']
    if 'linewidth' not in kwds:
        # for the tiny violin plots used, is best
        # to use a thin lindwidth.
        kwds['linewidth'] = 0.5

    # set by default the violin plot cut=0 to limit the extend
    # of the violin plot as this produces better plots that wont extend
    # to negative values for example. From seaborn.violin documentation:
    #
    # cut: Distance, in units of bandwidth size, to extend the density past
    # the extreme datapoints. Set to 0 to limit the violin range within
    # the range of the observed data (i.e., to have the same effect as
    # trim=True in ggplot.
    if 'cut' not in kwds:
        kwds['cut'] = 0
    if groupby is None or len(categories) <= 1:
        # dendrogram can only be computed  between groupby categories
        dendrogram = False

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
        )

        var_group_labels = dendro_data['var_group_labels']
        var_group_positions = dendro_data['var_group_positions']

        # reorder obs_tidy
        if dendro_data['var_names_idx_ordered'] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data['var_names_idx_ordered']]
            var_names = [
                var_names[x] for x in dendro_data['var_names_idx_ordered']
            ]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data['categories_idx_ordered']],
            ordered=True,
        )
        categories = [
            categories[x] for x in dendro_data['categories_idx_ordered']
        ]

        # set order to None, to avoid a different order in case this is given
        # by the user as the dendrogram order takes precedence.
        if order is not None:
            logg.warning(
                "The `order` and `dendrogram` parameters are both set. "
                "Categories will only be ordered by the order"
                "given by the dendrogram"
            )

    elif order is not None:
        if set(obs_tidy.index.categories) != set(order):
            logg.error(
                "Please check that the categories given by "
                "the `order` parameter match the categories that "
                "want to be reordered.\n\n"
                f"Mismatch: {set(obs_tidy.index.categories).difference(order)}\n\n"
                f"Given order categories: {order}\n\n"
                f"{groupby} categories: {list(obs_tidy.index.categories)}\n"
            )
            return

        else:
            obs_tidy.index = obs_tidy.index.reorder_categories(
                order, ordered=True
            )
    global count
    count = 0

    def rename_cols_to_int(value):
        global count
        count += 1
        return count

    # All columns should have a unique name, otherwise the
    # pd.melt object that is passed to seaborn will merge non-unique columns.
    # Here, I simply rename the columns using a count from 1..n using the
    # mapping function `rename_cols_to_int` to solve the problem.
    obs_tidy.rename(rename_cols_to_int, axis='columns', inplace=True)

    if not swap_axes:
        # plot image in which x = var_names and y = groupby categories

        dendro_width = 1.4 if dendrogram else 0

        if figsize is None:
            height = len(categories) * 0.2 + 3
            width = len(var_names) * 0.2 + 1 + dendro_width
        else:
            width, height = figsize

        num_rows = len(categories)
        height_ratios = None
        if has_var_groups:
            # add some space in case 'brackets' want to be plotted on top of the image
            num_rows += 2  # +2 to add the row for the brackets and a spacer
            height_ratios = [0.2, 0.05] + [
                float(height) / len(categories)
            ] * len(categories)
            categories = [None, None] + list(categories)
        fig = pl.figure(figsize=(width, height))

        # define a layout of nrows = len(categories) rows x 2 columns
        # each row is one violin plot. Second column is reserved for dendrogram (if any)
        # if var_group_positions is defined, a new row is added
        axs = gridspec.GridSpec(
            nrows=num_rows,
            ncols=2,
            height_ratios=height_ratios,
            width_ratios=[width, dendro_width],
            wspace=0.08,
        )
        axs_list = []
        if dendrogram:
            first_plot_idx = 1 if has_var_groups else 0
            dendro_ax = fig.add_subplot(axs[first_plot_idx:, 1])
            _plot_dendrogram(
                dendro_ax, adata, groupby, dendrogram_key=dendrogram
            )
            axs_list.append(dendro_ax)
        ax0 = None
        if is_color_like(row_palette):
            row_colors = [row_palette] * len(categories)
        else:
            row_colors = sns.color_palette(
                row_palette, n_colors=len(categories)
            )
        # Iterate in reverse to start on the bottom plot.
        # This facilitates adding the brackets plot (if
        # needed) by sharing the x axis with a previous
        # violin plot.
        for idx in range(num_rows)[::-1]:
            category = categories[idx]
            if has_var_groups and idx <= 1:
                # if var_group_positions is given, axs[0] and axs[1] are the location for the
                # brackets and a spacer (axs[1])
                if idx == 0:
                    brackets_ax = fig.add_subplot(axs[0], sharex=ax0)
                    _plot_gene_groups_brackets(
                        brackets_ax,
                        group_positions=var_group_positions,
                        group_labels=var_group_labels,
                        rotation=var_group_rotation,
                    )

                continue

            df = pd.melt(
                obs_tidy[obs_tidy.index == category],
                value_vars=obs_tidy.columns,
            )
            if ax0 is None:
                ax = fig.add_subplot(axs[idx, 0])
                ax0 = ax

            else:
                ax = fig.add_subplot(axs[idx, 0])

            axs_list.append(ax)

            ax = sns.violinplot(
                'variable',
                y='value',
                data=df,
                inner=None,
                orient='vertical',
                scale=scale,
                ax=ax,
                color=row_colors[idx],
                **kwds,
            )

            if stripplot:
                ax = sns.stripplot(
                    'variable',
                    y='value',
                    data=df,
                    jitter=jitter,
                    color='black',
                    size=size,
                    ax=ax,
                )

            # remove the grids because in such a compact plot are unnecessary
            ax.grid(False)

            ax.tick_params(
                axis='y',
                left=False,
                right=True,
                labelright=True,
                labelleft=False,
                labelsize='x-small',
                length=1,
                pad=1,
            )

            ax.set_ylabel(
                category,
                rotation=0,
                fontsize='small',
                labelpad=8,
                ha='right',
                va='center',
            )
            ax.set_xlabel('')
            if log:
                ax.set_yscale('log')

            if idx < num_rows - 1:
                # remove the xticks labels except for the last processed plot (first from bottom-up).
                # Because the plots share the x axis it is redundant and less compact to plot the
                # axis ticks and labels for each plot
                ax.set_xticklabels([])
                ax.tick_params(
                    axis='x',
                    bottom=False,
                    top=False,
                    labeltop=False,
                    labelbottom=False,
                )
            else:
                ax.set_xticklabels(var_names)

        ax0.tick_params(axis='x', labelrotation=90, labelsize='small')

    else:
        # plot image in which x = group by and y = var_names
        dendro_height = 3 if dendrogram else 0
        vargroups_width = 0.45 if has_var_groups else 0
        if figsize is None:
            height = len(var_names) * 0.3 + dendro_height
            width = len(categories) * 0.4 + vargroups_width
        else:
            width, height = figsize

        fig = pl.figure(figsize=(width, height))

        # define a layout of nrows = var_names x 1 columns
        # if plot dendrogram a row is added
        # each row is one violin plot.
        num_rows = len(var_names) + 1  # +1 to account for dendrogram
        height_ratios = [dendro_height] + ([1] * len(var_names))

        axs = gridspec.GridSpec(
            nrows=num_rows,
            ncols=2,
            height_ratios=height_ratios,
            wspace=0.2,
            width_ratios=[width - vargroups_width, vargroups_width],
        )

        axs_list = []
        if dendrogram:
            dendro_ax = fig.add_subplot(axs[0])
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                orientation='top',
                dendrogram_key=dendrogram,
            )
            axs_list.append(dendro_ax)
        first_ax = None
        if is_color_like(row_palette):
            row_colors = [row_palette] * len(var_names)
        else:
            row_colors = sns.color_palette(row_palette, n_colors=len(var_names))
        for idx, y in enumerate(obs_tidy.columns):
            ax_idx = idx + 1  # +1 to account that idx 0 is the dendrogram
            if first_ax is None:
                ax = fig.add_subplot(axs[ax_idx, 0])
                first_ax = ax
            else:
                ax = fig.add_subplot(axs[ax_idx, 0])
            axs_list.append(ax)
            ax = sns.violinplot(
                x=obs_tidy.index,
                y=y,
                data=obs_tidy,
                inner=None,
                orient='vertical',
                scale=scale,
                ax=ax,
                color=row_colors[idx],
                **kwds,
            )
            if stripplot:
                ax = sns.stripplot(
                    x=obs_tidy.index,
                    y=y,
                    data=obs_tidy,
                    jitter=jitter,
                    color='black',
                    size=size,
                    ax=ax,
                )

            ax.set_ylabel(
                var_names[idx],
                rotation=0,
                fontsize='small',
                labelpad=8,
                ha='right',
                va='center',
            )
            # remove the grids because in such a compact plot are unnecessary
            ax.grid(False)
            ax.tick_params(
                axis='y',
                right=True,
                labelright=True,
                left=False,
                labelleft=False,
                labelrotation=0,
                labelsize='x-small',
            )
            ax.tick_params(axis='x', labelsize='small')

            # remove the xticks labels except for the last processed plot (first from bottom-up).
            # Because the plots share the x axis it is redundant and less compact to plot the
            # axis for each plot
            if idx < len(var_names) - 1:
                ax.tick_params(
                    labelbottom=False, labeltop=False, bottom=False, top=False
                )
                ax.set_xlabel('')
            if log:
                ax.set_yscale('log')

            if max([len(x) for x in categories]) > 1:
                ax.tick_params(axis='x', labelrotation=90)

        if has_var_groups:
            start = 1 if dendrogram else 0
            gene_groups_ax = fig.add_subplot(axs[start:, 1])
            arr = []
            for idx, pos in enumerate(var_group_positions):
                arr += [idx] * (pos[1] + 1 - pos[0])
            _plot_gene_groups_brackets(
                gene_groups_ax,
                var_group_positions,
                var_group_labels,
                left_adjustment=0.3,
                right_adjustment=0.7,
                orientation='right',
            )
            gene_groups_ax.set_ylim(len(var_names), 0)
            axs_list.append(gene_groups_ax)

    # remove the spacing between subplots
    pl.subplots_adjust(wspace=0, hspace=0)

    _utils.savefig_or_show('stacked_violin', show=show, save=save)

    return axs_list


@_doc_params(
    show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args
)
def heatmap(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    var_group_rotation: Optional[float] = None,
    layer: Optional[str] = None,
    standard_scale: Optional[Literal['var', 'obs']] = None,
    swap_axes: bool = False,
    show_gene_labels: Optional[bool] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwds,
):
    """\
    Heatmap of the expression values of genes.

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
    **kwds
        Are passed to :func:`matplotlib.pyplot.imshow`.

    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`

    Examples
    -------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.heatmap(adata, markers, groupby='bulk_labels', dendrogram=True, swap_axes=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.heatmap(adata, markers, groupby='bulk_labels', dendrogram=True)

    See also
    --------
    rank_genes_groups_heatmap: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """
    if use_raw is None and adata.raw is not None:
        use_raw = True

    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    if standard_scale == 'obs':
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    if groupby is None or len(categories) <= 1:
        categorical = False
        # dendrogram can only be computed  between groupby categories
        dendrogram = False
    else:
        categorical = True
        # get categories colors:
        if groupby + "_colors" in adata.uns:
            groupby_colors = adata.uns[groupby + "_colors"]
        else:
            groupby_colors = None

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
        )

        var_group_labels = dendro_data['var_group_labels']
        var_group_positions = dendro_data['var_group_positions']

        # reorder obs_tidy
        if dendro_data['var_names_idx_ordered'] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data['var_names_idx_ordered']]
            var_names = [
                var_names[x] for x in dendro_data['var_names_idx_ordered']
            ]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data['categories_idx_ordered']],
            ordered=True,
        )

        # reorder groupby colors
        if groupby_colors is not None:
            groupby_colors = [
                groupby_colors[x] for x in dendro_data['categories_idx_ordered']
            ]

    if show_gene_labels is None:
        if len(var_names) <= 50:
            show_gene_labels = True
        else:
            show_gene_labels = False
            logg.warning(
                'Gene labels are not shown when more than 50 genes are visualized. '
                'To show gene labels set `show_gene_labels=True`'
            )
    if categorical:
        obs_tidy = obs_tidy.sort_index()

    colorbar_width = 0.2

    if not swap_axes:
        # define a layout of 2 rows x 4 columns
        # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
        # second row is for main content. This second row is divided into three axes:
        #   first ax is for the categories defined by `groupby`
        #   second ax is for the heatmap
        #   third ax is for the dendrogram
        #   fourth ax is for colorbar

        dendro_width = 1 if dendrogram else 0
        groupby_width = 0.2 if categorical else 0
        if figsize is None:
            height = 6
            if show_gene_labels:
                heatmap_width = len(var_names) * 0.3
            else:
                heatmap_width = 8
            width = heatmap_width + dendro_width + groupby_width
        else:
            width, height = figsize
            heatmap_width = width - (dendro_width + groupby_width)

        if var_group_positions is not None and len(var_group_positions) > 0:
            # add some space in case 'brackets' want to be plotted on top of the image
            height_ratios = [0.15, height]
        else:
            height_ratios = [0, height]

        width_ratios = [
            groupby_width,
            heatmap_width,
            dendro_width,
            colorbar_width,
        ]
        fig = pl.figure(figsize=(width, height))

        axs = gridspec.GridSpec(
            nrows=2,
            ncols=4,
            width_ratios=width_ratios,
            wspace=0.15 / width,
            hspace=0.13 / height,
            height_ratios=height_ratios,
        )

        heatmap_ax = fig.add_subplot(axs[1, 1])
        im = heatmap_ax.imshow(obs_tidy.values, aspect='auto', **kwds)
        heatmap_ax.set_ylim(obs_tidy.shape[0] - 0.5, -0.5)
        heatmap_ax.set_xlim(-0.5, obs_tidy.shape[1] - 0.5)
        heatmap_ax.tick_params(axis='y', left=False, labelleft=False)
        heatmap_ax.set_ylabel('')
        heatmap_ax.grid(False)

        # sns.heatmap(obs_tidy, yticklabels="auto", ax=heatmap_ax, cbar_ax=heatmap_cbar_ax, **kwds)
        if show_gene_labels:
            heatmap_ax.tick_params(axis='x', labelsize='small')
            heatmap_ax.set_xticks(np.arange(len(var_names)))
            heatmap_ax.set_xticklabels(var_names, rotation=90)
        else:
            heatmap_ax.tick_params(axis='x', labelbottom=False, bottom=False)

        # plot colorbar
        _plot_colorbar(im, fig, axs[1, 3])

        if categorical:
            groupby_ax = fig.add_subplot(axs[1, 0])
            ticks, labels, groupby_cmap, norm = _plot_categories_as_colorblocks(
                groupby_ax, obs_tidy, colors=groupby_colors, orientation='left'
            )

            # add lines to main heatmap
            line_positions = (
                np.cumsum(obs_tidy.index.value_counts(sort=False))[:-1] - 0.5
            )
            heatmap_ax.hlines(
                line_positions,
                -0.73,
                len(var_names) - 0.5,
                lw=0.6,
                zorder=10,
                clip_on=False,
            )

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[1, 2], sharey=heatmap_ax)
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                ticks=ticks,
                dendrogram_key=dendrogram,
            )

        # plot group legends on top of heatmap_ax (if given)
        if var_group_positions is not None and len(var_group_positions) > 0:
            gene_groups_ax = fig.add_subplot(axs[0, 1], sharex=heatmap_ax)
            _plot_gene_groups_brackets(
                gene_groups_ax,
                group_positions=var_group_positions,
                group_labels=var_group_labels,
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
            if show_gene_labels:
                heatmap_height = len(var_names) * 0.18
            else:
                heatmap_height = 4
            width = 10
            height = heatmap_height + dendro_height + groupby_height
        else:
            width, height = figsize
            heatmap_height = height - (dendro_height + groupby_height)

        height_ratios = [dendro_height, heatmap_height, groupby_height]

        if var_group_positions is not None and len(var_group_positions) > 0:
            # add some space in case 'brackets' want to be plotted on top of the image
            width_ratios = [width, 0.14, colorbar_width]
        else:
            width_ratios = [width, 0, colorbar_width]

        fig = pl.figure(figsize=(width, height))
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

        im = heatmap_ax.imshow(obs_tidy.T.values, aspect='auto', **kwds)
        heatmap_ax.set_xlim(0, obs_tidy.shape[0])
        heatmap_ax.set_ylim(obs_tidy.shape[1] - 0.5, -0.5)
        heatmap_ax.tick_params(axis='x', bottom=False, labelbottom=False)
        heatmap_ax.set_xlabel('')
        heatmap_ax.grid(False)
        if show_gene_labels:
            heatmap_ax.tick_params(axis='y', labelsize='small', length=1)
            heatmap_ax.set_yticks(np.arange(len(var_names)))
            heatmap_ax.set_yticklabels(var_names, rotation=0)
        else:
            heatmap_ax.tick_params(axis='y', labelleft=False, left=False)

        if categorical:
            groupby_ax = fig.add_subplot(axs[2, 0])
            ticks, labels, groupby_cmap, norm = _plot_categories_as_colorblocks(
                groupby_ax,
                obs_tidy,
                colors=groupby_colors,
                orientation='bottom',
            )
            # add lines to main heatmap
            line_positions = (
                np.cumsum(obs_tidy.index.value_counts(sort=False))[:-1] - 0.5
            )
            heatmap_ax.vlines(
                line_positions,
                -0.5,
                len(var_names) + 0.35,
                lw=0.6,
                zorder=10,
                clip_on=False,
            )

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[0, 0], sharex=heatmap_ax)
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                dendrogram_key=dendrogram,
                ticks=ticks,
                orientation='top',
            )

        # plot group legends next to the heatmap_ax (if given)
        if var_group_positions is not None and len(var_group_positions) > 0:
            gene_groups_ax = fig.add_subplot(axs[1, 1])
            arr = []
            for idx, pos in enumerate(var_group_positions):
                arr += [idx] * (pos[1] + 1 - pos[0])

            gene_groups_ax.imshow(
                np.matrix(arr).T, aspect='auto', cmap=groupby_cmap, norm=norm
            )
            gene_groups_ax.axis('off')

        # plot colorbar
        _plot_colorbar(im, fig, axs[1, 2])

    _utils.savefig_or_show('heatmap', show=show, save=save)

    return axs


@_doc_params(
    show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args
)
def dotplot(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    expression_cutoff: float = 0.0,
    mean_only_expressed: bool = False,
    color_map: str = 'Reds',
    dot_max: Optional[float] = None,
    dot_min: Optional[float] = None,
    standard_scale: Literal['var', 'group'] = None,
    smallest_dot: float = 0.0,
    figsize: Optional[Tuple[float, float]] = None,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    var_group_rotation: Optional[float] = None,
    layer: Optional[str] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    **kwds,
):
    """\
    Makes a *dot plot* of the expression values of `var_names`.

    For each var_name and each `groupby` category a dot is plotted.
    Each dot represents two values: mean expression within each category
    (visualized by color) and fraction of cells expressing the `var_name` in the
    category (visualized by the size of the dot). If `groupby` is not given,
    the dotplot assumes that all data belongs to a single category.

    .. note::
       A gene is considered expressed if the expression value in the `adata` (or
       `adata.raw`) is above the specified threshold which is zero by default.

    An example of dotplot usage is to visualize, for multiple marker genes,
    the mean value and the percentage of cells expressing the gene
    accross multiple clusters.

    Parameters
    ----------
    {common_plot_args}
    expression_cutoff
        Expression cutoff that is used for binarizing the gene expression and
        determining the fraction of cells expressing given genes. A gene is
        expressed only if the expression value is greater than this threshold.
    mean_only_expressed
        If True, gene expression is averaged only over the cells
        expressing the given genes.
    color_map
        String denoting matplotlib color map.
    dot_max
        If none, the maximum dot size is set to the maximum fraction value found
        (e.g. 0.6). If given, the value should be a number between 0 and 1.
        All fractions larger than dot_max are clipped to this value.
    dot_min
        If none, the minimum dot size is set to 0. If given,
        the value should be a number between 0 and 1.
        All fractions smaller than dot_min are clipped to this value.
    standard_scale
        Whether or not to standardize that dimension between 0 and 1,
        meaning for each variable or group,
        subtract the minimum and divide each by its maximum.
    smallest_dot
        If none, the smallest dot has size 0.
        All expression levels with `dot_min` are plotted with this size.

    {show_save_ax}
    **kwds
        Are passed to :func:`matplotlib.pyplot.scatter`.

    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`

    Examples
    -------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    See also
    --------
    rank_genes_groups_dotplot: to plot marker genes identified using the
    :func:`~scanpy.tl.rank_genes_groups` function.
    """
    if use_raw is None and adata.raw is not None:
        use_raw = True
    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )
    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        num_categories,
        layer=layer,
        gene_symbols=gene_symbols,
    )

    # for if category defined by groupby (if any) compute for each var_name
    # 1. the fraction of cells in the category having a value >expression_cutoff
    # 2. the mean value over the category

    # 1. compute fraction of cells having value > expression_cutoff
    # transform obs_tidy into boolean matrix using the expression_cutoff
    obs_bool = obs_tidy > expression_cutoff

    # compute the sum per group which in the boolean matrix this is the number
    # of values >expression_cutoff, and divide the result by the total number of
    # values in the group (given by `count()`)
    fraction_obs = (
        obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
    )

    # 2. compute mean value
    if mean_only_expressed:
        mean_obs = obs_tidy.mask(~obs_bool).groupby(level=0).mean().fillna(0)
    else:
        mean_obs = obs_tidy.groupby(level=0).mean()

    if standard_scale == 'group':
        mean_obs = mean_obs.sub(mean_obs.min(1), axis=0)
        mean_obs = mean_obs.div(mean_obs.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        mean_obs -= mean_obs.min(0)
        mean_obs = (mean_obs / mean_obs.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    dendro_width = 0.8 if dendrogram else 0
    colorbar_width = 0.2
    colorbar_width_spacer = 0.5
    size_legend_width = 0.25
    if figsize is None:
        height = len(categories) * 0.3 + 1  # +1 for labels
        # if the number of categories is small (eg 1 or 2) use
        # a larger height
        height = max([1.5, height])
        heatmap_width = len(var_names) * 0.35
        width = (
            heatmap_width
            + colorbar_width
            + size_legend_width
            + dendro_width
            + colorbar_width_spacer
        )
    else:
        width, height = figsize
        heatmap_width = width - (
            colorbar_width
            + size_legend_width
            + dendro_width
            + colorbar_width_spacer
        )

    # colorbar ax width should not change with differences in the width of the image
    # otherwise can become too small

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.5, 10]
    else:
        height_ratios = [0, 10.5]

    # define a layout of 2 rows x 5 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
    # second row is for main content. This second row
    # is divided into 4 axes:
    #   first ax is for the main figure
    #   second ax is for dendrogram (if present)
    #   third ax is for the color bar legend
    #   fourth ax is for an spacer that avoids the ticks
    #             from the color bar to be hidden beneath the size lengend axis
    #   fifth ax is to plot the size legend
    fig = pl.figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        nrows=2,
        ncols=5,
        wspace=0.02,
        hspace=0.04,
        width_ratios=[
            heatmap_width,
            dendro_width,
            colorbar_width,
            colorbar_width_spacer,
            size_legend_width,
        ],
        height_ratios=height_ratios,
    )
    if len(categories) < 4:
        # when few categories are shown, the colorbar and size legend
        # need to be larger than the main plot, otherwise they would look
        # compressed. For this, the dotplot ax is split into two:
        axs2 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=axs[1, 0],
            height_ratios=[len(categories) * 0.3, 1],
        )
        dot_ax = fig.add_subplot(axs2[0])
    else:
        dot_ax = fig.add_subplot(axs[1, 0])

    color_legend = fig.add_subplot(axs[1, 2])

    if groupby is None or len(categories) <= 1:
        # dendrogram can only be computed  between groupby categories
        dendrogram = False

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
        )

        var_group_labels = dendro_data['var_group_labels']
        var_group_positions = dendro_data['var_group_positions']

        # reorder matrix
        if dendro_data['var_names_idx_ordered'] is not None:
            # reorder columns (usually genes) if needed. This only happens when
            # var_group_positions and var_group_labels is set
            mean_obs = mean_obs.iloc[:, dendro_data['var_names_idx_ordered']]
            fraction_obs = fraction_obs.iloc[
                :, dendro_data['var_names_idx_ordered']
            ]

        # reorder rows (categories) to match the dendrogram order
        mean_obs = mean_obs.iloc[dendro_data['categories_idx_ordered'], :]
        fraction_obs = fraction_obs.iloc[
            dendro_data['categories_idx_ordered'], :
        ]

        y_ticks = range(mean_obs.shape[0])
        dendro_ax = fig.add_subplot(axs[1, 1], sharey=dot_ax)
        _plot_dendrogram(
            dendro_ax, adata, groupby, dendrogram_key=dendrogram, ticks=y_ticks
        )

    # to keep the size_legen of about the same height, irrespective
    # of the number of categories, the fourth ax is subdivided into two parts
    size_legend_height = min(1.3, height)
    # wspace is proportional to the width but a constant value is
    # needed such that the spacing is the same for thinner or wider images.
    wspace = 10.5 / width
    axs3 = gridspec.GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=axs[1, 4],
        wspace=wspace,
        height_ratios=[
            size_legend_height / height,
            (height - size_legend_height) / height,
        ],
    )
    # make scatter plot in which
    # x = var_names
    # y = groupby category
    # size = fraction
    # color = mean expression

    y, x = np.indices(mean_obs.shape)
    y = y.flatten()
    x = x.flatten()
    frac = fraction_obs.values.flatten()
    mean_flat = mean_obs.values.flatten()
    cmap = pl.get_cmap(color_map)
    if dot_max is None:
        dot_max = np.ceil(max(frac) * 10) / 10
    else:
        if dot_max < 0 or dot_max > 1:
            raise ValueError("`dot_max` value has to be between 0 and 1")
    if dot_min is None:
        dot_min = 0
    else:
        if dot_min < 0 or dot_min > 1:
            raise ValueError("`dot_min` value has to be between 0 and 1")

    if dot_min != 0 or dot_max != 1:
        # clip frac between dot_min and  dot_max
        frac = np.clip(frac, dot_min, dot_max)
        old_range = dot_max - dot_min
        # re-scale frac between 0 and 1
        frac = (frac - dot_min) / old_range

    size = (frac * 10) ** 2
    size += smallest_dot
    import matplotlib.colors

    normalize = matplotlib.colors.Normalize(
        vmin=kwds.get('vmin'), vmax=kwds.get('vmax')
    )
    colors = cmap(normalize(mean_flat))
    dot_ax.scatter(
        x,
        y,
        color=colors,
        s=size,
        cmap=cmap,
        norm=None,
        edgecolor='none',
        **kwds,
    )
    y_ticks = range(mean_obs.shape[0])
    dot_ax.set_yticks(y_ticks)
    dot_ax.set_yticklabels([mean_obs.index[idx] for idx in y_ticks])

    x_ticks = range(mean_obs.shape[1])
    dot_ax.set_xticks(x_ticks)
    dot_ax.set_xticklabels(
        [mean_obs.columns[idx] for idx in x_ticks], rotation=90
    )
    dot_ax.tick_params(axis='both', labelsize='small')
    dot_ax.grid(False)
    dot_ax.set_xlim(-0.5, len(var_names) + 0.5)
    dot_ax.set_ylabel(groupby)

    # to be consistent with the heatmap plot, is better to
    # invert the order of the y-axis, such that the first group is on
    # top
    ymin, ymax = dot_ax.get_ylim()
    dot_ax.set_ylim(ymax + 0.5, ymin - 0.5)

    dot_ax.set_xlim(-1, len(var_names))

    # plot group legends on top of dot_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)
        _plot_gene_groups_brackets(
            gene_groups_ax,
            group_positions=var_group_positions,
            group_labels=var_group_labels,
            rotation=var_group_rotation,
        )

    # plot colorbar
    import matplotlib.colorbar

    matplotlib.colorbar.ColorbarBase(color_legend, cmap=cmap, norm=normalize)

    # for the dot size legend, use step between dot_max and dot_min
    # based on how different they are.
    diff = dot_max - dot_min
    if 0.3 < diff <= 0.6:
        step = 0.1
    elif diff <= 0.3:
        step = 0.05
    else:
        step = 0.2
    # a descending range that is afterwards inverted is used
    # to guarantee that dot_max is in the legend.
    fracs_legends = np.arange(dot_max, dot_min, step * -1)[::-1]
    if dot_min != 0 or dot_max != 1:
        fracs_values = (fracs_legends - dot_min) / old_range
    else:
        fracs_values = fracs_legends
    size = (fracs_values * 10) ** 2
    size += smallest_dot
    color = [
        cmap(normalize(value))
        for value in np.repeat(max(mean_flat) * 0.7, len(size))
    ]

    # plot size bar
    size_legend = fig.add_subplot(axs3[0])

    size_legend.scatter(
        np.repeat(0, len(size)), range(len(size)), s=size, color=color
    )
    size_legend.set_yticks(range(len(size)))
    labels = ["{:.0%}".format(x) for x in fracs_legends]
    if dot_max < 1:
        labels[-1] = ">" + labels[-1]
    size_legend.set_yticklabels(labels)
    size_legend.set_yticklabels(["{:.0%}".format(x) for x in fracs_legends])

    size_legend.tick_params(
        axis='y', left=False, labelleft=False, labelright=True
    )

    # remove x ticks and labels
    size_legend.tick_params(axis='x', bottom=False, labelbottom=False)

    # remove surrounding lines
    size_legend.spines['right'].set_visible(False)
    size_legend.spines['top'].set_visible(False)
    size_legend.spines['left'].set_visible(False)
    size_legend.spines['bottom'].set_visible(False)
    size_legend.grid(False)

    ymin, ymax = size_legend.get_ylim()
    size_legend.set_ylim(ymin, ymax + 0.5)

    _utils.savefig_or_show('dotplot', show=show, save=save)
    return axs


@_doc_params(
    show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args
)
def matrixplot(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    figsize: Optional[Tuple[float, float]] = None,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    var_group_rotation: Optional[float] = None,
    layer: Optional[str] = None,
    standard_scale: Literal['var', 'group'] = None,
    swap_axes: bool = False,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    **kwds,
):
    """\
    Creates a heatmap of the mean expression values per cluster of each var_names
    If groupby is not given, the matrixplot assumes that all data belongs to a single
    category.

    Parameters
    ----------
    {common_plot_args}
    standard_scale
        Whether or not to standardize that dimension between 0 and 1, meaning for each variable or group,
        subtract the minimum and divide each by its maximum.
    {show_save_ax}
    **kwds
        Are passed to :func:`matplotlib.pyplot.pcolor`.

    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.matrixplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.matrixplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    See also
    --------
    rank_genes_groups_matrixplot: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """

    if use_raw is None and adata.raw is not None:
        use_raw = True
    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )
    if groupby is None or len(categories) <= 1:
        # dendrogram can only be computed  between groupby categories
        dendrogram = False

    mean_obs = obs_tidy.groupby(level=0).mean()

    if standard_scale == 'group':
        mean_obs = mean_obs.sub(mean_obs.min(1), axis=0)
        mean_obs = mean_obs.div(mean_obs.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        mean_obs -= mean_obs.min(0)
        mean_obs = (mean_obs / mean_obs.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
        )

        var_group_labels = dendro_data['var_group_labels']
        var_group_positions = dendro_data['var_group_positions']

        # reorder matrix
        if dendro_data['var_names_idx_ordered'] is not None:
            # reorder columns (usually genes) if needed. This only happens when
            # var_group_positions and var_group_labels is set
            mean_obs = mean_obs.iloc[:, dendro_data['var_names_idx_ordered']]

        # reorder rows (categories) to match the dendrogram order
        mean_obs = mean_obs.iloc[dendro_data['categories_idx_ordered'], :]

    colorbar_width = 0.2

    if not swap_axes:
        dendro_width = 0.8 if dendrogram else 0
        if figsize is None:
            height = len(categories) * 0.2 + 1  # +1 for labels
            heatmap_width = len(var_names) * 0.32
            width = heatmap_width + dendro_width + colorbar_width
        else:
            width, height = figsize
            heatmap_width = width - (dendro_width + colorbar_width)

        if var_group_positions is not None and len(var_group_positions) > 0:
            # add some space in case 'brackets' want to be plotted on top of the image
            height_ratios = [0.5, 10]
            height += 0.5
        else:
            height_ratios = [0, 10.5]

        # define a layout of 2 rows x 3 columns
        # first row is for 'brackets' (if no brackets needed, the height of this row is zero)
        # second row is for main content. This second row
        # is divided into three axes:
        #   first ax is for the main matrix figure
        #   second ax is for the dendrogram
        #   third ax is for the color bar legend
        fig = pl.figure(figsize=(width, height))
        axs = gridspec.GridSpec(
            nrows=2,
            ncols=3,
            wspace=0.02,
            hspace=0.04,
            width_ratios=[heatmap_width, dendro_width, colorbar_width],
            height_ratios=height_ratios,
        )

        matrix_ax = fig.add_subplot(axs[1, 0])
        y_ticks = np.arange(mean_obs.shape[0]) + 0.5
        matrix_ax.set_yticks(y_ticks)
        matrix_ax.set_yticklabels(
            [mean_obs.index[idx] for idx in range(mean_obs.shape[0])]
        )

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[1, 1], sharey=matrix_ax)
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                dendrogram_key=dendrogram,
                ticks=y_ticks,
            )

        pc = matrix_ax.pcolor(mean_obs, edgecolor='gray', **kwds)

        # invert y axis to show categories ordered from top to bottom
        matrix_ax.set_ylim(mean_obs.shape[0], 0)

        x_ticks = np.arange(mean_obs.shape[1]) + 0.5
        matrix_ax.set_xticks(x_ticks)
        matrix_ax.set_xticklabels(
            [mean_obs.columns[idx] for idx in range(mean_obs.shape[1])],
            rotation=90,
        )
        matrix_ax.tick_params(axis='both', labelsize='small')
        matrix_ax.grid(False)
        matrix_ax.set_xlim(-0.5, len(var_names) + 0.5)
        matrix_ax.set_ylabel(groupby)
        matrix_ax.set_xlim(0, mean_obs.shape[1])

        # plot group legends on top of matrix_ax (if given)
        if var_group_positions is not None and len(var_group_positions) > 0:
            gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=matrix_ax)
            _plot_gene_groups_brackets(
                gene_groups_ax,
                group_positions=var_group_positions,
                group_labels=var_group_labels,
                rotation=var_group_rotation,
                left_adjustment=0.2,
                right_adjustment=0.8,
            )

        # plot colorbar
        _plot_colorbar(pc, fig, axs[1, 2])
    else:
        dendro_height = 0.5 if dendrogram else 0
        if var_group_positions is not None and len(var_group_positions) > 0:
            # add some space in case 'color blocks' want to be plotted on the right of the image
            vargroups_width = 0.4
        else:
            vargroups_width = 0

        if figsize is None:
            heatmap_height = len(var_names) * 0.2
            height = dendro_height + heatmap_height + 1  # +1 for labels
            heatmap_width = len(categories) * 0.3
            width = heatmap_width + vargroups_width + colorbar_width
        else:
            width, height = figsize
            heatmap_width = width - (vargroups_width + colorbar_width)
            heatmap_height = height - dendro_height

        # define a layout of 2 rows x 3 columns
        # first row is for 'dendrogram' (if no dendrogram is plotted, the height of this row is zero)
        # second row is for main content. This row
        # is divided into three axes:
        #   first ax is for the main matrix figure
        #   second ax is for the groupby categories (eg. brackets)
        #   third ax is for the color bar legend
        fig = pl.figure(figsize=(width, height))
        axs = gridspec.GridSpec(
            nrows=2,
            ncols=3,
            wspace=0.05,
            hspace=0.005,
            width_ratios=[heatmap_width, vargroups_width, colorbar_width],
            height_ratios=[dendro_height, heatmap_height],
        )

        mean_obs = mean_obs.T
        matrix_ax = fig.add_subplot(axs[1, 0])
        pc = matrix_ax.pcolor(mean_obs, edgecolor='gray', **kwds)
        y_ticks = np.arange(mean_obs.shape[0]) + 0.5
        matrix_ax.set_yticks(y_ticks)
        matrix_ax.set_yticklabels(
            [mean_obs.index[idx] for idx in range(mean_obs.shape[0])]
        )

        x_ticks = np.arange(mean_obs.shape[1]) + 0.5
        matrix_ax.set_xticks(x_ticks)
        matrix_ax.set_xticklabels(
            [mean_obs.columns[idx] for idx in range(mean_obs.shape[1])],
            rotation=90,
        )
        matrix_ax.tick_params(axis='both', labelsize='small')
        matrix_ax.grid(False)
        matrix_ax.set_xlim(0, len(categories))
        matrix_ax.set_xlabel(groupby)
        # invert y axis to show var_names ordered from top to bottom
        matrix_ax.set_ylim(mean_obs.shape[0], 0)

        if dendrogram:
            dendro_ax = fig.add_subplot(axs[0, 0], sharex=matrix_ax)
            _plot_dendrogram(
                dendro_ax,
                adata,
                groupby,
                dendrogram_key=dendrogram,
                ticks=x_ticks,
                orientation='top',
            )

        # plot group legends on top of matrix_ax (if given)
        if var_group_positions is not None and len(var_group_positions) > 0:
            gene_groups_ax = fig.add_subplot(axs[1, 1], sharey=matrix_ax)
            _plot_gene_groups_brackets(
                gene_groups_ax,
                group_positions=var_group_positions,
                group_labels=var_group_labels,
                rotation=var_group_rotation,
                left_adjustment=0.2,
                right_adjustment=0.8,
                orientation='right',
            )

        # plot colorbar
        _plot_colorbar(pc, fig, axs[1, 2])

    _utils.savefig_or_show('matrixplot', show=show, save=save)
    return axs


@_doc_params(
    show_save_ax=doc_show_save_ax, common_plot_args=doc_common_plot_args
)
def tracksplot(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: str,
    use_raw: Optional[bool] = None,
    log: bool = False,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    layer: Optional[str] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwds,
):
    """\
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
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.tracksplot(adata, markers, 'bulk_labels', dendrogram=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.heatmap(adata, markers, groupby='bulk_labels', dendrogram=True)

    .. currentmodule:: scanpy

    See also
    --------
    pl.rank_genes_groups_tracksplot: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """

    if (
        groupby not in adata.obs_keys()
        or adata.obs[groupby].dtype.name != 'category'
    ):
        raise ValueError(
            'groupby has to be a valid categorical observation. '
            f'Given value: {groupby}, valid categorical observations: '
            f'{[x for x in adata.obs_keys() if adata.obs[x].dtype.name == "category"]}'
        )

    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        None,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    # get categories colors:
    if groupby + "_colors" not in adata.uns:
        from ._utils import _set_default_colors_for_categorical_obs

        _set_default_colors_for_categorical_obs(adata, groupby)

    groupby_colors = adata.uns[groupby + "_colors"]

    if dendrogram:
        # compute dendrogram if needed and reorder
        # rows and columns to match leaves order.
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
        )
        # reorder obs_tidy
        if dendro_data['var_names_idx_ordered'] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data['var_names_idx_ordered']]
            var_names = [
                var_names[x] for x in dendro_data['var_names_idx_ordered']
            ]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data['categories_idx_ordered']],
            ordered=True,
        )
        categories = [
            categories[x] for x in dendro_data['categories_idx_ordered']
        ]

        groupby_colors = [
            groupby_colors[x] for x in dendro_data['categories_idx_ordered']
        ]

    obs_tidy = obs_tidy.sort_index()

    # obtain the start and end of each category and make
    # a list of ranges that will be used to plot a different
    # color
    cumsum = [0] + list(np.cumsum(obs_tidy.index.value_counts(sort=False)))
    x_values = [(x, y) for x, y in zip(cumsum[:-1], cumsum[1:])]

    dendro_height = 1 if dendrogram else 0

    groupby_height = 0.24
    # +2 because of dendrogram on top and categories at bottom
    num_rows = len(var_names) + 2
    if figsize is None:
        width = 12
        track_height = 0.25
    else:
        width, height = figsize
        track_height = (height - (dendro_height + groupby_height)) / len(
            var_names
        )

    height_ratios = (
        [dendro_height] + [track_height] * len(var_names) + [groupby_height]
    )
    height = sum(height_ratios)

    obs_tidy = obs_tidy.T

    fig = pl.figure(figsize=(width, height))
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
        for cat_idx, category in enumerate(categories):
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
            ax.tick_params(
                labelbottom=False, labeltop=False, bottom=False, top=False
            )
            ax.set_xlabel('')
        if log:
            ax.set_yscale('log')
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.grid(False)
        ymin, ymax = ax.get_ylim()
        ymax = int(ymax)
        ax.set_yticks([ymax])
        tt = ax.set_yticklabels([str(ymax)], ha='left', va='top')
        ax.spines['right'].set_position(('axes', 1.01))
        ax.tick_params(
            axis='y',
            labelsize='x-small',
            right=True,
            left=False,
            length=2,
            which='both',
            labelright=True,
            labelleft=False,
            direction='in',
        )
        ax.set_ylabel(
            var, rotation=0, fontsize='small', ha='right', va='bottom'
        )
        ax.yaxis.set_label_coords(-0.005, 0.1)
    ax.set_xlim(0, x_end)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)

    # the ax to plot the groupby categories is split to add a small space
    # between the rest of the plot and the categories
    axs2 = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=axs[num_rows - 1, 0], height_ratios=[1, 1]
    )

    groupby_ax = fig.add_subplot(axs2[1])

    ticks, labels, groupby_cmap, norm = _plot_categories_as_colorblocks(
        groupby_ax, obs_tidy.T, colors=groupby_colors, orientation='bottom'
    )
    # add lines to plot
    overlay_ax = fig.add_subplot(axs[1:-1, 0], sharex=first_ax)
    line_positions = np.cumsum(obs_tidy.T.index.value_counts(sort=False))[:-1]
    overlay_ax.vlines(line_positions, 0, 1, lw=0.5, linestyle="--")
    overlay_ax.axis('off')
    overlay_ax.set_ylim(0, 1)

    if dendrogram:
        dendro_ax = fig.add_subplot(axs[0], sharex=first_ax)
        _plot_dendrogram(
            dendro_ax,
            adata,
            groupby,
            dendrogram_key=dendrogram,
            orientation='top',
            ticks=ticks,
        )
        axs_list.append(dendro_ax)

    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[1:-1, 1])
        arr = []
        for idx, pos in enumerate(var_group_positions):
            arr += [idx] * (pos[1] + 1 - pos[0])

        gene_groups_ax.imshow(
            np.matrix(arr).T, aspect='auto', cmap=groupby_cmap, norm=norm
        )
        gene_groups_ax.axis('off')
        axs_list.append(gene_groups_ax)

    _utils.savefig_or_show('tracksplot', show=show, save=save)
    return axs_list


@_doc_params(show_save_ax=doc_show_save_ax)
def dendrogram(
    adata: AnnData,
    groupby: str,
    *,
    dendrogram_key: Optional[str] = None,
    orientation: Literal['top', 'bottom', 'left', 'right'] = 'top',
    remove_labels: bool = False,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    ax: Optional[Axes] = None,
):
    """\
    Plots a dendrogram of the categories defined in `groupby`.

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
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, 'bulk_labels')
    >>> sc.pl.dendrogram(adata, 'bulk_labels')
    """
    if ax is None:
        _, ax = pl.subplots()
    _plot_dendrogram(
        ax,
        adata,
        groupby,
        dendrogram_key=dendrogram_key,
        remove_labels=remove_labels,
        orientation=orientation,
    )
    _utils.savefig_or_show('dendrogram', show=show, save=save)
    return ax


@_doc_params(show_save_ax=doc_show_save_ax)
def correlation_matrix(
    adata: AnnData,
    groupby: str,
    show_correlation_numbers: bool = False,
    dendrogram: Union[bool, str, None] = None,
    figsize: Optional[Tuple[float, float]] = None,
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    ax: Optional[Axes] = None,
    **kwds,
) -> Union[Axes, List[Axes]]:
    """\
    Plots the correlation matrix computed as part of `sc.tl.dendrogram`.

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
    **kwds
        Only if `show_correlation` is True:
        Are passed to :func:`matplotlib.pyplot.pcolormesh` when plotting the
        correlation heatmap. Useful values to pas are `vmax`, `vmin` and `cmap`.

    Returns
    -------

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, 'bulk_labels')
    >>> sc.pl.correlation_matrix(adata, 'bulk_labels')
    """

    dendrogram_key = _get_dendrogram_key(adata, dendrogram, groupby)

    index = adata.uns[dendrogram_key]['categories_idx_ordered']
    corr_matrix = adata.uns[dendrogram_key]['correlation_matrix']
    # reorder matrix columns according to the dendrogram
    if dendrogram is None:
        dendrogram = ax is None
    if dendrogram:
        if ax is not None:
            raise ValueError(
                'Can only plot dendrogram when not plotting to an axis'
            )
        assert (len(index)) == corr_matrix.shape[0]
        corr_matrix = corr_matrix[index, :]
        corr_matrix = corr_matrix[:, index]
        labels = list(adata.obs[groupby].cat.categories)
        labels = np.array(labels).astype('str')[index]
    else:
        labels = adata.obs[groupby].cat.categories
    num_rows = corr_matrix.shape[0]
    colorbar_height = 0.2
    if dendrogram:
        dendrogram_width = 1.8
    else:
        dendrogram_width = 0
    if figsize is None:
        corr_matrix_height = num_rows * 0.6
        height = corr_matrix_height + colorbar_height
        width = corr_matrix_height + dendrogram_width
    else:
        width, height = figsize
        corr_matrix_height = height - colorbar_height

    fig = pl.figure(figsize=(width, height)) if ax is None else None
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
            orientation='left',
            ticks=np.arange(corr_matrix.shape[0]) + 0.5,
        )
        axs.append(dendro_ax)
    # define some default pcolormesh parameters
    if 'edge_color' not in kwds:
        if corr_matrix.shape[0] > 30:
            # when there are too many rows it is better to remove
            # the black lines surrounding the boxes in the heatmap
            kwds['edgecolors'] = 'none'
        else:
            kwds['edgecolors'] = 'black'
            kwds['linewidth'] = 0.01
    if 'vmax' not in kwds and 'vmin' not in kwds:
        kwds['vmax'] = 1
        kwds['vmin'] = -1
    if 'cmap' not in kwds:
        # by default use a divergent color map
        kwds['cmap'] = 'bwr'

    img_mat = corr_matrix_ax.pcolormesh(corr_matrix, **kwds)
    corr_matrix_ax.set_xlim(0, num_rows)
    corr_matrix_ax.set_ylim(0, num_rows)

    corr_matrix_ax.yaxis.tick_right()
    corr_matrix_ax.set_yticks(np.arange(corr_matrix.shape[0]) + 0.5)
    corr_matrix_ax.set_yticklabels(labels)

    corr_matrix_ax.xaxis.set_tick_params(labeltop=True)
    corr_matrix_ax.xaxis.set_tick_params(labelbottom=False)
    corr_matrix_ax.set_xticks(np.arange(corr_matrix.shape[0]) + 0.5)
    corr_matrix_ax.set_xticklabels(labels, rotation=45, ha='left')

    for ax_name in 'xy':
        corr_matrix_ax.tick_params(
            axis=ax_name, which='both', bottom=False, top=False
        )

    if show_correlation_numbers:
        for row, col in product(range(num_rows), repeat=2):
            corr_matrix_ax.text(
                row + 0.5,
                col + 0.5,
                f"{corr_matrix[row, col]:.2f}",
                ha='center',
                va='center',
            )

    axs.append(corr_matrix_ax)

    if ax is None:  # Plot colorbar
        colormap_ax = fig.add_subplot(gs[3])
        cobar = pl.colorbar(img_mat, cax=colormap_ax, orientation='horizontal')
        cobar.solids.set_edgecolor("face")
        axs.append(colormap_ax)

    show = settings.autoshow if show is None else show
    _utils.savefig_or_show('correlation_matrix', show=show, save=save)
    if ax is None and not show:
        return axs


def _prepare_dataframe(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    layer=None,
    gene_symbols: Optional[str] = None,
):
    """
    Given the anndata object, prepares a data frame in which the row index are the categories
    defined by group by and the columns correspond to var_names.

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
        Use `raw` attribute of `adata` if present.
    log
        Use the log of the values
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
    from scipy.sparse import issparse

    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None:
        use_raw = True
    if isinstance(var_names, str):
        var_names = [var_names]

    if groupby is not None:
        if groupby not in adata.obs_keys():
            raise ValueError(
                'groupby has to be a valid observation. '
                f'Given {groupby}, valid observations: {adata.obs_keys()}'
            )

    if gene_symbols is not None and gene_symbols in adata.var.columns:
        # translate gene_symbols to var_names
        # slow method but gives a meaningful error if no gene symbol is found:
        translated_var_names = []
        for symbol in var_names:
            if symbol not in adata.var[gene_symbols].values:
                logg.error(
                    f"Gene symbol {symbol!r} not found in given "
                    f"gene_symbols column: {gene_symbols!r}"
                )
                return
            translated_var_names.append(
                adata.var[adata.var[gene_symbols] == symbol].index[0]
            )
        symbols = var_names
        var_names = translated_var_names
    if layer is not None:
        if layer not in adata.layers.keys():
            raise KeyError(
                f'Selected layer: {layer} is not in the layers list. '
                f'The list of valid layers is: {adata.layers.keys()}'
            )
        matrix = adata[:, var_names].layers[layer]
    elif use_raw:
        matrix = adata.raw[:, var_names].X
    else:
        matrix = adata[:, var_names].X

    if issparse(matrix):
        matrix = matrix.toarray()
    if log:
        matrix = np.log1p(matrix)

    obs_tidy = pd.DataFrame(matrix, columns=var_names)
    if groupby is None:
        groupby = ''
        categorical = pd.Series(np.repeat('', len(obs_tidy))).astype('category')
    else:
        if not is_categorical_dtype(adata.obs[groupby]):
            # if the groupby column is not categorical, turn it into one
            # by subdividing into  `num_categories` categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]

    obs_tidy.set_index(categorical, groupby, inplace=True)
    if gene_symbols is not None:
        # translate the column names to the symbol names
        obs_tidy.rename(
            columns=dict(
                [(var_names[x], symbols[x]) for x in range(len(var_names))]
            ),
            inplace=True,
        )
    categories = obs_tidy.index.categories

    return categories, obs_tidy


def _plot_gene_groups_brackets(
    gene_groups_ax: Axes,
    group_positions: Iterable[Tuple[int, int]],
    group_labels: Sequence[str],
    left_adjustment: float = -0.3,
    right_adjustment: float = 0.3,
    rotation: Optional[float] = None,
    orientation: Literal['top', 'right'] = 'top',
):
    """\
    Draws brackets that represent groups of genes on the give axis.
    For best results, this axis is located on top of an image whose
    x axis contains gene names.

    The gene_groups_ax should share the x axis with the main ax.

    Eg: gene_groups_ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)

    This function is used by dotplot, heatmap etc.

    Parameters
    ----------
    gene_groups_ax
        In this axis the gene marks are drawn
    group_positions
        Each item in the list, should contain the start and end position that the
        bracket should cover.
        Eg. [(0, 4), (5, 8)] means that there are two brackets, one for the var_names (eg genes)
        in positions 0-4 and other for positions 5-8
    group_labels
        List of group labels
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
    import matplotlib.patches as patches
    from matplotlib.path import Path

    # get the 'brackets' coordinates as lists of start and end positions

    left = [x[0] + left_adjustment for x in group_positions]
    right = [x[1] + right_adjustment for x in group_positions]

    # verts and codes are used by PathPatch to make the brackets
    verts = []
    codes = []
    if orientation == 'top':
        # rotate labels if any of them is longer than 4 characters
        if rotation is None and group_labels:
            if max([len(x) for x in group_labels]) > 4:
                rotation = 90
            else:
                rotation = 0
        for idx in range(len(left)):
            verts.append((left[idx], 0))  # lower-left
            verts.append((left[idx], 0.6))  # upper-left
            verts.append((right[idx], 0.6))  # upper-right
            verts.append((right[idx], 0))  # lower-right

            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            try:
                group_x_center = left[idx] + float(right[idx] - left[idx]) / 2
                gene_groups_ax.text(
                    group_x_center,
                    1.1,
                    group_labels[idx],
                    ha='center',
                    va='bottom',
                    rotation=rotation,
                )
            except:
                pass
    else:
        top = left
        bottom = right
        for idx in range(len(top)):
            verts.append((0, top[idx]))  # upper-left
            verts.append((0.15, top[idx]))  # upper-right
            verts.append((0.15, bottom[idx]))  # lower-right
            verts.append((0, bottom[idx]))  # lower-left

            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            try:
                diff = bottom[idx] - top[idx]
                group_y_center = top[idx] + float(diff) / 2
                if diff * 2 < len(group_labels[idx]):
                    # cut label to fit available space
                    group_labels[idx] = group_labels[idx][: int(diff * 2)] + "."
                gene_groups_ax.text(
                    0.6,
                    group_y_center,
                    group_labels[idx],
                    ha='right',
                    va='center',
                    rotation=270,
                    fontsize='small',
                )
            except Exception as e:
                print('problems {}'.format(e))
                pass

    path = Path(verts, codes)

    patch = patches.PathPatch(path, facecolor='none', lw=1.5)

    gene_groups_ax.add_patch(patch)
    gene_groups_ax.spines['right'].set_visible(False)
    gene_groups_ax.spines['top'].set_visible(False)
    gene_groups_ax.spines['left'].set_visible(False)
    gene_groups_ax.spines['bottom'].set_visible(False)
    gene_groups_ax.grid(False)

    # remove y ticks
    gene_groups_ax.tick_params(axis='y', left=False, labelleft=False)
    # remove x ticks and labels
    gene_groups_ax.tick_params(
        axis='x', bottom=False, labelbottom=False, labeltop=False
    )


def _reorder_categories_after_dendrogram(
    adata: AnnData,
    groupby,
    dendrogram,
    var_names=None,
    var_group_labels=None,
    var_group_positions=None,
):
    """\
    Function used by plotting functions that need to reorder the the groupby
    observations based on the dendrogram results.

    The function checks if a dendrogram has already been precomputed.
    If not, `sc.tl.dendrogram` is run with default parameters.

    The results found in `.uns[dendrogram_key]` are used to reorder
    `var_group_labels` and `var_group_positions`.


    Returns
    -------
    dictionary with keys:
    'categories_idx_ordered', 'var_group_names_idx_ordered',
    'var_group_labels', and 'var_group_positions'
    """

    key = _get_dendrogram_key(adata, dendrogram, groupby)

    dendro_info = adata.uns[key]
    if groupby != dendro_info['groupby']:
        raise ValueError(
            "Incompatible observations. The precomputed dendrogram contains "
            f"information for the observation: '{groupby}' while the plot is "
            f"made for the observation: '{dendro_info['groupby']}. "
            "Please run `sc.tl.dendrogram` using the right observation.'"
        )

    categories = adata.obs[dendro_info['groupby']].cat.categories

    # order of groupby categories
    categories_idx_ordered = dendro_info['categories_idx_ordered']

    if len(categories) != len(categories_idx_ordered):
        raise ValueError(
            "Incompatible observations. Dendrogram data has "
            f"{len(categories_idx_ordered)} categories but current groupby "
            f"observation {groupby!r} contains {len(categories)} categories. "
            "Most likely the underlying groupby observation changed after the "
            "initial computation of `sc.tl.dendrogram`. "
            "Please run `sc.tl.dendrogram` again.'"
        )

    # reorder var_groups (if any)
    if var_names is not None:
        var_names_idx_ordered = list(range(len(var_names)))

    if var_group_positions:
        if list(var_group_labels) == list(categories):
            positions_ordered = []
            labels_ordered = []
            position_start = 0
            var_names_idx_ordered = []
            for idx in categories_idx_ordered:
                position = var_group_positions[idx]
                _var_names = var_names[position[0] : position[1] + 1]
                var_names_idx_ordered.extend(
                    range(position[0], position[1] + 1)
                )
                positions_ordered.append(
                    (position_start, position_start + len(_var_names) - 1)
                )
                position_start += len(_var_names)
                labels_ordered.append(var_group_labels[idx])
            var_group_labels = labels_ordered
            var_group_positions = positions_ordered
        else:
            logg.warning(
                "Groups are not reordered because the `groupby` categories "
                "and the `var_group_labels` are different.\n"
                f"categories: {_format_first_three_categories(categories)}\n"
                f"var_group_labels: {_format_first_three_categories(var_group_labels)}"
            )
    else:
        var_names_idx_ordered = None

    return dict(
        categories_idx_ordered=categories_idx_ordered,
        var_names_idx_ordered=var_names_idx_ordered,
        var_group_labels=var_group_labels,
        var_group_positions=var_group_positions,
    )


def _format_first_three_categories(categories):
    categories = list(categories)
    if len(categories) > 3:
        categories = categories[:3] + ['etc.']
    return ', '.join(categories)


def _get_dendrogram_key(adata, dendrogram_key, groupby):
    # the `dendrogram_key` can be a bool an NoneType or the name of the
    # dendrogram key. By default the name of the dendrogram key is 'dendrogram'
    if not isinstance(dendrogram_key, str):
        dendrogram_key = f'dendrogram_{groupby}'

    if dendrogram_key not in adata.uns:
        from ..tools._dendrogram import dendrogram

        logg.warning(
            f"dendrogram data not found (using key={dendrogram_key}). "
            "Running `sc.tl.dendrogram` with default parameters. For fine "
            "tuning it is recommended to run `sc.tl.dendrogram` independently."
        )
        dendrogram(adata, groupby, key_added=dendrogram_key)

    if 'dendrogram_info' not in adata.uns[dendrogram_key]:
        raise ValueError(
            f"The given dendrogram key ({dendrogram_key!r}) does not contain "
            "valid dendrogram information."
        )

    return dendrogram_key


def _plot_dendrogram(
    dendro_ax: Axes,
    adata: AnnData,
    groupby: str,
    dendrogram_key: Optional[str] = None,
    orientation: Literal['top', 'bottom', 'left', 'right'] = 'right',
    remove_labels: bool = True,
    ticks: Optional[Collection[float]] = None,
):
    """\
    Plots a dendrogram on the given ax using the precomputed dendrogram
    information stored in `.uns[dendrogram_key]`
    """

    dendrogram_key = _get_dendrogram_key(adata, dendrogram_key, groupby)

    def translate_pos(pos_list, new_ticks, old_ticks):
        """\
        transforms the dendrogram coordinates to a given new position.
        The xlabel_pos and orig_ticks should be of the same
        length.

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
        ...     [0,  1,  2, 3 ],
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

    dendro_info = adata.uns[dendrogram_key]['dendrogram_info']
    leaves = dendro_info["ivl"]
    icoord = np.array(dendro_info['icoord'])
    dcoord = np.array(dendro_info['dcoord'])

    orig_ticks = np.arange(5, len(leaves) * 10 + 5, 10).astype(float)
    # check that ticks has the same length as orig_ticks
    if ticks is not None and len(orig_ticks) != len(ticks):
        logg.warning(
            "ticks argument does not have the same size as orig_ticks. "
            "The argument will be ignored"
        )
        ticks = None

    for xs, ys in zip(icoord, dcoord):
        if ticks is not None:
            xs = translate_pos(xs, ticks, orig_ticks)
        if orientation in ['right', 'left']:
            xs, ys = ys, xs
        dendro_ax.plot(xs, ys, color='#555555')

    dendro_ax.tick_params(bottom=False, top=False, left=False, right=False)
    ticks = ticks if ticks is not None else orig_ticks
    if orientation in ['right', 'left']:
        dendro_ax.set_yticks(ticks)
        dendro_ax.set_yticklabels(leaves, fontsize='small', rotation=0)
        dendro_ax.tick_params(labelbottom=False, labeltop=False)
        if orientation == 'left':
            xmin, xmax = dendro_ax.get_xlim()
            dendro_ax.set_xlim(xmax, xmin)
            dendro_ax.tick_params(labelleft=False, labelright=True)
    else:
        dendro_ax.set_xticks(ticks)
        dendro_ax.set_xticklabels(leaves, fontsize='small', rotation=90)
        dendro_ax.tick_params(labelleft=False, labelright=False)
        if orientation == 'bottom':
            ymin, ymax = dendro_ax.get_ylim()
            dendro_ax.set_ylim(ymax, ymin)
            dendro_ax.tick_params(labeltop=True, labelbottom=False)

    if remove_labels:
        dendro_ax.tick_params(
            labelbottom=False, labeltop=False, labelleft=False, labelright=False
        )

    dendro_ax.grid(False)

    dendro_ax.spines['right'].set_visible(False)
    dendro_ax.spines['top'].set_visible(False)
    dendro_ax.spines['left'].set_visible(False)
    dendro_ax.spines['bottom'].set_visible(False)


def _plot_categories_as_colorblocks(
    groupby_ax: Axes,
    obs_tidy: pd.DataFrame,
    colors=None,
    orientation: Literal['top', 'bottom', 'left', 'right'] = 'left',
    cmap_name: str = 'tab20',
):
    """\
    Plots categories as colored blocks. If orientation is 'left', the categories
    are plotted vertically, otherwise they are plotted horizontally.

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
    from matplotlib.colors import ListedColormap, BoundaryNorm

    if colors is None:
        groupby_cmap = pl.get_cmap(cmap_name)
    else:
        groupby_cmap = ListedColormap(colors, groupby + '_cmap')
    norm = BoundaryNorm(np.arange(groupby_cmap.N + 1) - 0.5, groupby_cmap.N)

    # determine groupby label positions such that they appear
    # centered next/below to the color code rectangle assigned to the category
    value_sum = 0
    ticks = []  # list of centered position of the labels
    labels = []
    label2code = {}  # dictionary of numerical values asigned to each label
    for code, (label, value) in enumerate(
        obs_tidy.index.value_counts(sort=False).iteritems()
    ):
        ticks.append(value_sum + (value / 2))
        labels.append(label)
        value_sum += value
        label2code[label] = code

    groupby_ax.grid(False)

    if orientation == 'left':
        groupby_ax.imshow(
            np.matrix([label2code[lab] for lab in obs_tidy.index]).T,
            aspect='auto',
            cmap=groupby_cmap,
            norm=norm,
        )
        if len(labels) > 1:
            groupby_ax.set_yticks(ticks)
            groupby_ax.set_yticklabels(labels)

        # remove y ticks
        groupby_ax.tick_params(axis='y', left=False, labelsize='small')
        # remove x ticks and labels
        groupby_ax.tick_params(axis='x', bottom=False, labelbottom=False)

        # remove surrounding lines
        groupby_ax.spines['right'].set_visible(False)
        groupby_ax.spines['top'].set_visible(False)
        groupby_ax.spines['left'].set_visible(False)
        groupby_ax.spines['bottom'].set_visible(False)

        groupby_ax.set_ylabel(groupby)
    else:
        groupby_ax.imshow(
            np.matrix([label2code[lab] for lab in obs_tidy.index]),
            aspect='auto',
            cmap=groupby_cmap,
            norm=norm,
        )
        if len(labels) > 1:
            groupby_ax.set_xticks(ticks)
            if max([len(x) for x in labels]) < 3:
                # if the labels are small do not rotate them
                rotation = 0
            else:
                rotation = 90
            groupby_ax.set_xticklabels(labels, rotation=rotation)

        # remove x ticks
        groupby_ax.tick_params(axis='x', bottom=False, labelsize='small')
        # remove y ticks and labels
        groupby_ax.tick_params(axis='y', left=False, labelleft=False)

        # remove surrounding lines
        groupby_ax.spines['right'].set_visible(False)
        groupby_ax.spines['top'].set_visible(False)
        groupby_ax.spines['left'].set_visible(False)
        groupby_ax.spines['bottom'].set_visible(False)

        groupby_ax.set_xlabel(groupby)

    return ticks, labels, groupby_cmap, norm


def _plot_colorbar(mappable, fig, subplot_spec, max_cbar_height: float = 4.0):
    """
    Plots a vertical color bar based on mappable.
    The height of the colorbar is min(figure-height, max_cmap_height)

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
    pl.colorbar(mappable, cax=heatmap_cbar_ax)
    return heatmap_cbar_ax


def _check_var_names_type(var_names, var_group_labels, var_group_positions):
    """
    checks if var_names is a dict. Is this is the cases, then set the
    correct values for var_group_labels and var_group_positions

    Returns
    -------
    var_names, var_group_labels, var_group_positions

    """
    if isinstance(var_names, cabc.Mapping):
        if var_group_labels is not None or var_group_positions is not None:
            logg.warning(
                "`var_names` is a dictionary. This will reset the current "
                "value of `var_group_labels` and `var_group_positions`."
            )
        var_group_labels = []
        _var_names = []
        var_group_positions = []
        start = 0
        for label, vars_list in var_names.items():
            if isinstance(vars_list, str):
                vars_list = [vars_list]
            # use list() in case var_list is a numpy array or pandas series
            _var_names.extend(list(vars_list))
            var_group_labels.append(label)
            var_group_positions.append((start, start + len(vars_list) - 1))
            start += len(vars_list)
        var_names = _var_names

    elif isinstance(var_names, str):
        var_names = [var_names]

    return var_names, var_group_labels, var_group_positions
