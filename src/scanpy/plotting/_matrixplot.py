from __future__ import annotations

from dataclasses import KW_ONLY, InitVar, dataclass
from typing import TYPE_CHECKING, ClassVar, cast

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .. import logging as logg
from .._compat import old_positionals
from .._settings import settings
from .._utils import _doc_params, _empty
from ._baseplot_class import BasePlot, doc_common_groupby_plot_args
from ._docs import (
    doc_common_plot_args,
    doc_show_save_ax,
    doc_vboundnorm,
)
from ._utils import DefaultProxy, check_colornorm, fix_kwds, savefig_or_show

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import Literal, Self

    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.colors import Normalize

    from .._utils import Empty
    from ._baseplot_class import _VarNames
    from ._utils import ColorLike, _AxesSubplot


@_doc_params(common_plot_args=doc_common_plot_args)
@dataclass
class MatrixPlot(BasePlot):
    """\
    Allows the visualization of values using a color map.

    Parameters
    ----------
    {common_plot_args}
    title
        Title for the figure.
    expression_cutoff
        Expression cutoff that is used for binarizing the gene expression and
        determining the fraction of cells expressing given genes. A gene is
        expressed only if the expression value is greater than this threshold.
    mean_only_expressed
        If True, gene expression is averaged only over the cells
        expressing the given genes.
    standard_scale
        Whether or not to standardize that dimension between 0 and 1,
        meaning for each variable or group,
        subtract the minimum and divide each by its maximum.
    values_df
        Optionally, a dataframe with the values to plot can be given. The
        index should be the grouby categories and the columns the genes names.

    kwds
        Are passed to :func:`matplotlib.pyplot.scatter`.

    See also
    --------
    :func:`~scanpy.pl.matrixplot`: Simpler way to call MatrixPlot but with less options.
    :func:`~scanpy.pl.rank_genes_groups_matrixplot`: to plot marker genes identified
        using the :func:`~scanpy.tl.rank_genes_groups` function.

    Examples
    --------

    Simple visualization of the average expression of a few genes grouped by
    the category 'bulk_labels'.

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.MatrixPlot(adata, markers, groupby='bulk_labels').show()

    Same visualization but passing var_names as dict, which adds a grouping of
    the genes on top of the image:

    .. plot::
        :context: close-figs

        markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
        sc.pl.MatrixPlot(adata, markers, groupby='bulk_labels').show()
    """

    DEFAULT_SAVE_PREFIX: ClassVar[str] = "matrixplot_"

    _: KW_ONLY
    colorbar_title: str = "Mean expression\nin group"
    # default style parameters
    cmap = None  # aka: rcParams["image.cmap"]
    values_df: pd.DataFrame | None = None
    standard_scale: InitVar[Literal["var", "group"] | None] = None
    edge_color: ColorLike = "gray"
    edge_lw: float = 0.1

    # deprecated default class variables
    DEFAULT_EDGE_COLOR: ClassVar[DefaultProxy[ColorLike]] = DefaultProxy("edge_color")
    DEFAULT_EDGE_LW: ClassVar[DefaultProxy[float]] = DefaultProxy("edge_lw")

    def __post_init__(
        self,
        dendrogram: str | None,
        with_swapped_axes: bool,
        standard_scale: Literal["var", "group"] | None,
    ):
        super().__post_init__(
            dendrogram=dendrogram, with_swapped_axes=with_swapped_axes
        )
        if self.values_df is not None:
            return

        # compute mean value
        self.values_df = cast(
            pd.DataFrame,
            self.obs_tidy.groupby(level=0, observed=True)
            .mean()
            .loc[
                self.categories_order
                if self.categories_order is not None
                else self.categories
            ],
        )

        if standard_scale == "group":
            self.values_df = self.values_df.sub(self.values_df.min(1), axis=0)
            self.values_df = self.values_df.div(self.values_df.max(1), axis=0).fillna(0)
        elif standard_scale == "var":
            self.values_df -= self.values_df.min(0)
            self.values_df = (self.values_df / self.values_df.max(0)).fillna(0)
        elif standard_scale is not None:
            logg.warning("Unknown type for standard_scale, ignored")

    def style(
        self,
        cmap: str | None | Empty = _empty,
        edge_color: ColorLike | Empty = _empty,
        edge_lw: float | Empty = _empty,
    ) -> Self:
        """\
        Modifies plot visual parameters.

        Parameters
        ----------
        cmap
            String denoting matplotlib color map.
        edge_color
            Edge color between the squares of matrix plot. Default is gray
        edge_lw
            Edge line width.

        Returns
        -------
        :class:`~scanpy.pl.MatrixPlot`

        Examples
        -------

        .. plot::
            :context: close-figs

            import scanpy as sc

            adata = sc.datasets.pbmc68k_reduced()
            markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']

        Change color map and turn off edges:


        .. plot::
            :context: close-figs

            (
                sc.pl.MatrixPlot(adata, markers, groupby='bulk_labels')
                .style(cmap='Blues', edge_color='none')
                .show()
            )

        """

        if cmap is not _empty:
            self.cmap = cmap
        if edge_color is not _empty:
            self.edge_color = edge_color
        if edge_lw is not _empty:
            self.edge_lw = edge_lw

        return self

    def _mainplot(self, ax: Axes):
        # work on a copy of the dataframes. This is to avoid changes
        # on the original data frames after repetitive calls to the
        # MatrixPlot object, for example once with swap_axes and other without

        _color_df = self.values_df.copy()
        if self.var_names_idx_order is not None:
            _color_df = _color_df.iloc[:, self.var_names_idx_order]

        if self.categories_order is not None:
            _color_df = _color_df.loc[self.categories_order, :]

        if self.are_axes_swapped:
            _color_df = _color_df.T
        cmap = plt.get_cmap(self.kwds.get("cmap", self.cmap))
        if "cmap" in self.kwds:
            del self.kwds["cmap"]
        normalize = check_colornorm(
            self.vboundnorm.vmin,
            self.vboundnorm.vmax,
            self.vboundnorm.vcenter,
            self.vboundnorm.norm,
        )

        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(1.5)

        kwds = fix_kwds(
            self.kwds,
            cmap=cmap,
            edgecolor=self.edge_color,
            linewidth=self.edge_lw,
            norm=normalize,
        )
        _ = ax.pcolor(_color_df, **kwds)

        y_labels = _color_df.index
        x_labels = _color_df.columns

        y_ticks = np.arange(len(y_labels)) + 0.5
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels)

        x_ticks = np.arange(len(x_labels)) + 0.5
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, rotation=90, ha="center", minor=False)

        ax.tick_params(axis="both", labelsize="small")
        ax.grid(visible=False)

        # to be consistent with the heatmap plot, is better to
        # invert the order of the y-axis, such that the first group is on
        # top
        ax.set_ylim(len(y_labels), 0)
        ax.set_xlim(0, len(x_labels))

        return normalize


@old_positionals(
    "use_raw",
    "log",
    "num_categories",
    "figsize",
    "dendrogram",
    "title",
    "cmap",
    "colorbar_title",
    "gene_symbols",
    "var_group_positions",
    "var_group_labels",
    "var_group_rotation",
    "layer",
    "standard_scale",
    # 17 positionals are enough for backwards compatibility
)
@_doc_params(
    show_save_ax=doc_show_save_ax,
    common_plot_args=doc_common_plot_args,
    groupby_plots_args=doc_common_groupby_plot_args,
    vminmax=doc_vboundnorm,
)
def matrixplot(
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str | Sequence[str],
    *,
    use_raw: bool | None = None,
    log: bool = False,
    num_categories: int = 7,
    figsize: tuple[float, float] | None = None,
    dendrogram: str | None = None,
    title: str | None = None,
    cmap: str | None = MatrixPlot.DEFAULT_COLORMAP,
    colorbar_title: str | None = MatrixPlot.DEFAULT_COLOR_LEGEND_TITLE,
    gene_symbols: str | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
    var_group_labels: Sequence[str] | None = None,
    var_group_rotation: float | None = None,
    layer: str | None = None,
    standard_scale: Literal["var", "group"] | None = None,
    values_df: pd.DataFrame | None = None,
    swap_axes: bool = False,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: _AxesSubplot | None = None,
    return_fig: bool | None = False,
    vmin: float | None = None,
    vmax: float | None = None,
    vcenter: float | None = None,
    norm: Normalize | None = None,
    **kwds,
) -> MatrixPlot | dict[str, Axes] | None:
    """\
    Creates a heatmap of the mean expression values per group of each var_names.

    This function provides a convenient interface to the :class:`~scanpy.pl.MatrixPlot`
    class. If you need more flexibility, you should use :class:`~scanpy.pl.MatrixPlot`
    directly.

    Parameters
    ----------
    {common_plot_args}
    {groupby_plots_args}
    {show_save_ax}
    {vminmax}
    kwds
        Are passed to :func:`matplotlib.pyplot.pcolor`.

    Returns
    -------
    If `return_fig` is `True`, returns a :class:`~scanpy.pl.MatrixPlot` object,
    else if `show` is false, return axes dict

    See also
    --------
    :class:`~scanpy.pl.MatrixPlot`: The MatrixPlot class can be used to to control
        several visual parameters not available in this function.
    :func:`~scanpy.pl.rank_genes_groups_matrixplot`: to plot marker genes
        identified using the :func:`~scanpy.tl.rank_genes_groups` function.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.matrixplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    .. plot::
        :context: close-figs

        markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
        sc.pl.matrixplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Get Matrix object for fine tuning:

    .. plot::
        :context: close-figs

        mp = sc.pl.matrixplot(adata, markers, 'bulk_labels', return_fig=True)
        mp.add_totals().style(edge_color='black').show()

    The axes used can be obtained using the get_axes() method

    .. plot::
        :context: close-figs

        axes_dict = mp.get_axes()
    """

    mp = MatrixPlot(
        adata,
        var_names,
        groupby=groupby,
        use_raw=use_raw,
        log=log,
        num_categories=num_categories,
        standard_scale=standard_scale,
        title=title,
        figsize=figsize,
        gene_symbols=gene_symbols,
        var_group_positions=var_group_positions,
        var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation,
        layer=layer,
        values_df=values_df,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        vcenter=vcenter,
        norm=norm,
        dendrogram=dendrogram,
        with_swapped_axes=swap_axes,
        kwds=kwds,
    )

    mp = mp.style(cmap=cmap).legend(title=colorbar_title)
    if return_fig:
        return mp
    mp.make_figure()
    savefig_or_show(MatrixPlot.DEFAULT_SAVE_PREFIX, show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return mp.get_axes()
