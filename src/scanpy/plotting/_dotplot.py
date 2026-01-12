from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from matplotlib import colormaps

from .. import logging as logg
from .._compat import old_positionals, warn
from .._settings import settings
from .._utils import _doc_params, _empty
from ._baseplot_class import BasePlot, doc_common_groupby_plot_args
from ._docs import doc_common_plot_args, doc_show_save_ax, doc_vboundnorm
from ._utils import (
    _create_white_to_color_gradient,
    _dk,
    check_colornorm,
    fix_kwds,
    make_grid_spec,
    savefig_or_show,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import Literal, Self

    import pandas as pd
    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.colors import Colormap, Normalize

    from .._utils import Empty
    from ._baseplot_class import _VarNames
    from ._utils import ColorLike, _AxesSubplot


@_doc_params(common_plot_args=doc_common_plot_args)
class DotPlot(BasePlot):
    """Allows the visualization of two values that are encoded as dot size and color.

    The size usually represents the fraction of cells (obs)
    that have a non-zero value for genes (var).

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
    across multiple clusters.

    Parameters
    ----------
    {common_plot_args}
    title
        Title for the figure
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
    kwds
        Are passed to :func:`matplotlib.pyplot.scatter`.

    See Also
    --------
    :func:`~scanpy.pl.dotplot`: Simpler way to call DotPlot but with less options.
    :func:`~scanpy.pl.rank_genes_groups_dotplot`: to plot marker
        genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ["C1QA", "PSAP", "CD79A", "CD79B", "CST3", "LYZ"]
    >>> sc.pl.DotPlot(adata, markers, groupby="bulk_labels").show()

    Using var_names as dict:

    >>> markers = {{"T-cell": "CD3D", "B-cell": "CD79A", "myeloid": "CST3"}}
    >>> sc.pl.DotPlot(adata, markers, groupby="bulk_labels").show()

    """

    DEFAULT_SAVE_PREFIX = "dotplot_"
    # default style parameters
    DEFAULT_COLORMAP = "Reds"
    DEFAULT_COLOR_ON = "dot"
    DEFAULT_DOT_MAX = None
    DEFAULT_DOT_MIN = None
    DEFAULT_SMALLEST_DOT = 0.0
    DEFAULT_LARGEST_DOT = 200.0
    DEFAULT_DOT_EDGECOLOR = "black"
    DEFAULT_DOT_EDGELW = 0.2
    DEFAULT_SIZE_EXPONENT = 1.5

    # default legend parameters
    DEFAULT_SIZE_LEGEND_TITLE = "Fraction of cells\nin group (%)"
    DEFAULT_COLOR_LEGEND_TITLE = "Mean expression\nin group"
    DEFAULT_LEGENDS_WIDTH = 1.5  # inches
    DEFAULT_PLOT_X_PADDING = 0.8  # a unit is the distance between two x-axis ticks
    DEFAULT_PLOT_Y_PADDING = 1.0  # a unit is the distance between two y-axis ticks

    @old_positionals(
        "use_raw",
        "log",
        "num_categories",
        "categories_order",
        "title",
        "figsize",
        "gene_symbols",
        "var_group_positions",
        "var_group_labels",
        "var_group_rotation",
        "layer",
        "expression_cutoff",
        "mean_only_expressed",
        "standard_scale",
        "dot_color_df",
        "dot_size_df",
        "ax",
        "vmin",
        "vmax",
        "vcenter",
        "norm",
    )
    def __init__(  # noqa: PLR0913
        self,
        adata: AnnData,
        var_names: _VarNames | Mapping[str, _VarNames],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        expression_cutoff: float = 0.0,
        mean_only_expressed: bool = False,
        standard_scale: Literal["var", "group"] | None = None,
        dot_color_df: pd.DataFrame | None = None,
        dot_size_df: pd.DataFrame | None = None,
        ax: _AxesSubplot | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        group_colors: Mapping[str, ColorLike] | None = None,
        **kwds,
    ) -> None:
        BasePlot.__init__(
            self,
            adata,
            var_names,
            groupby,
            use_raw=use_raw,
            log=log,
            num_categories=num_categories,
            categories_order=categories_order,
            title=title,
            figsize=figsize,
            gene_symbols=gene_symbols,
            var_group_positions=var_group_positions,
            var_group_labels=var_group_labels,
            var_group_rotation=var_group_rotation,
            layer=layer,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            vcenter=vcenter,
            norm=norm,
            **kwds,
        )

        # Set default style parameters
        self.cmap = self.DEFAULT_COLORMAP
        self.dot_max = self.DEFAULT_DOT_MAX
        self.dot_min = self.DEFAULT_DOT_MIN
        self.smallest_dot = self.DEFAULT_SMALLEST_DOT
        self.largest_dot = self.DEFAULT_LARGEST_DOT
        self.color_on = self.DEFAULT_COLOR_ON
        self.size_exponent = self.DEFAULT_SIZE_EXPONENT
        self.grid = False
        self.plot_x_padding = self.DEFAULT_PLOT_X_PADDING
        self.plot_y_padding = self.DEFAULT_PLOT_Y_PADDING

        self.dot_edge_color = self.DEFAULT_DOT_EDGECOLOR
        self.dot_edge_lw = self.DEFAULT_DOT_EDGELW

        # set legend defaults
        self.color_legend_title = self.DEFAULT_COLOR_LEGEND_TITLE
        self.size_title = self.DEFAULT_SIZE_LEGEND_TITLE
        self.legends_width = self.DEFAULT_LEGENDS_WIDTH
        self.show_size_legend = True
        self.show_colorbar = True

        # Store parameters needed by helper methods and prepare the dot data.
        self.standard_scale = standard_scale
        self.expression_cutoff = expression_cutoff
        self.mean_only_expressed = mean_only_expressed
        self.dot_color_df, self.dot_size_df = self._prepare_dot_data(
            dot_color_df, dot_size_df
        )
        self.group_cmaps = self._prepare_group_cmaps(group_colors)

    def _prepare_group_cmaps(
        self, group_colors: Mapping[str, ColorLike] | None
    ) -> dict[str, Colormap] | None:
        if group_colors is None:
            return None
        group_cmaps = {}
        missing_groups = []
        for group in self.dot_color_df.index:
            if group in group_colors:
                group_cmaps[group] = _create_white_to_color_gradient(
                    group_colors[group]
                )
            else:
                group_cmaps[group] = self.cmap
                missing_groups.append(group)
        if missing_groups:
            warn(
                f"The following groups will use the default colormap as no "
                f"specific colors were assigned: {missing_groups}",
                UserWarning,
            )
        return group_cmaps

    def _prepare_dot_data(
        self, dot_color_df: pd.DataFrame | None, dot_size_df: pd.DataFrame | None
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate the dataframes for dot size and color.

        Refactored to helper to satisfy complexity checks.
        """
        # for if category defined by groupby (if any) compute for each var_name
        # 1. the fraction of cells in the category having a value >expression_cutoff
        # 2. the mean value over the category

        # 1. compute fraction of cells having value > expression_cutoff
        # transform obs_tidy into boolean matrix using the expression_cutoff
        obs_bool = self.obs_tidy > self.expression_cutoff

        # compute the sum per group which in the boolean matrix this is the number
        # of values >expression_cutoff, and divide the result by the total number of
        # values in the group (given by `count()`)
        if dot_size_df is None:
            dot_size_df = (
                obs_bool.groupby(level=0, observed=True).sum()
                / obs_bool.groupby(level=0, observed=True).count()
            )

        if dot_color_df is None:
            # 2. compute mean expression value value
            if self.mean_only_expressed:
                dot_color_df = (
                    self.obs_tidy
                    .mask(~obs_bool)
                    .groupby(level=0, observed=True)
                    .mean()
                    .fillna(0)
                )
            else:
                dot_color_df = self.obs_tidy.groupby(level=0, observed=True).mean()

            if self.standard_scale == "group":
                dot_color_df = dot_color_df.sub(dot_color_df.min(1), axis=0)
                dot_color_df = dot_color_df.div(dot_color_df.max(1), axis=0).fillna(0)
            elif self.standard_scale == "var":
                dot_color_df -= dot_color_df.min(0)
                dot_color_df = (dot_color_df / dot_color_df.max(0)).fillna(0)
            elif self.standard_scale is None:
                pass
            else:
                logg.warning("Unknown type for standard_scale, ignored")
        else:
            # check that both matrices have the same shape
            if dot_color_df.shape != dot_size_df.shape:
                logg.error(
                    "the given dot_color_df data frame has a different shape than "
                    "the data frame used for the dot size. Both data frames need "
                    "to have the same index and columns"
                )

            # Because genes (columns) can be duplicated (e.g. when the
            # same gene is reported as marker gene in two clusters)
            # they need to be removed first,
            # otherwise, the duplicated genes are further duplicated when reordering
            # Eg. A df with columns ['a', 'b', 'a'] after reordering columns
            # with df[['a', 'a', 'b']], results in a df with columns:
            # ['a', 'a', 'a', 'a', 'b']

            unique_var_names, unique_idx = np.unique(
                dot_color_df.columns, return_index=True
            )
            # remove duplicate columns
            if len(unique_var_names) != len(self.var_names):
                dot_color_df = dot_color_df.iloc[:, unique_idx]

            # get the same order for rows and columns in the dot_color_df
            # using the order from the doc_size_df
            dot_color_df = dot_color_df.loc[dot_size_df.index][dot_size_df.columns]

        # Use self.categories_order if set, else self.categories
        order = (
            self.categories_order
            if self.categories_order is not None
            else self.categories
        )
        dot_color_df, dot_size_df = (
            df.loc[order] for df in (dot_color_df, dot_size_df)
        )

        return dot_color_df, dot_size_df

    @old_positionals(
        "cmap",
        "color_on",
        "dot_max",
        "dot_min",
        "smallest_dot",
        "largest_dot",
        "dot_edge_color",
        "dot_edge_lw",
        "size_exponent",
        "grid",
        "x_padding",
        "y_padding",
    )
    def style(  # noqa: PLR0913
        self,
        *,
        cmap: Colormap | str | None | Empty = _empty,
        color_on: Literal["dot", "square"] | Empty = _empty,
        dot_max: float | None | Empty = _empty,
        dot_min: float | None | Empty = _empty,
        smallest_dot: float | Empty = _empty,
        largest_dot: float | Empty = _empty,
        dot_edge_color: ColorLike | None | Empty = _empty,
        dot_edge_lw: float | None | Empty = _empty,
        size_exponent: float | Empty = _empty,
        grid: bool | Empty = _empty,
        x_padding: float | Empty = _empty,
        y_padding: float | Empty = _empty,
    ) -> Self:
        r"""Modify plot visual parameters.

        Parameters
        ----------
        cmap
            String denoting matplotlib color map.
        color_on
            By default the color map is applied to the color of the ``"dot"``.
            Optionally, the colormap can be applied to a ``"square"`` behind the dot,
            in which case the dot is transparent and only the edge is shown.
        dot_max
            If ``None``, the maximum dot size is set to the maximum fraction value found (e.g. 0.6).
            If given, the value should be a number between 0 and 1.
            All fractions larger than dot_max are clipped to this value.
        dot_min
            If ``None``, the minimum dot size is set to 0.
            If given, the value should be a number between 0 and 1.
            All fractions smaller than dot_min are clipped to this value.
        smallest_dot
            All expression fractions with `dot_min` are plotted with this size.
        largest_dot
            All expression fractions with `dot_max` are plotted with this size.
        dot_edge_color
            Dot edge color.
            When `color_on='dot'`, ``None`` means no edge.
            When `color_on='square'`, ``None`` means that
            the edge color is white for darker colors and black for lighter background square colors.
        dot_edge_lw
            Dot edge line width.
            When `color_on='dot'`, ``None`` means no edge.
            When `color_on='square'`, ``None`` means a line width of 1.5.
        size_exponent
            Dot size is computed as:
            fraction  ** size exponent and afterwards scaled to match the
            `smallest_dot` and `largest_dot` size parameters.
            Using a different size exponent changes the relative sizes of the dots
            to each other.
        grid
            Set to true to show grid lines. By default grid lines are not shown.
            Further configuration of the grid lines can be achieved directly on the
            returned ax.
        x_padding
            Space between the plot left/right borders and the dots center. A unit
            is the distance between the x ticks. Only applied when color_on = dot
        y_padding
            Space between the plot top/bottom borders and the dots center. A unit is
            the distance between the y ticks. Only applied when color_on = dot

        Returns
        -------
        :class:`~scanpy.pl.DotPlot`

        Examples
        --------
        >>> import scanpy as sc
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']

        Change color map and apply it to the square behind the dot

        >>> sc.pl.DotPlot(adata, markers, groupby='bulk_labels') \
        ...     .style(cmap='RdBu_r', color_on='square').show()

        Add edge to dots and plot a grid

        >>> sc.pl.DotPlot(adata, markers, groupby='bulk_labels') \
        ...     .style(dot_edge_color='black', dot_edge_lw=1, grid=True) \
        ...     .show()

        """
        super().style(cmap=cmap)

        if dot_max is not _empty:
            self.dot_max = dot_max
        if dot_min is not _empty:
            self.dot_min = dot_min
        if smallest_dot is not _empty:
            self.smallest_dot = smallest_dot
        if largest_dot is not _empty:
            self.largest_dot = largest_dot
        if color_on is not _empty:
            self.color_on = color_on
        if size_exponent is not _empty:
            self.size_exponent = size_exponent
        if dot_edge_color is not _empty:
            self.dot_edge_color = dot_edge_color
        if dot_edge_lw is not _empty:
            self.dot_edge_lw = dot_edge_lw
        if grid is not _empty:
            self.grid = grid
        if x_padding is not _empty:
            self.plot_x_padding = x_padding
        if y_padding is not _empty:
            self.plot_y_padding = y_padding

        return self

    @old_positionals(
        "show",
        "show_size_legend",
        "show_colorbar",
        "size_title",
        "colorbar_title",
        "width",
    )
    def legend(
        self,
        *,
        show: bool | None = True,
        show_size_legend: bool | None = True,
        show_colorbar: bool | None = True,
        size_title: str | None = DEFAULT_SIZE_LEGEND_TITLE,
        colorbar_title: str | None = DEFAULT_COLOR_LEGEND_TITLE,
        width: float | None = DEFAULT_LEGENDS_WIDTH,
    ) -> Self:
        r"""Configure dot size and the colorbar legends.

        Parameters
        ----------
        show
            Set to `False` to hide the default plot of the legends. This sets the
            legend width to zero, which will result in a wider main plot.
        show_size_legend
            Set to `False` to hide the dot size legend
        show_colorbar
            Set to `False` to hide the colorbar legend
        size_title
            Title for the dot size legend. Use ``\n`` to add line breaks. Appears on top
            of dot sizes
        colorbar_title
            Title for the color bar. Use ``\n`` to add line breaks. Appears on top of the
            color bar
        width
            Width of the legends area. The unit is the same as in matplotlib (inches).

        Returns
        -------
        :class:`~scanpy.pl.DotPlot`

        Examples
        --------
        Set color bar title:

        >>> import scanpy as sc
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = {"T-cell": "CD3D", "B-cell": "CD79A", "myeloid": "CST3"}
        >>> dp = sc.pl.DotPlot(adata, markers, groupby="bulk_labels")
        >>> dp.legend(colorbar_title="log(UMI counts + 1)").show()

        """
        if not show:
            # turn of legends by setting width to 0
            self.legends_width = 0
        else:
            self.color_legend_title = colorbar_title
            self.size_title = size_title
            self.legends_width = width
            self.show_size_legend = show_size_legend
            self.show_colorbar = show_colorbar

        return self

    def _plot_size_legend(self, size_legend_ax: Axes):
        # for the dot size legend, use step between dot_max and dot_min
        # based on how different they are.
        diff = self.dot_max - self.dot_min
        if 0.3 < diff <= 0.6:
            step = 0.1
        elif diff <= 0.3:
            step = 0.05
        else:
            step = 0.2
        # a descending range that is afterwards inverted is used
        # to guarantee that dot_max is in the legend.
        size_range = np.arange(self.dot_max, self.dot_min, step * -1)[::-1]
        if self.dot_min != 0 or self.dot_max != 1:
            dot_range = self.dot_max - self.dot_min
            size_values = (size_range - self.dot_min) / dot_range
        else:
            size_values = size_range

        size = size_values**self.size_exponent
        size = size * (self.largest_dot - self.smallest_dot) + self.smallest_dot

        # plot size bar
        size_legend_ax.scatter(
            np.arange(len(size)) + 0.5,
            np.repeat(0, len(size)),
            s=size,
            color="gray",
            edgecolor="black",
            linewidth=self.dot_edge_lw,
            zorder=100,
        )
        size_legend_ax.set_xticks(np.arange(len(size)) + 0.5)
        labels = [f"{np.round((x * 100), decimals=0).astype(int)}" for x in size_range]
        size_legend_ax.set_xticklabels(labels, fontsize="small")

        # remove y ticks and labels
        size_legend_ax.tick_params(
            axis="y", left=False, labelleft=False, labelright=False
        )

        # remove surrounding lines
        size_legend_ax.spines["right"].set_visible(False)
        size_legend_ax.spines["top"].set_visible(False)
        size_legend_ax.spines["left"].set_visible(False)
        size_legend_ax.spines["bottom"].set_visible(False)
        size_legend_ax.grid(visible=False)

        ymax = size_legend_ax.get_ylim()[1]
        size_legend_ax.set_ylim(-1.05 - self.largest_dot * 0.003, 4)
        size_legend_ax.set_title(self.size_title, y=ymax + 0.45, size="small")

        xmin, xmax = size_legend_ax.get_xlim()
        size_legend_ax.set_xlim(xmin - 0.15, xmax + 0.5)

    def _plot_legend(self, legend_ax, return_ax_dict, normalize):
        # to maintain the fixed height size of the legends, a
        # spacer of variable height is added at the bottom.
        # The structure for the legends is:
        # first row: variable space to keep the other rows of
        #            the same size (avoid stretching)
        # second row: legend for dot size
        # third row: spacer to avoid color and size legend titles to overlap
        # fourth row: colorbar

        # Define base heights for legend components as a fraction of figure height
        cbar_legend_height = self.min_figure_height * 0.08
        size_legend_height = self.min_figure_height * 0.27
        spacer_height = self.min_figure_height * 0.3

        # If group_colors is used, dynamically calculate the total height needed for all colorbars
        if self.group_cmaps is not None:
            per_cbar_height = (
                self.min_figure_height * 0.12
            )  # Use a slightly larger height for better spacing
            n_cbars = len(self.dot_color_df.index)
            cbar_legend_height = per_cbar_height * n_cbars

        # Calculate the height of the top spacer to push content down
        top_spacer_height = (
            self.height - size_legend_height - cbar_legend_height - spacer_height
        )
        top_spacer_height = max(top_spacer_height, 0)  # prevent negative height

        # Create the 4-row GridSpec for the legend area
        height_ratios = [
            top_spacer_height,
            size_legend_height,
            spacer_height,
            cbar_legend_height,
        ]
        fig, legend_gs = make_grid_spec(
            legend_ax, nrows=4, ncols=1, height_ratios=height_ratios
        )
        # Hide the frame of the main legend container axis for a cleaner look
        legend_ax.set_axis_off()

        # Plot size legend into the second row of the grid
        if self.show_size_legend:
            size_legend_ax = fig.add_subplot(legend_gs[1])
            self._plot_size_legend(size_legend_ax)
            return_ax_dict["size_legend_ax"] = size_legend_ax

        # Plot colorbar(s) into the fourth row of the grid
        if self.show_colorbar:
            if self.group_cmaps is None:
                color_legend_ax = fig.add_subplot(legend_gs[3])
                self._plot_colorbar(color_legend_ax, normalize)
                return_ax_dict["color_legend_ax"] = color_legend_ax
            else:
                self._plot_stacked_colorbars(fig, legend_gs[3], normalize)
                return_ax_dict["color_legend_ax"] = legend_ax

    def _plot_stacked_colorbars(self, fig, colorbar_area_spec, normalize):
        """Plot the stacked colorbars legend when using group_colors."""
        import matplotlib as mpl
        import matplotlib.colorbar
        from matplotlib.cm import ScalarMappable

        plotted_groups = self.dot_color_df.index
        groups_to_plot = list(plotted_groups)
        n_cbars = len(groups_to_plot)

        # Create a sub-grid just for the colorbars
        # Create an empty column to keep colorbars at 3/4 of legend width (1.5 like default with dp.legend_width = 2.0)
        colorbar_gs = colorbar_area_spec.subgridspec(
            n_cbars, 2, hspace=0.6, width_ratios=[3, 1]
        )

        # Create a dedicated normalizer for the legend
        vmin = self.dot_color_df.values.min()
        vmax = self.dot_color_df.values.max()
        legend_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        for i, group_name in enumerate(groups_to_plot):
            ax = fig.add_subplot(
                colorbar_gs[i, 0]
            )  # Place the colorbar Axes in the first, wider column

            # self.group_cmaps[group_name] is already a Colormap object (or string from fallback)
            cmap = self.group_cmaps[group_name]
            if isinstance(cmap, str):
                cmap = colormaps.get_cmap(cmap)

            mappable = ScalarMappable(norm=legend_norm, cmap=cmap)

            cb = matplotlib.colorbar.Colorbar(
                ax, mappable=mappable, orientation="horizontal"
            )
            cb.ax.xaxis.set_tick_params(labelsize="small")

            ax.text(
                1.1,
                0.5,
                group_name,
                ha="left",
                va="center",
                transform=ax.transAxes,
                fontsize="small",
            )

            if i == 0:
                cb.ax.set_title(self.color_legend_title, fontsize="small")

            if i < n_cbars - 1:
                cb.ax.xaxis.set_ticklabels([])
                cb.ax.xaxis.set_ticks([])

    def _mainplot(self, ax: Axes):
        # work on a copy of the dataframes. This is to avoid changes
        # on the original data frames after repetitive calls to the
        # DotPlot object, for example once with swap_axes and other without

        _color_df = self.dot_color_df.copy()
        _size_df = self.dot_size_df.copy()
        if self.var_names_idx_order is not None:
            _color_df = _color_df.iloc[:, self.var_names_idx_order]
            _size_df = _size_df.iloc[:, self.var_names_idx_order]

        if self.categories_order is not None:
            _color_df = _color_df.loc[self.categories_order, :]
            _size_df = _size_df.loc[self.categories_order, :]

        if self.are_axes_swapped:
            _size_df = _size_df.T
            _color_df = _color_df.T
        self.cmap = self.kwds.pop("cmap", self.cmap)

        normalize, dot_min, dot_max = self._dotplot(
            _size_df,
            _color_df,
            ax,
            cmap=self.cmap,
            group_cmaps=self.group_cmaps,
            color_on=self.color_on,
            dot_max=self.dot_max,
            dot_min=self.dot_min,
            standard_scale=self.standard_scale,
            edge_color=self.dot_edge_color,
            edge_lw=self.dot_edge_lw,
            smallest_dot=self.smallest_dot,
            largest_dot=self.largest_dot,
            size_exponent=self.size_exponent,
            grid=self.grid,
            x_padding=self.plot_x_padding,
            y_padding=self.plot_y_padding,
            vmin=self.vboundnorm.vmin,
            vmax=self.vboundnorm.vmax,
            vcenter=self.vboundnorm.vcenter,
            norm=self.vboundnorm.norm,
            are_axes_swapped=self.are_axes_swapped,
            **self.kwds,
        )

        self.dot_min, self.dot_max = dot_min, dot_max
        return normalize

    @staticmethod
    def _dotplot(  # noqa: PLR0912, PLR0913, PLR0915
        dot_size: pd.DataFrame,
        dot_color: pd.DataFrame,
        dot_ax: Axes,
        *,
        cmap: Colormap | str | None,
        group_cmaps: Mapping[str, str] | None,
        color_on: Literal["dot", "square"],
        dot_max: float | None,
        dot_min: float | None,
        standard_scale: Literal["var", "group"] | None,
        smallest_dot: float,
        largest_dot: float,
        size_exponent: float,
        edge_color: ColorLike | None,
        edge_lw: float | None,
        grid: bool,
        x_padding: float,
        y_padding: float,
        vmin: float | None,
        vmax: float | None,
        vcenter: float | None,
        norm: Normalize | None,
        are_axes_swapped: bool | None,
        **kwds,
    ):
        """Make a *dot plot* given two data frames.

        One containing the dot size and other containing the dot color.
        The indices and columns of the data frame are used to label the output image.

        The dots are plotted using :func:`matplotlib.pyplot.scatter`. Thus, additional
        arguments can be passed.

        Parameters
        ----------
        dot_size
            Data frame containing the dot_size.
        dot_color
            Data frame containing the dot_color, should have the same,
            shape, columns and indices as dot_size.
        dot_ax
            matplotlib axis
        cmap
        color_on
        dot_max
        dot_min
        standard_scale
        smallest_dot
        edge_color
        edge_lw
        grid
        x_padding
        y_padding
            See `style`
        kwds
            Are passed to :func:`matplotlib.pyplot.scatter`.

        Returns
        -------
        matplotlib.colors.Normalize, dot_min, dot_max

        """
        assert dot_size.shape == dot_color.shape, (
            "please check that dot_size and dot_color dataframes have the same shape"
        )

        assert list(dot_size.index) == list(dot_color.index), (
            "please check that dot_size and dot_color dataframes have the same index"
        )

        assert list(dot_size.columns) == list(dot_color.columns), (
            "please check that the dot_size "
            "and dot_color dataframes have the same columns"
        )
        match are_axes_swapped, standard_scale:
            case True, "group":
                group_axis = 1
            case True, "var":
                group_axis = 0
            case _, "group":
                group_axis = 0
            case _, "var":
                group_axis = 1
        if standard_scale is not None:
            dot_color = dot_color.sub(
                dot_color.min((group_axis + 1) % 2), axis=group_axis
            )
            dot_color = dot_color.div(
                dot_color.max((group_axis + 1) % 2), axis=group_axis
            ).fillna(0)
        # make scatter plot in which
        # x = var_names
        # y = groupby category
        # size = fraction
        # color = mean expression

        # +0.5 in y and x to set the dot center at 0.5 multiples
        # this facilitates dendrogram and totals alignment for
        # matrixplot, dotplot and stackec_violin using the same coordinates.
        y, x = np.indices(dot_color.shape)
        y = y.flatten() + 0.5
        x = x.flatten() + 0.5
        frac = dot_size.values.flatten()
        mean_flat = dot_color.values.flatten()
        cmap = colormaps.get_cmap(cmap)
        if dot_max is None:
            dot_max = np.ceil(max(frac) * 10) / 10
        elif dot_max < 0 or dot_max > 1:
            msg = "`dot_max` value has to be between 0 and 1"
            raise ValueError(msg)
        if dot_min is None:
            dot_min = 0
        elif dot_min < 0 or dot_min > 1:
            msg = "`dot_min` value has to be between 0 and 1"
            raise ValueError(msg)

        if dot_min != 0 or dot_max != 1:
            # clip frac between dot_min and  dot_max
            frac = np.clip(frac, dot_min, dot_max)
            old_range = dot_max - dot_min
            # re-scale frac between 0 and 1
            frac = (frac - dot_min) / old_range

        size = frac**size_exponent
        # rescale size to match smallest_dot and largest_dot
        size = size * (largest_dot - smallest_dot) + smallest_dot
        normalize = check_colornorm(vmin, vmax, vcenter, norm)

        if group_cmaps is None:
            # Plotting logic for single colormap
            if color_on == "square":
                if edge_color is None:
                    from seaborn.utils import relative_luminance

                    # use either black or white for the edge color
                    # depending on the luminance of the background
                    # square color
                    edge_color = []
                    for color_value in cmap(normalize(mean_flat)):
                        lum = relative_luminance(color_value)
                        edge_color.append(".15" if lum > 0.408 else "w")

                edge_lw = 1.5 if edge_lw is None else edge_lw

                # first make a heatmap similar to `sc.pl.matrixplot`
                # (squares with the asigned colormap). Circles will be plotted
                # on top
                dot_ax.pcolor(dot_color.values, cmap=cmap, norm=normalize)
                for axis in ["top", "bottom", "left", "right"]:
                    dot_ax.spines[axis].set_linewidth(1.5)
                # Create a temporary kwargs dict for this group's scatter call
                # to avoid modifying the original kwds dictionary within the loop.
                kwds_scatter = fix_kwds(
                    kwds,
                    s=size,
                    linewidth=edge_lw,
                    facecolor="none",
                    edgecolor=edge_color,
                )
                dot_ax.scatter(x, y, **kwds_scatter)
            else:
                edge_color = "none" if edge_color is None else edge_color
                edge_lw = 0.0 if edge_lw is None else edge_lw
                color = cmap(normalize(mean_flat))
                kwds_scatter = fix_kwds(
                    kwds,
                    s=size,
                    color=color,
                    linewidth=edge_lw,
                    edgecolor=edge_color,
                )
                dot_ax.scatter(x, y, **kwds_scatter)
        else:
            # Plotting logic for group-specific colormaps
            groups_iter = dot_color.columns if are_axes_swapped else dot_color.index
            n_vars = dot_color.shape[0] if are_axes_swapped else dot_color.shape[1]
            n_groups = len(groups_iter)

            # Here we loop through each group and plot it with its own cmap
            for group_idx, group_name in enumerate(groups_iter):
                group_cmap = group_cmaps[group_name]
                # Handle fallback case where group_cmap might be a string
                if isinstance(group_cmap, str):
                    group_cmap = colormaps.get_cmap(group_cmap)

                # Slice the flattened data arrays correctly depending on orientation
                if not are_axes_swapped:
                    # Slicing data for a whole row
                    indices = slice(group_idx * n_vars, (group_idx + 1) * n_vars)
                else:
                    # Slicing data for a whole column
                    indices = slice(group_idx, None, n_groups)

                x_group = x[indices]
                y_group = y[indices]
                size_group = size[indices]
                mean_group = mean_flat[indices]

                color = group_cmap(normalize(mean_group))
                kwds_scatter = fix_kwds(
                    kwds,
                    s=size_group,
                    color=color,
                    linewidth=edge_lw,
                    edgecolor=edge_color,
                )
                dot_ax.scatter(x_group, y_group, **kwds_scatter)

        y_ticks = np.arange(dot_color.shape[0]) + 0.5
        dot_ax.set_yticks(y_ticks)
        dot_ax.set_yticklabels(
            [dot_color.index[idx] for idx, _ in enumerate(y_ticks)], minor=False
        )

        x_ticks = np.arange(dot_color.shape[1]) + 0.5
        dot_ax.set_xticks(x_ticks)
        dot_ax.set_xticklabels(
            [dot_color.columns[idx] for idx, _ in enumerate(x_ticks)],
            rotation=90,
            ha="center",
            minor=False,
        )
        dot_ax.tick_params(axis="both", labelsize="small")
        dot_ax.grid(visible=False)

        # to be consistent with the heatmap plot, is better to
        # invert the order of the y-axis, such that the first group is on
        # top
        dot_ax.set_ylim(dot_color.shape[0], 0)
        dot_ax.set_xlim(0, dot_color.shape[1])

        if color_on == "dot":
            # add padding to the x and y lims when the color is not in the square
            # default y range goes from 0.5 to num cols + 0.5
            # and default x range goes from 0.5 to num rows + 0.5, thus
            # the padding needs to be corrected.
            x_padding = x_padding - 0.5
            y_padding = y_padding - 0.5
            dot_ax.set_ylim(dot_color.shape[0] + y_padding, -y_padding)

            dot_ax.set_xlim(-x_padding, dot_color.shape[1] + x_padding)

        if grid:
            dot_ax.grid(visible=True, color="gray", linewidth=0.1)
            dot_ax.set_axisbelow(True)

        return normalize, dot_min, dot_max


@old_positionals(
    "use_raw",
    "log",
    "num_categories",
    "expression_cutoff",
    "mean_only_expressed",
    "cmap",
    "dot_max",
    "dot_min",
    "standard_scale",
    "smallest_dot",
    "title",
    "colorbar_title",
    "size_title",
    # No need to have backwards compat for > 16 positional parameters
)
@_doc_params(
    show_save_ax=doc_show_save_ax,
    common_plot_args=doc_common_plot_args,
    groupby_plots_args=doc_common_groupby_plot_args,
    vminmax=doc_vboundnorm,
)
def dotplot(  # noqa: PLR0913
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str | Sequence[str],
    *,
    use_raw: bool | None = None,
    log: bool = False,
    num_categories: int = 7,
    categories_order: Sequence[str] | None = None,
    expression_cutoff: float = 0.0,
    mean_only_expressed: bool = False,
    standard_scale: Literal["var", "group"] | None = None,
    title: str | None = None,
    colorbar_title: str | None = DotPlot.DEFAULT_COLOR_LEGEND_TITLE,
    size_title: str | None = DotPlot.DEFAULT_SIZE_LEGEND_TITLE,
    figsize: tuple[float, float] | None = None,
    dendrogram: bool | str = False,
    gene_symbols: str | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
    var_group_labels: Sequence[str] | None = None,
    var_group_rotation: float | None = None,
    layer: str | None = None,
    swap_axes: bool | None = False,
    dot_color_df: pd.DataFrame | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    ax: _AxesSubplot | None = None,
    return_fig: bool | None = False,
    vmin: float | None = None,
    vmax: float | None = None,
    vcenter: float | None = None,
    norm: Normalize | None = None,
    # Style parameters
    cmap: Colormap | str | None = DotPlot.DEFAULT_COLORMAP,
    group_colors: Mapping[str, ColorLike] | None = None,
    dot_max: float | None = DotPlot.DEFAULT_DOT_MAX,
    dot_min: float | None = DotPlot.DEFAULT_DOT_MIN,
    smallest_dot: float = DotPlot.DEFAULT_SMALLEST_DOT,
    **kwds,
) -> DotPlot | dict | None:
    r"""Make a *dot plot* of the expression values of `var_names`.

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
    across  multiple clusters.

    This function provides a convenient interface to the :class:`~scanpy.pl.DotPlot`
    class. If you need more flexibility, you should use :class:`~scanpy.pl.DotPlot`
    directly.

    Parameters
    ----------
    {common_plot_args}
    {groupby_plots_args}
    size_title
        Title for the size legend. New line character (\n) can be used.
    expression_cutoff
        Expression cutoff that is used for binarizing the gene expression and
        determining the fraction of cells expressing given genes. A gene is
        expressed only if the expression value is greater than this threshold.
    mean_only_expressed
        If True, gene expression is averaged only over the cells
        expressing the given genes.
    group_colors
        A mapping of group names to colors.
        e.g. `{{'T-cell': 'blue', 'B-cell': '#aa40fc'}}`.
        Colors can be specified as any valid matplotlib color.
        If `group_colors` is used, a colormap is generated from white
        to the given color for each group.
        If a group is not present in the dictionary, the value of `cmap`
        is used.
    dot_max
        If ``None``, the maximum dot size is set to the maximum fraction value found
        (e.g. 0.6). If given, the value should be a number between 0 and 1.
        All fractions larger than dot_max are clipped to this value.
    dot_min
        If ``None``, the minimum dot size is set to 0. If given,
        the value should be a number between 0 and 1.
        All fractions smaller than dot_min are clipped to this value.
    smallest_dot
        All expression levels with `dot_min` are plotted with this size.
    {show_save_ax}
    {vminmax}
    kwds
        Are passed to :func:`matplotlib.pyplot.scatter`.

    Returns
    -------
    If `return_fig` is `True`, returns a :class:`~scanpy.pl.DotPlot` object,
    else if `show` is false, return axes dict

    See Also
    --------
    :class:`~scanpy.pl.DotPlot`: The DotPlot class can be used to to control
        several visual parameters not available in this function.
    :func:`~scanpy.pl.rank_genes_groups_dotplot`: to plot marker genes
        identified using the :func:`~scanpy.tl.rank_genes_groups` function.

    Examples
    --------
    Create a dot plot using the given markers and the PBMC example dataset grouped by
    the category 'bulk_labels'.

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    .. plot::
        :context: close-figs

        markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
        sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)

    Get DotPlot object for fine tuning

    .. plot::
        :context: close-figs

        dp = sc.pl.dotplot(adata, markers, 'bulk_labels', return_fig=True)
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).show()

    The axes used can be obtained using the get_axes() method

    .. code-block:: python

        axes_dict = dp.get_axes()
        print(axes_dict)

    """
    # backwards compatibility: previous version of dotplot used `color_map`
    # instead of `cmap`
    cmap = kwds.pop("color_map", cmap)

    # Warn if both cmap and group_colors are specified
    if group_colors is not None and cmap != DotPlot.DEFAULT_COLORMAP:
        warn(
            "Both `cmap` and `group_colors` are specified. "
            "`group_colors` takes precedence for the specified groups.",
            UserWarning,
        )

    dp = DotPlot(
        adata,
        var_names,
        groupby,
        use_raw=use_raw,
        log=log,
        num_categories=num_categories,
        categories_order=categories_order,
        expression_cutoff=expression_cutoff,
        mean_only_expressed=mean_only_expressed,
        standard_scale=standard_scale,
        title=title,
        figsize=figsize,
        gene_symbols=gene_symbols,
        var_group_positions=var_group_positions,
        var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation,
        layer=layer,
        dot_color_df=dot_color_df,
        group_colors=group_colors,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        vcenter=vcenter,
        norm=norm,
        **kwds,
    )

    if dendrogram:
        dp.add_dendrogram(dendrogram_key=_dk(dendrogram))
    if swap_axes:
        dp.swap_axes()

    dp = dp.style(
        cmap=cmap,
        dot_max=dot_max,
        dot_min=dot_min,
        smallest_dot=smallest_dot,
        dot_edge_lw=kwds.pop("linewidth", _empty),
    ).legend(
        colorbar_title=colorbar_title, size_title=size_title, width=2.0
    )  # Width 2.0 to avoid size legend circles to overlap

    if return_fig:
        return dp
    else:
        dp.make_figure()
        savefig_or_show(DotPlot.DEFAULT_SAVE_PREFIX, show=show, save=save)
        show = settings.autoshow if show is None else show
        if not show:
            return dp.get_axes()
