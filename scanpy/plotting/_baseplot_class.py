"""BasePlot for dotplot, matrixplot and stacked_violin
"""
import collections.abc as cabc
from collections import namedtuple
from typing import Optional, Union, Mapping  # Special
from typing import Sequence, Iterable  # ABCs
from typing import Tuple  # Classes

import numpy as np
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib import pyplot as pl
from matplotlib import gridspec
from matplotlib.colors import Normalize
from warnings import warn

from .. import logging as logg
from .._compat import Literal
from ._utils import make_grid_spec, check_colornorm
from ._utils import ColorLike, _AxesSubplot
from ._anndata import _plot_dendrogram, _get_dendrogram_key, _prepare_dataframe

_VarNames = Union[str, Sequence[str]]

doc_common_groupby_plot_args = """\
title
    Title for the figure
colorbar_title
    Title for the color bar. New line character (\\n) can be used.
cmap
    String denoting matplotlib color map.
standard_scale
    Whether or not to standardize the given dimension between 0 and 1, meaning for
    each variable or group, subtract the minimum and divide each by its maximum.
swap_axes
     By default, the x axis contains `var_names` (e.g. genes) and the y axis
     the `groupby` categories. By setting `swap_axes` then x are the
     `groupby` categories and y the `var_names`.
return_fig
    Returns :class:`DotPlot` object. Useful for fine-tuning
    the plot. Takes precedence over `show=False`.
"""


class BasePlot(object):
    """\
    Generic class for the visualization of AnnData categories and
    selected `var` (features or genes).

    Takes care of the visual location of a main plot, additional plots
    in the margins (e.g. dendrogram, margin totals) and legends. Also
    understand how to adapt the visual parameter if the plot is rotated

    Classed based on BasePlot implement their own _mainplot() method.

    The BasePlot works by method chaining. For example:
    BasePlot(adata, ...).legend(title='legend').style(cmap='binary').show()
    """

    DEFAULT_SAVE_PREFIX = 'baseplot_'
    MIN_FIGURE_HEIGHT = 2.5
    DEFAULT_CATEGORY_HEIGHT = 0.35
    DEFAULT_CATEGORY_WIDTH = 0.37

    # gridspec parameter. Sets the space between mainplot, dendrogram and legend
    DEFAULT_WSPACE = 0

    DEFAULT_COLORMAP = 'winter'
    DEFAULT_LEGENDS_WIDTH = 1.5
    DEFAULT_COLOR_LEGEND_TITLE = 'Expression\nlevel in group'

    MAX_NUM_CATEGORIES = 500  # maximum number of categories allowed to be plotted

    def __init__(
        self,
        adata: AnnData,
        var_names: Union[_VarNames, Mapping[str, _VarNames]],
        groupby: Union[str, Sequence[str]],
        use_raw: Optional[bool] = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Optional[Sequence[str]] = None,
        title: Optional['str'] = None,
        figsize: Optional[Tuple[float, float]] = None,
        gene_symbols: Optional[str] = None,
        var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
        var_group_labels: Optional[Sequence[str]] = None,
        var_group_rotation: Optional[float] = None,
        layer: Optional[str] = None,
        ax: Optional[_AxesSubplot] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        vcenter: Optional[float] = None,
        norm: Optional[Normalize] = None,
        **kwds,
    ):
        self.var_names = var_names
        self.var_group_labels = var_group_labels
        self.var_group_positions = var_group_positions
        self.var_group_rotation = var_group_rotation
        self.width, self.height = figsize if figsize is not None else (None, None)

        self.has_var_groups = (
            True
            if var_group_positions is not None and len(var_group_positions) > 0
            else False
        )

        self._update_var_groups()

        self.categories, self.obs_tidy = _prepare_dataframe(
            adata,
            self.var_names,
            groupby,
            use_raw,
            log,
            num_categories,
            layer=layer,
            gene_symbols=gene_symbols,
        )
        if len(self.categories) > self.MAX_NUM_CATEGORIES:
            warn(
                f"Over {self.MAX_NUM_CATEGORIES} categories found. "
                "Plot would be very large."
            )

        if categories_order is not None:
            if set(self.obs_tidy.index.categories) != set(categories_order):
                logg.error(
                    "Please check that the categories given by "
                    "the `order` parameter match the categories that "
                    "want to be reordered.\n\n"
                    "Mismatch: "
                    f"{set(self.obs_tidy.index.categories).difference(categories_order)}\n\n"
                    f"Given order categories: {categories_order}\n\n"
                    f"{groupby} categories: {list(self.obs_tidy.index.categories)}\n"
                )
                return

        self.adata = adata
        self.groupby = [groupby] if isinstance(groupby, str) else groupby
        self.log = log
        self.kwds = kwds

        VBoundNorm = namedtuple('VBoundNorm', ['vmin', 'vmax', 'vcenter', 'norm'])
        self.vboundnorm = VBoundNorm(vmin=vmin, vmax=vmax, vcenter=vcenter, norm=norm)

        # set default values for legend
        self.color_legend_title = self.DEFAULT_COLOR_LEGEND_TITLE
        self.legends_width = self.DEFAULT_LEGENDS_WIDTH

        # set style defaults
        self.cmap = self.DEFAULT_COLORMAP

        # style default parameters
        self.are_axes_swapped = False
        self.categories_order = categories_order
        self.var_names_idx_order = None

        self.wspace = self.DEFAULT_WSPACE

        # minimum height required for legends to plot properly
        self.min_figure_height = self.MIN_FIGURE_HEIGHT

        self.fig_title = title

        self.group_extra_size = 0
        self.plot_group_extra = None
        # after .render() is called the fig value is assigned and ax_dict
        # contains a dictionary of the axes used in the plot
        self.fig = None
        self.ax_dict = None
        self.ax = ax

    def swap_axes(self, swap_axes: Optional[bool] = True):
        """
        Plots a transposed image.

        By default, the x axis contains `var_names` (e.g. genes) and the y
        axis the `groupby` categories. By setting `swap_axes` then x are
        the `groupby` categories and y the `var_names`.

        Parameters
        ----------
        swap_axes
            Boolean to turn on (True) or off (False) 'swap_axes'. Default True


        Returns
        -------
        BasePlot

        """
        self.DEFAULT_CATEGORY_HEIGHT, self.DEFAULT_CATEGORY_WIDTH = (
            self.DEFAULT_CATEGORY_WIDTH,
            self.DEFAULT_CATEGORY_HEIGHT,
        )

        self.are_axes_swapped = swap_axes
        return self

    def add_dendrogram(
        self,
        show: Optional[bool] = True,
        dendrogram_key: Optional[str] = None,
        size: Optional[float] = 0.8,
    ):
        """\
        Show dendrogram based on the hierarchical clustering between the `groupby`
        categories. Categories are reordered to match the dendrogram order.

        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `sc.tl.dendrogram` has not been called previously the function is called
        with default parameters.

        The dendrogram is by default shown on the right side of the plot or on top
        if the axes are swapped.

        `var_names` are reordered to produce a more pleasing output if:
            * The data contains `var_groups`
            * the `var_groups` match the categories.

        The previous conditions happen by default when using Plot
        to show the results from :func:`~scanpy.tl.rank_genes_groups` (aka gene markers), by
        calling `scanpy.tl.rank_genes_groups_(plot_name)`.


        Parameters
        ----------
        show
            Boolean to turn on (True) or off (False) 'add_dendrogram'
        dendrogram_key
            Needed if `sc.tl.dendrogram` saved the dendrogram using a key different
            than the default name.
        size
            size of the dendrogram. Corresponds to width when dendrogram shown on
            the right of the plot, or height when shown on top. The unit is the same
            as in matplotlib (inches).

        Returns
        -------
        BasePlot


        Examples
        --------
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        >>> sc.pl.BasePlot(adata, markers, groupby='bulk_labels').add_dendrogram().show()

        """

        if not show:
            self.plot_group_extra = None
            return self

        if self.groupby is None or len(self.categories) <= 2:
            # dendrogram can only be computed  between groupby categories
            logg.warning(
                "Dendrogram not added. Dendrogram is added only "
                "when the number of categories to plot > 2"
            )
            return self

        self.group_extra_size = size

        # to correctly plot the dendrogram the categories need to be ordered
        # according to the dendrogram ordering.
        self._reorder_categories_after_dendrogram(dendrogram_key)

        dendro_ticks = np.arange(len(self.categories)) + 0.5

        self.group_extra_size = size
        self.plot_group_extra = {
            'kind': 'dendrogram',
            'width': size,
            'dendrogram_key': dendrogram_key,
            'dendrogram_ticks': dendro_ticks,
        }
        return self

    def add_totals(
        self,
        show: Optional[bool] = True,
        sort: Literal['ascending', 'descending'] = None,
        size: Optional[float] = 0.8,
        color: Optional[Union[ColorLike, Sequence[ColorLike]]] = None,
    ):
        """\
        Show barplot for the number of cells in in `groupby` category.

        The barplot is by default shown on the right side of the plot or on top
        if the axes are swapped.


        Parameters
        ----------
        show
            Boolean to turn on (True) or off (False) 'add_totals'
        sort
            Set to either 'ascending' or 'descending' to reorder the categories
            by cell number
        size
            size of the barplot. Corresponds to width when shown on
            the right of the plot, or height when shown on top. The unit is the same
            as in matplotlib (inches).
        color
            Color for the bar plots or list of colors for each of the bar plots.
            By default, each bar plot uses the colors assigned in
            `adata.uns[{groupby}_colors]`.


        Returns
        -------
        BasePlot


        Examples
        --------
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        >>> sc.pl.BasePlot(adata, markers, groupby='bulk_labels').add_totals().show()
        """
        self.group_extra_size = size

        if not show:
            # hide totals
            self.plot_group_extra = None
            self.group_extra_size = 0
            return self

        _sort = True if sort is not None else False
        _ascending = True if sort == 'ascending' else False
        counts_df = self.obs_tidy.index.value_counts(sort=_sort, ascending=_ascending)

        if _sort:
            self.categories_order = counts_df.index

        self.plot_group_extra = {
            'kind': 'group_totals',
            'width': size,
            'sort': sort,
            'counts_df': counts_df,
            'color': color,
        }
        return self

    def style(self, cmap: Optional[str] = DEFAULT_COLORMAP):
        """\
        Set visual style parameters

        Parameters
        ----------
        cmap
            colormap

        Returns
        -------
        BasePlot
        """

        self.cmap = cmap

    def legend(
        self,
        show: Optional[bool] = True,
        title: Optional[str] = DEFAULT_COLOR_LEGEND_TITLE,
        width: Optional[float] = DEFAULT_LEGENDS_WIDTH,
    ):
        """\
        Configure legend parameters

        Parameters
        ----------
        show
            Set to 'False' to hide the default plot of the legend. This sets the
            legend width to zero which will result in a wider main plot.
        title
            Legend title. Appears on top of the color bar. Use '\\n' to add line breaks.
        width
            Width of the legend. The unit is the same as in matplotlib (inches)

        Returns
        -------
        BasePlot


        Examples
        --------

        Set legend title:

        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        >>> dp = sc.pl.BasePlot(adata, markers, groupby='bulk_labels')
        >>> dp.legend(colorbar_title='log(UMI counts + 1)').show()
        """

        if not show:
            # turn of legends by setting width to 0
            self.legends_width = 0
        else:
            self.color_legend_title = title
            self.legends_width = width

        return self

    def get_axes(self):
        if self.ax_dict is None:
            self.make_figure()
        return self.ax_dict

    def _plot_totals(
        self, total_barplot_ax: Axes, orientation: Literal['top', 'right']
    ):
        """
        Makes the bar plot for totals
        """
        params = self.plot_group_extra
        counts_df = params['counts_df']
        if self.categories_order is not None:
            counts_df = counts_df.loc[self.categories_order]
        if params['color'] is None:
            if f'{self.groupby}_colors' in self.adata.uns:
                color = self.adata.uns[f'{self.groupby}_colors']
            else:
                color = 'salmon'
        else:
            color = params['color']

        if orientation == 'top':
            counts_df.plot(
                kind="bar",
                color=color,
                position=0.5,
                ax=total_barplot_ax,
                edgecolor="black",
                width=0.65,
            )
            # add numbers to the top of the bars
            max_y = max([p.get_height() for p in total_barplot_ax.patches])

            for p in total_barplot_ax.patches:
                p.set_x(p.get_x() + 0.5)
                if p.get_height() >= 1000:
                    display_number = f'{np.round(p.get_height()/1000, decimals=1)}k'
                else:
                    display_number = np.round(p.get_height(), decimals=1)
                total_barplot_ax.annotate(
                    display_number,
                    (p.get_x() + p.get_width() / 2.0, (p.get_height() + max_y * 0.05)),
                    ha="center",
                    va="top",
                    xytext=(0, 10),
                    fontsize="x-small",
                    textcoords="offset points",
                )
            # for k in total_barplot_ax.spines.keys():
            #     total_barplot_ax.spines[k].set_visible(False)
            total_barplot_ax.set_ylim(0, max_y * 1.4)

        elif orientation == 'right':
            counts_df.plot(
                kind="barh",
                color=color,
                position=-0.3,
                ax=total_barplot_ax,
                edgecolor="black",
                width=0.65,
            )

            # add numbers to the right of the bars
            max_x = max([p.get_width() for p in total_barplot_ax.patches])
            for p in total_barplot_ax.patches:
                if p.get_width() >= 1000:
                    display_number = f'{np.round(p.get_width()/1000, decimals=1)}k'
                else:
                    display_number = np.round(p.get_width(), decimals=1)
                total_barplot_ax.annotate(
                    display_number,
                    ((p.get_width()), p.get_y() + p.get_height()),
                    ha="center",
                    va="top",
                    xytext=(10, 10),
                    fontsize="x-small",
                    textcoords="offset points",
                )
            total_barplot_ax.set_xlim(0, max_x * 1.4)

        total_barplot_ax.grid(False)
        total_barplot_ax.axis("off")

    def _plot_colorbar(self, color_legend_ax: Axes, normalize):
        """
        Plots a horizontal colorbar given the ax an normalize values

        Parameters
        ----------
        color_legend_ax
        normalize

        Returns
        -------
        None, updates color_legend_ax

        """
        cmap = pl.get_cmap(self.cmap)

        import matplotlib.colorbar
        from matplotlib.cm import ScalarMappable

        mappable = ScalarMappable(norm=normalize, cmap=cmap)

        matplotlib.colorbar.Colorbar(
            color_legend_ax, mappable=mappable, orientation='horizontal'
        )

        color_legend_ax.set_title(self.color_legend_title, fontsize='small')

        color_legend_ax.xaxis.set_tick_params(labelsize='small')

    def _plot_legend(self, legend_ax, return_ax_dict, normalize):

        # to maintain the fixed height size of the legends, a
        # spacer of variable height is added at top and bottom.
        # The structure for the legends is:
        # first row: variable space to keep the first rows of the same size
        # second row: size legend

        legend_height = self.min_figure_height * 0.08
        height_ratios = [
            self.height - legend_height,
            legend_height,
        ]
        fig, legend_gs = make_grid_spec(
            legend_ax, nrows=2, ncols=1, height_ratios=height_ratios
        )

        color_legend_ax = fig.add_subplot(legend_gs[1])

        self._plot_colorbar(color_legend_ax, normalize)
        return_ax_dict['color_legend_ax'] = color_legend_ax

    def _mainplot(self, ax):
        y_labels = self.categories
        x_labels = self.var_names

        if self.var_names_idx_order is not None:
            x_labels = [x_labels[x] for x in self.var_names_idx_order]

        if self.categories_order is not None:
            y_labels = self.categories_order

        if self.are_axes_swapped:
            x_labels, y_labels = y_labels, x_labels
            ax.set_xlabel(self.groupby)
        else:
            ax.set_ylabel(self.groupby)

        y_ticks = np.arange(len(y_labels)) + 0.5
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels)

        x_ticks = np.arange(len(x_labels)) + 0.5
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, rotation=90, ha='center', minor=False)

        ax.tick_params(axis='both', labelsize='small')
        ax.grid(False)

        # to be consistent with the heatmap plot, is better to
        # invert the order of the y-axis, such that the first group is on
        # top
        ax.set_ylim(len(y_labels), 0)
        ax.set_xlim(0, len(x_labels))

        return check_colornorm(
            self.vboundnorm.vmin,
            self.vboundnorm.vmax,
            self.vboundnorm.vcenter,
            self.vboundnorm.norm,
        )

    def make_figure(self):
        """
        Renders the image but does not call :func:`matplotlib.pyplot.show`. Useful
        when several plots are put together into one figure.

        See also
        --------
        `show()`: Renders and shows the plot.
        `savefig()`: Saves the plot.

        Examples
        --------

        >>> import matplotlib.pyplot as plt
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        >>> fig, (ax0, ax1) = plt.subplots(1, 2)
        >>> sc.pl.MatrixPlot(adata, markers, groupby='bulk_labels', ax=ax0)\
        ...               .style(cmap='Blues', edge_color='none').make_figure()
        >>> sc.pl.DotPlot(adata, markers, groupby='bulk_labels', ax=ax1).make_figure()
        """

        category_height = self.DEFAULT_CATEGORY_HEIGHT
        category_width = self.DEFAULT_CATEGORY_WIDTH

        if self.height is None:
            mainplot_height = len(self.categories) * category_height
            mainplot_width = (
                len(self.var_names) * category_width + self.group_extra_size
            )
            if self.are_axes_swapped:
                mainplot_height, mainplot_width = mainplot_width, mainplot_height

            height = mainplot_height + 1  # +1 for labels

            # if the number of categories is small use
            # a larger height, otherwise the legends do not fit
            self.height = max([self.min_figure_height, height])
            self.width = mainplot_width + self.legends_width
        else:
            self.min_figure_height = self.height
            mainplot_height = self.height

            mainplot_width = self.width - (self.legends_width + self.group_extra_size)

        return_ax_dict = {}
        # define a layout of 1 rows x 2 columns
        #   first ax is for the main figure.
        #   second ax is to plot legends
        legends_width_spacer = 0.7 / self.width

        self.fig, gs = make_grid_spec(
            self.ax or (self.width, self.height),
            nrows=1,
            ncols=2,
            wspace=legends_width_spacer,
            width_ratios=[mainplot_width + self.group_extra_size, self.legends_width],
        )

        if self.has_var_groups:
            # add some space in case 'brackets' want to be plotted on top of the image
            if self.are_axes_swapped:
                var_groups_height = category_height
            else:
                var_groups_height = category_height / 2

        else:
            var_groups_height = 0

        mainplot_width = mainplot_width - self.group_extra_size
        spacer_height = self.height - var_groups_height - mainplot_height
        if not self.are_axes_swapped:
            height_ratios = [spacer_height, var_groups_height, mainplot_height]
            width_ratios = [mainplot_width, self.group_extra_size]

        else:
            height_ratios = [spacer_height, self.group_extra_size, mainplot_height]
            width_ratios = [mainplot_width, var_groups_height]
            # gridspec is the same but rows and columns are swapped

        if self.fig_title is not None and self.fig_title.strip() != '':
            # for the figure title use the ax that contains
            # all the main graphical elements (main plot, dendrogram etc)
            # otherwise the title may overlay with the figure.
            # also, this puts the title centered on the main figure and not
            # centered between the main figure and the legends
            _ax = self.fig.add_subplot(gs[0, 0])
            _ax.axis('off')
            _ax.set_title(self.fig_title)

        # the main plot is divided into three rows and two columns
        # first row is an spacer that is adjusted in case the
        #           legends need more height than the main plot
        # second row is for brackets (if needed),
        # third row is for mainplot and dendrogram/totals (legend goes in gs[0,1]
        # defined earlier)
        mainplot_gs = gridspec.GridSpecFromSubplotSpec(
            nrows=3,
            ncols=2,
            wspace=self.wspace,
            hspace=0.0,
            subplot_spec=gs[0, 0],
            width_ratios=width_ratios,
            height_ratios=height_ratios,
        )
        main_ax = self.fig.add_subplot(mainplot_gs[2, 0])
        return_ax_dict['mainplot_ax'] = main_ax
        if not self.are_axes_swapped:
            if self.plot_group_extra is not None:
                group_extra_ax = self.fig.add_subplot(mainplot_gs[2, 1], sharey=main_ax)
                group_extra_orientation = 'right'
            if self.has_var_groups:
                gene_groups_ax = self.fig.add_subplot(mainplot_gs[1, 0], sharex=main_ax)
                var_group_orientation = 'top'
        else:
            if self.plot_group_extra:
                group_extra_ax = self.fig.add_subplot(mainplot_gs[1, 0], sharex=main_ax)
                group_extra_orientation = 'top'
            if self.has_var_groups:
                gene_groups_ax = self.fig.add_subplot(mainplot_gs[2, 1], sharey=main_ax)
                var_group_orientation = 'right'

        if self.plot_group_extra is not None:
            if self.plot_group_extra['kind'] == 'dendrogram':
                _plot_dendrogram(
                    group_extra_ax,
                    self.adata,
                    self.groupby,
                    dendrogram_key=self.plot_group_extra['dendrogram_key'],
                    ticks=self.plot_group_extra['dendrogram_ticks'],
                    orientation=group_extra_orientation,
                )
            if self.plot_group_extra['kind'] == 'group_totals':
                self._plot_totals(group_extra_ax, group_extra_orientation)

            return_ax_dict['group_extra_ax'] = group_extra_ax

        # plot group legends on top or left of main_ax (if given)
        if self.has_var_groups:
            self._plot_var_groups_brackets(
                gene_groups_ax,
                group_positions=self.var_group_positions,
                group_labels=self.var_group_labels,
                rotation=self.var_group_rotation,
                left_adjustment=0.2,
                right_adjustment=0.7,
                orientation=var_group_orientation,
            )
            return_ax_dict['gene_group_ax'] = gene_groups_ax

        # plot the mainplot
        normalize = self._mainplot(main_ax)

        # code from pandas.plot in add_totals adds
        # minor ticks that need to be removed
        main_ax.yaxis.set_tick_params(which='minor', left=False, right=False)
        main_ax.xaxis.set_tick_params(which='minor', top=False, bottom=False, length=0)
        main_ax.set_zorder(100)
        if self.legends_width > 0:
            legend_ax = self.fig.add_subplot(gs[0, 1])
            self._plot_legend(legend_ax, return_ax_dict, normalize)

        self.ax_dict = return_ax_dict

    def show(self, return_axes: Optional[bool] = None):
        """
        Show the figure

        Parameters
        ----------
        return_axes
             If true return a dictionary with the figure axes. When return_axes is true
             then :func:`matplotlib.pyplot.show` is not called.

        Returns
        -------
        If `return_axes=True`: Dict of :class:`matplotlib.axes.Axes`. The dict key
        indicates the type of ax (eg. `mainplot_ax`)

        See also
        --------
        `render()`: Renders the plot but does not call :func:`matplotlib.pyplot.show`
        `savefig()`: Saves the plot.

        Examples
        -------
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        >>> sc.pl.Plot(adata, markers, groupby='bulk_labels').show()

        """

        self.make_figure()

        if return_axes:
            return self.ax_dict
        else:
            pl.show()

    def savefig(self, filename: str, bbox_inches: Optional[str] = 'tight', **kwargs):
        """
        Save the current figure

        Parameters
        ----------
        filename
            Figure filename. Figure *format* is taken from the file ending unless
            the parameter `format` is given.
        bbox_inches
            By default is set to 'tight' to avoid cropping of the legends.
        kwargs
            Passed to :func:`matplotlib.pyplot.savefig`

        See also
        --------
        `render()`: Renders the plot but does not call :func:`matplotlib.pyplot.show`
        `show()`: Renders and shows the plot

        Examples
        -------
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        >>> sc.pl.BasePlot(adata, markers, groupby='bulk_labels').savefig('plot.pdf')

        """
        self.make_figure()
        pl.savefig(filename, bbox_inches=bbox_inches, **kwargs)

    def _reorder_categories_after_dendrogram(self, dendrogram):
        """\
        Function used by plotting functions that need to reorder the the groupby
        observations based on the dendrogram results.

        The function checks if a dendrogram has already been precomputed.
        If not, `sc.tl.dendrogram` is run with default parameters.

        The results found in `.uns[dendrogram_key]` are used to reorder
        `var_group_labels` and `var_group_positions`.


        Returns
        -------
        None internally updates
        'categories_idx_ordered', 'var_group_names_idx_ordered',
        'var_group_labels' and 'var_group_positions'
        """

        def _format_first_three_categories(_categories):
            """used to clean up warning message"""
            _categories = list(_categories)
            if len(_categories) > 3:
                _categories = _categories[:3] + ['etc.']
            return ', '.join(_categories)

        key = _get_dendrogram_key(self.adata, dendrogram, self.groupby)

        dendro_info = self.adata.uns[key]
        if self.groupby != dendro_info['groupby']:
            raise ValueError(
                "Incompatible observations. The precomputed dendrogram contains "
                f"information for the observation: '{self.groupby}' while the plot is "
                f"made for the observation: '{dendro_info['groupby']}. "
                "Please run `sc.tl.dendrogram` using the right observation.'"
            )

        # order of groupby categories
        categories_idx_ordered = dendro_info['categories_idx_ordered']
        categories_ordered = dendro_info['categories_ordered']

        if len(self.categories) != len(categories_idx_ordered):
            raise ValueError(
                "Incompatible observations. Dendrogram data has "
                f"{len(categories_idx_ordered)} categories but current groupby "
                f"observation {self.groupby!r} contains {len(self.categories)} categories. "
                "Most likely the underlying groupby observation changed after the "
                "initial computation of `sc.tl.dendrogram`. "
                "Please run `sc.tl.dendrogram` again.'"
            )

        # reorder var_groups (if any)
        if self.var_names is not None:
            var_names_idx_ordered = list(range(len(self.var_names)))

        if self.has_var_groups:
            if set(self.var_group_labels) == set(self.categories):
                positions_ordered = []
                labels_ordered = []
                position_start = 0
                var_names_idx_ordered = []
                for cat_name in categories_ordered:
                    idx = self.var_group_labels.index(cat_name)
                    position = self.var_group_positions[idx]
                    _var_names = self.var_names[position[0] : position[1] + 1]
                    var_names_idx_ordered.extend(range(position[0], position[1] + 1))
                    positions_ordered.append(
                        (position_start, position_start + len(_var_names) - 1)
                    )
                    position_start += len(_var_names)
                    labels_ordered.append(self.var_group_labels[idx])
                self.var_group_labels = labels_ordered
                self.var_group_positions = positions_ordered
            else:
                logg.warning(
                    "Groups are not reordered because the `groupby` categories "
                    "and the `var_group_labels` are different.\n"
                    f"categories: {_format_first_three_categories(self.categories)}\n"
                    "var_group_labels: "
                    f"{_format_first_three_categories(self.var_group_labels)}"
                )

        if var_names_idx_ordered is not None:
            var_names_ordered = [self.var_names[x] for x in var_names_idx_ordered]
        else:
            var_names_ordered = None

        self.categories_idx_ordered = categories_idx_ordered
        self.categories_order = dendro_info['categories_ordered']
        self.var_names_idx_order = var_names_idx_ordered
        self.var_names_ordered = var_names_ordered

    @staticmethod
    def _plot_var_groups_brackets(
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
            for idx, (left_coor, right_coor) in enumerate(zip(left, right)):
                verts.append((left_coor, 0))  # lower-left
                verts.append((left_coor, 0.6))  # upper-left
                verts.append((right_coor, 0.6))  # upper-right
                verts.append((right_coor, 0))  # lower-right

                codes.append(Path.MOVETO)
                codes.append(Path.LINETO)
                codes.append(Path.LINETO)
                codes.append(Path.LINETO)

                group_x_center = left[idx] + float(right[idx] - left[idx]) / 2
                gene_groups_ax.text(
                    group_x_center,
                    1.1,
                    group_labels[idx],
                    ha='center',
                    va='bottom',
                    rotation=rotation,
                )
        else:
            top = left
            bottom = right
            for idx, (top_coor, bottom_coor) in enumerate(zip(top, bottom)):
                verts.append((0, top_coor))  # upper-left
                verts.append((0.4, top_coor))  # upper-right
                verts.append((0.4, bottom_coor))  # lower-right
                verts.append((0, bottom_coor))  # lower-left

                codes.append(Path.MOVETO)
                codes.append(Path.LINETO)
                codes.append(Path.LINETO)
                codes.append(Path.LINETO)

                diff = bottom[idx] - top[idx]
                group_y_center = top[idx] + float(diff) / 2
                if diff * 2 < len(group_labels[idx]):
                    # cut label to fit available space
                    group_labels[idx] = group_labels[idx][: int(diff * 2)] + "."
                gene_groups_ax.text(
                    1.1,
                    group_y_center,
                    group_labels[idx],
                    ha='right',
                    va='center',
                    rotation=270,
                    fontsize='small',
                )

        path = Path(verts, codes)

        patch = patches.PathPatch(path, facecolor='none', lw=1.5)

        gene_groups_ax.add_patch(patch)
        gene_groups_ax.grid(False)
        gene_groups_ax.axis('off')
        # remove y ticks
        gene_groups_ax.tick_params(axis='y', left=False, labelleft=False)
        # remove x ticks and labels
        gene_groups_ax.tick_params(
            axis='x', bottom=False, labelbottom=False, labeltop=False
        )

    def _update_var_groups(self):
        """
        checks if var_names is a dict. Is this is the cases, then set the
        correct values for var_group_labels and var_group_positions

        updates var_names, var_group_labels, var_group_positions

        Returns
        -------
        None

        """
        if isinstance(self.var_names, cabc.Mapping):
            if self.has_var_groups:
                logg.warning(
                    "`var_names` is a dictionary. This will reset the current "
                    "values of `var_group_labels` and `var_group_positions`."
                )
            var_group_labels = []
            _var_names = []
            var_group_positions = []
            start = 0
            for label, vars_list in self.var_names.items():
                if isinstance(vars_list, str):
                    vars_list = [vars_list]
                # use list() in case var_list is a numpy array or pandas series
                _var_names.extend(list(vars_list))
                var_group_labels.append(label)
                var_group_positions.append((start, start + len(vars_list) - 1))
                start += len(vars_list)
            self.var_names = _var_names
            self.var_group_labels = var_group_labels
            self.var_group_positions = var_group_positions
            self.has_var_groups = True

        elif isinstance(self.var_names, str):
            self.var_names = [self.var_names]
