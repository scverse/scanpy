"""Set the default matplotlib.rcParams."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl
from cycler import cycler
from matplotlib import rcParams

from . import palettes

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    # Collected from the print_* functions in matplotlib.backends
    type _Format = (
        Literal["png", "jpg", "tif", "tiff"]  # noqa: PYI030
        | Literal["pdf", "ps", "eps", "svg", "svgz", "pgf"]
        | Literal["raw", "rgba"]
    )

# Module-level state (used by plotting functions throughout the legacy module)
FRAMEON: bool = True
VECTOR_FRIENDLY: bool = False


def _is_run_from_ipython() -> bool:
    import builtins

    return getattr(builtins, "__IPYTHON__", False)


def set_figure_params(  # noqa: PLR0913
    *,
    scanpy: bool = True,
    dpi: int = 80,
    dpi_save: int = 150,
    frameon: bool = True,
    vector_friendly: bool = True,
    fontsize: int = 14,
    figsize: int | None = None,
    color_map: str | None = None,
    format: _Format = "pdf",
    facecolor: str | None = None,
    transparent: bool = False,
    ipython_format: str | Iterable[str] = "retina",
) -> None:
    """Set resolution/size, styling and format of figures.

    Parameters
    ----------
    scanpy
        Init default values for :obj:`matplotlib.rcParams` suited for Scanpy.
    dpi
        Resolution of rendered figures – this influences the size of figures in notebooks.
    dpi_save
        Resolution of saved figures.
        This should typically be higher to achieve publication quality.
    frameon
        Add frames and axes labels to scatter plots.
    vector_friendly
        Plot scatter plots using `png` backend even when exporting as `pdf` or `svg`.
    fontsize
        Set the fontsize for several `rcParams` entries. Ignored if `scanpy=False`.
    figsize
        Set `rcParams['figure.figsize']`.
    color_map
        Convenience method for setting the default color map. Ignored if `scanpy=False`.
    format
        This sets the default format for saving figures: `file_format_figs`.
    facecolor
        Sets backgrounds via `rcParams['figure.facecolor'] = facecolor` and
        `rcParams['axes.facecolor'] = facecolor`.
    transparent
        Save figures with transparent background.
        Sets `rcParams['savefig.transparent']`.
    ipython_format
        Only concerns the notebook/IPython environment; see
        `matplotlib_inline.backend_inline.set_matplotlib_formats
        <https://github.com/ipython/matplotlib-inline/blob/b93777db35267acefe6e37d14214360362d2e8b2/matplotlib_inline/backend_inline.py#L280-L281>`_
        for details.

    """
    global FRAMEON, VECTOR_FRIENDLY  # noqa: PLW0603

    if _is_run_from_ipython():
        # No docs yet: https://github.com/ipython/matplotlib-inline/issues/12
        from matplotlib_inline.backend_inline import set_matplotlib_formats

        if isinstance(ipython_format, str):
            ipython_format = [ipython_format]

        set_matplotlib_formats(*ipython_format)

    VECTOR_FRIENDLY = vector_friendly

    from scanpy._settings import settings

    settings.file_format_figs = format

    if dpi is not None:
        rcParams["figure.dpi"] = dpi
    if dpi_save is not None:
        rcParams["savefig.dpi"] = dpi_save
    if transparent is not None:
        rcParams["savefig.transparent"] = transparent
    if facecolor is not None:
        rcParams["figure.facecolor"] = facecolor
        rcParams["axes.facecolor"] = facecolor
    if scanpy:
        set_rcParams_scanpy(fontsize=fontsize, color_map=color_map)
    if figsize is not None:
        rcParams["figure.figsize"] = figsize

    FRAMEON = frameon


def set_rcParams_scanpy(fontsize=14, color_map=None) -> None:  # noqa: N802
    """Set matplotlib.rcParams to Scanpy defaults.

    Call this through :func:`scanpy.set_figure_params`.
    """
    # figure
    rcParams["figure.figsize"] = (4, 4)
    rcParams["figure.subplot.left"] = 0.18
    rcParams["figure.subplot.right"] = 0.96
    rcParams["figure.subplot.bottom"] = 0.15
    rcParams["figure.subplot.top"] = 0.91

    rcParams["lines.linewidth"] = 1.5  # the line width of the frame
    rcParams["lines.markersize"] = 6
    rcParams["lines.markeredgewidth"] = 1

    # font
    rcParams["font.sans-serif"] = [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ]
    rcParams["font.size"] = fontsize
    rcParams["legend.fontsize"] = 0.92 * fontsize
    rcParams["axes.titlesize"] = fontsize
    rcParams["axes.labelsize"] = fontsize

    # legend
    rcParams["legend.numpoints"] = 1
    rcParams["legend.scatterpoints"] = 1
    rcParams["legend.handlelength"] = 0.5
    rcParams["legend.handletextpad"] = 0.4

    # color cycle
    rcParams["axes.prop_cycle"] = cycler(color=palettes.default_20)

    # lines
    rcParams["axes.linewidth"] = 0.8
    rcParams["axes.edgecolor"] = "black"
    rcParams["axes.facecolor"] = "white"

    # ticks
    rcParams["xtick.color"] = "k"
    rcParams["ytick.color"] = "k"
    rcParams["xtick.labelsize"] = fontsize
    rcParams["ytick.labelsize"] = fontsize

    # axes grid
    rcParams["axes.grid"] = True
    rcParams["grid.color"] = ".8"

    # color map
    rcParams["image.cmap"] = rcParams["image.cmap"] if color_map is None else color_map


def set_rcParams_defaults() -> None:  # noqa: N802
    """Reset `matplotlib.rcParams` to defaults."""
    rcParams.update(mpl.rcParamsDefault)
