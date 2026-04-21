from __future__ import annotations

import sys
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Annotated, Literal, TextIO

import pydantic_settings
from pydantic import AfterValidator

from .. import logging
from .._compat import deprecated
from ..logging import _RootLogger, _set_log_file, _set_log_level
from .presets import Default, Preset
from .verbosity import Verbosity

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .verbosity import _VerbosityName

    # Collected from the print_* functions in matplotlib.backends
    type _Format = (
        Literal["png", "jpg", "tif", "tiff"]  # noqa: PYI030
        | Literal["pdf", "ps", "eps", "svg", "svgz", "pgf"]
        | Literal["raw", "rgba"]
    )


__all__ = ["AnnDataFileFormat", "Default", "Preset", "Verbosity"]

AnnDataFileFormat = Literal["h5ad", "zarr"]


def _default_logfile() -> TextIO:
    return sys.stdout if _is_run_from_ipython() else sys.stderr


def _is_run_from_ipython() -> bool:
    """Determine whether we're currently in IPython."""
    import builtins

    return getattr(builtins, "__IPYTHON__", False)


# `type` is only here because of https://github.com/astral-sh/ruff/issues/20225
class Settings(pydantic_settings.BaseSettings):
    def model_post_init(self, context: object) -> None:
        # logging
        self._verbosity = Verbosity.warning
        self._root_logger = _RootLogger(logging.WARNING)
        self._logfile = _default_logfile()
        self._logpath = None
        _set_log_level(self)
        _set_log_file(self)

        # figure
        self._frameon = True
        self._vector_friendly = False
        self._low_resolution_warning = True
        self._start = self._previous_time = time()
        self._previous_memory_usage = -1

    preset: Annotated[Preset, AfterValidator(Preset.check)] = Preset.ScanpyV1
    """Preset to use."""

    # logging
    _verbosity: Verbosity
    _root_logger: _RootLogger
    _logfile: TextIO
    _logpath: Path | None

    # rest
    N_PCS: int = 50
    """Default number of principal components to use."""

    plot_suffix: str = ""
    """Global suffix that is appended to figure filenames."""

    file_format_data: AnnDataFileFormat = "h5ad"
    """File format for saving AnnData objects."""

    file_format_figs: str = "pdf"
    """File format for saving figures.

    For example `'png'`, `'pdf'` or `'svg'`. Many other formats work as well (see
    :func:`matplotlib.pyplot.savefig`).
    """

    autosave: bool = False
    """Automatically save figures in :attr:`~scanpy.settings.figdir` (default `False`).

    Do not show plots/figures interactively.
    """

    autoshow: bool = True
    """Automatically show figures if `autosave == False` (default `True`).

    There is no need to call the matplotlib pl.show() in this case.
    """

    writedir: Path = Path("./write")
    """Directory where the function scanpy.write writes to by default."""

    cachedir: Path = Path("./cache")
    """Directory for cache files (default `'./cache/'`)."""

    datasetdir: Annotated[Path, AfterValidator(Path.resolve)] = Path("./data")
    """Directory for example :mod:`~scanpy.datasets` (default `'./data/'`)."""

    figdir: Path = Path("./figures")
    r"""Directory for `autosave`\ ing figures (default `'./figures/'`)."""

    cache_compression: Literal["lzf", "gzip"] | None = "lzf"
    """Compression for `sc.read(..., cache=True)` (default `'lzf'`)."""

    max_memory: float = 15.0
    """Maximum memory usage in Gigabyte.

    Is currently not well respected…
    """

    n_jobs: int = 1
    """Default number of jobs/ CPUs to use for parallel computing.

    Set to `-1` in order to use all available cores.
    Not all algorithms support special behavior for numbers < `-1`,
    so make sure to leave this setting as >= `-1`.
    """

    categories_to_ignore: list[str] = ["N/A", "dontknow", "no_gate", "?"]
    """Categories that are omitted in plotting etc."""

    _frameon: bool
    """See set_figure_params."""
    _vector_friendly: bool
    """Set to true if you want to include pngs in svgs and pdfs."""
    _low_resolution_warning: bool
    """Print warning when saving a figure with low resolution."""
    _start: float
    """Time when the settings module is first imported."""
    _previous_time: float
    """Variable for timing program parts."""
    _previous_memory_usage: int
    """Stores the previous memory usage."""

    @property
    def verbosity(self) -> Verbosity:
        """Verbosity level (default :attr:`Verbosity.warning`)."""
        return self._verbosity

    @verbosity.setter
    def _set_verbosity(self, verbosity: Verbosity | _VerbosityName | int) -> None:
        try:
            self._verbosity = (
                Verbosity[verbosity.lower()]
                if isinstance(verbosity, str)
                else Verbosity(verbosity)
            )
        except KeyError:
            msg = (
                f"Cannot set verbosity to {verbosity}. "
                f"Accepted string values are: {Verbosity.__members__.keys()}"
            )
            raise ValueError(msg) from None
        _set_log_level(self)

    @property
    def logpath(self) -> Path | None:
        """The file path `logfile` was set to."""
        return self._logpath

    @logpath.setter
    def logpath(self, logpath: Path | str | None) -> None:
        if logpath is None:
            self.logfile = None
            self._logpath = None
            return
        # set via “file object” branch of logfile.setter
        self.logfile = Path(logpath).open("a")  # noqa: SIM115
        self._logpath = Path(logpath)

    @property
    def logfile(self) -> TextIO:
        """The open file to write logs to.

        Set it to a :class:`~pathlib.Path` or :class:`str` to open a new one.
        The default `None` corresponds to :obj:`sys.stdout` in jupyter notebooks
        and to :obj:`sys.stderr` otherwise.

        For backwards compatibility, setting it to `''` behaves like setting it to `None`.
        """
        return self._logfile

    @logfile.setter
    def logfile(self, logfile: Path | str | TextIO | None) -> None:
        if not logfile:  # "" or None
            logfile = _default_logfile()
        if isinstance(logfile, Path | str):
            self.logpath = logfile
            return
        self._logfile = logfile
        self._logpath = None
        _set_log_file(self)

    # --------------------------------------------------------------------------------
    # Functions
    # --------------------------------------------------------------------------------

    @deprecated("Use `scanpy.set_figure_params` instead")
    def set_figure_params(self, *args, **kwargs) -> None:
        self._set_figure_params(*args, **kwargs)

    def _set_figure_params(  # noqa: PLR0913
        self,
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
        if _is_run_from_ipython():
            # No docs yet: https://github.com/ipython/matplotlib-inline/issues/12
            from matplotlib_inline.backend_inline import set_matplotlib_formats

            if isinstance(ipython_format, str):
                ipython_format = [ipython_format]

            set_matplotlib_formats(*ipython_format)

        from matplotlib import rcParams

        self._vector_friendly = vector_friendly
        self.file_format_figs = format
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
            from ..plotting._rcmod import set_rcParams_scanpy

            set_rcParams_scanpy(fontsize=fontsize, color_map=color_map)
        if figsize is not None:
            rcParams["figure.figsize"] = figsize
        self._frameon = frameon

    def __str__(self) -> str:
        return "\n".join(f"{k} = {getattr(self, k)!r}" for k in type(self).model_fields)

    def __hash__(self) -> int:
        return hash((id(self),))


settings = Settings()
