from __future__ import annotations

import inspect
import sys
from contextlib import contextmanager
from enum import IntEnum, StrEnum, auto
from functools import cached_property
from logging import getLevelNamesMapping
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Literal, get_args

from . import logging
from ._compat import old_positionals
from ._singleton import SingletonMeta
from .logging import _RootLogger, _set_log_file, _set_log_level

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable
    from typing import Any, ClassVar, Self, TextIO

    from .preprocessing import _highly_variable_genes as hvg

    # Collected from the print_* functions in matplotlib.backends
    _Format = (
        Literal["png", "jpg", "tif", "tiff"]  # noqa: PYI030
        | Literal["pdf", "ps", "eps", "svg", "svgz", "pgf"]
        | Literal["raw", "rgba"]
    )
    _VerbosityName = Literal["error", "warning", "info", "hint", "debug"]
    _LoggingLevelName = Literal["CRITICAL", "ERROR", "WARNING", "INFO", "HINT", "DEBUG"]


AnnDataFileFormat = Literal["h5ad", "zarr"]


class Preset(StrEnum):
    ScanpyV1 = auto()
    SeuratV5 = auto()

    @cached_property
    def highly_variable_genes(self) -> hvg.Flavor:
        match self:
            case Preset.ScanpyV1:
                return "seurat"
            case Preset.SeuratV5:
                return "seurat_v3"


_VERBOSITY_TO_LOGLEVEL: dict[int | _VerbosityName, _LoggingLevelName] = {
    "error": "ERROR",
    "warning": "WARNING",
    "info": "INFO",
    "hint": "HINT",
    "debug": "DEBUG",
}
_VERBOSITY_TO_LOGLEVEL.update(dict(enumerate(list(_VERBOSITY_TO_LOGLEVEL.values()))))


class Verbosity(IntEnum):
    """Logging verbosity levels."""

    error = 0
    warning = 1
    info = 2
    hint = 3
    debug = 4

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Verbosity):
            return self is other
        if isinstance(other, int):
            return self.value == other
        if isinstance(other, str):
            return self.name == other
        return NotImplemented

    @property
    def level(self) -> int:
        m = getLevelNamesMapping()
        return m[_VERBOSITY_TO_LOGLEVEL[self.name]]

    @contextmanager
    def override(
        self, verbosity: Verbosity | _VerbosityName | int
    ) -> Generator[Verbosity, None, None]:
        """Temporarily override verbosity."""
        settings.verbosity = verbosity
        yield self
        settings.verbosity = self


# backwards compat
Verbosity.warn = Verbosity.warning


def _type_check(var: Any, varname: str, types: type | tuple[type, ...]) -> None:
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = f"{', '.join(type_names[:-1])} or {type_names[-1]}"
    msg = f"{varname} must be of type {possible_types_str}"
    raise TypeError(msg)


class SettingsMeta(SingletonMeta):
    _preset: Preset
    # logging
    _root_logger: _RootLogger
    _logfile: TextIO | None
    _verbosity: Verbosity
    # rest
    N_PCS: int
    """Default number of principal components to use."""
    _plot_suffix: str
    _file_format_data: AnnDataFileFormat
    _file_format_figs: str
    _autosave: bool
    _autoshow: bool
    _writedir: Path
    _cachedir: Path
    _datasetdir: Path
    _figdir: Path
    _cache_compression: str | None
    _max_memory: float
    _n_jobs: int
    _categories_to_ignore: list[str]
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
    def preset(cls) -> Preset:
        """Preset to use."""
        return cls._preset

    @preset.setter
    def preset(cls, preset: Preset | str) -> None:
        cls._preset = Preset(preset)

    @property
    def verbosity(cls) -> Verbosity:
        """Verbosity level (default `warning`).

        Level 0: only show 'error' messages.
        Level 1: also show 'warning' messages.
        Level 2: also show 'info' messages.
        Level 3: also show 'hint' messages.
        Level 4: also show very detailed progress for 'debug'ging.
        """
        return cls._verbosity

    @verbosity.setter
    def verbosity(cls, verbosity: Verbosity | _VerbosityName | int) -> None:
        verbosity_str_options: list[_VerbosityName] = [
            v for v in _VERBOSITY_TO_LOGLEVEL if isinstance(v, str)
        ]
        if isinstance(verbosity, Verbosity):
            cls._verbosity = verbosity
        elif isinstance(verbosity, int):
            cls._verbosity = Verbosity(verbosity)
        elif isinstance(verbosity, str):
            verbosity = verbosity.lower()
            if verbosity not in verbosity_str_options:
                msg = (
                    f"Cannot set verbosity to {verbosity}. "
                    f"Accepted string values are: {verbosity_str_options}"
                )
                raise ValueError(msg)
            cls._verbosity = Verbosity(verbosity_str_options.index(verbosity))
        else:
            _type_check(verbosity, "verbosity", (str, int))
        _set_log_level(cls, _VERBOSITY_TO_LOGLEVEL[cls._verbosity.name])

    @property
    def plot_suffix(cls) -> str:
        """Global suffix that is appended to figure filenames."""
        return cls._plot_suffix

    @plot_suffix.setter
    def plot_suffix(cls, plot_suffix: str) -> None:
        _type_check(plot_suffix, "plot_suffix", str)
        cls._plot_suffix = plot_suffix

    @property
    def file_format_data(cls) -> AnnDataFileFormat:
        """File format for saving AnnData objects."""
        return cls._file_format_data

    @file_format_data.setter
    def file_format_data(cls, file_format: AnnDataFileFormat) -> None:
        _type_check(file_format, "file_format_data", str)
        if file_format not in (file_format_options := get_args(AnnDataFileFormat)):
            msg = (
                f"Cannot set file_format_data to {file_format}. "
                f"Must be one of {file_format_options}"
            )
            raise ValueError(msg)
        cls._file_format_data: AnnDataFileFormat = file_format

    @property
    def file_format_figs(cls) -> str:
        """File format for saving figures.

        For example `'png'`, `'pdf'` or `'svg'`. Many other formats work as well (see
        `matplotlib.pyplot.savefig`).
        """
        return cls._file_format_figs

    @file_format_figs.setter
    def file_format_figs(self, figure_format: str) -> None:
        _type_check(figure_format, "figure_format_data", str)
        self._file_format_figs = figure_format

    @property
    def autosave(cls) -> bool:
        """Automatically save figures in :attr:`~scanpy.settings.figdir` (default `False`).

        Do not show plots/figures interactively.
        """
        return cls._autosave

    @autosave.setter
    def autosave(cls, autosave: bool) -> None:
        _type_check(autosave, "autosave", bool)
        cls._autosave = autosave

    @property
    def autoshow(cls) -> bool:
        """Automatically show figures if `autosave == False` (default `True`).

        There is no need to call the matplotlib pl.show() in this case.
        """
        return cls._autoshow

    @autoshow.setter
    def autoshow(cls, autoshow: bool) -> None:
        _type_check(autoshow, "autoshow", bool)
        cls._autoshow = autoshow

    @property
    def writedir(cls) -> Path:
        """Directory where the function scanpy.write writes to by default."""
        return cls._writedir

    @writedir.setter
    def writedir(cls, writedir: Path | str) -> None:
        _type_check(writedir, "writedir", (str, Path))
        cls._writedir = Path(writedir)

    @property
    def cachedir(cls) -> Path:
        """Directory for cache files (default `'./cache/'`)."""
        return cls._cachedir

    @cachedir.setter
    def cachedir(cls, cachedir: Path | str) -> None:
        _type_check(cachedir, "cachedir", (str, Path))
        cls._cachedir = Path(cachedir)

    @property
    def datasetdir(cls) -> Path:
        """Directory for example :mod:`~scanpy.datasets` (default `'./data/'`)."""
        return cls._datasetdir

    @datasetdir.setter
    def datasetdir(cls, datasetdir: Path | str) -> None:
        _type_check(datasetdir, "datasetdir", (str, Path))
        cls._datasetdir = Path(datasetdir).resolve()

    @property
    def figdir(cls) -> Path:
        """Directory for saving figures (default `'./figures/'`)."""
        return cls._figdir

    @figdir.setter
    def figdir(cls, figdir: Path | str) -> None:
        _type_check(figdir, "figdir", (str, Path))
        cls._figdir = Path(figdir)

    @property
    def cache_compression(cls) -> str | None:
        """Compression for `sc.read(..., cache=True)` (default `'lzf'`).

        May be `'lzf'`, `'gzip'`, or `None`.
        """
        return cls._cache_compression

    @cache_compression.setter
    def cache_compression(cls, cache_compression: str | None) -> None:
        if cache_compression not in {"lzf", "gzip", None}:
            msg = (
                f"`cache_compression` ({cache_compression}) "
                "must be in {'lzf', 'gzip', None}"
            )
            raise ValueError(msg)
        cls._cache_compression = cache_compression

    @property
    def max_memory(cls) -> int | float:
        """Maximum memory usage in Gigabyte.

        Is currently not well respected…
        """
        return cls._max_memory

    @max_memory.setter
    def max_memory(cls, max_memory: float) -> None:
        _type_check(max_memory, "max_memory", (int, float))
        cls._max_memory = max_memory

    @property
    def n_jobs(cls) -> int:
        """Default number of jobs/ CPUs to use for parallel computing.

        Set to `-1` in order to use all available cores.
        Not all algorithms support special behavior for numbers < `-1`,
        so make sure to leave this setting as >= `-1`.
        """
        return cls._n_jobs

    @n_jobs.setter
    def n_jobs(cls, n_jobs: int) -> None:
        _type_check(n_jobs, "n_jobs", int)
        cls._n_jobs = n_jobs

    @property
    def logpath(cls) -> Path | None:
        """The file path `logfile` was set to."""
        return cls._logpath

    @logpath.setter
    def logpath(cls, logpath: Path | str | None) -> None:
        _type_check(logpath, "logfile", (str, Path))
        if logpath is None:
            cls._logfile = None
            cls._logpath = None
            return
        # set via “file object” branch of logfile.setter
        cls.logfile = Path(logpath).open("a")  # noqa: SIM115
        cls._logpath = Path(logpath)

    @property
    def logfile(cls) -> TextIO | None:
        """The open file to write logs to.

        Set it to a :class:`~pathlib.Path` or :class:`str` to open a new one.
        The default `None` corresponds to :obj:`sys.stdout` in jupyter notebooks
        and to :obj:`sys.stderr` otherwise.

        For backwards compatibility, setting it to `''` behaves like setting it to `None`.
        """
        return cls._logfile

    @logfile.setter
    def logfile(cls, logfile: Path | str | TextIO | None) -> None:
        if not hasattr(logfile, "write") and logfile:
            cls.logpath = logfile
        else:  # file object
            if not logfile:  # None or ''
                logfile = sys.stdout if cls._is_run_from_ipython() else sys.stderr
            cls._logfile = logfile
            cls._logpath = None
            _set_log_file(cls)

    @property
    def categories_to_ignore(cls) -> list[str]:
        """Categories that are omitted in plotting etc."""
        return cls._categories_to_ignore

    @categories_to_ignore.setter
    def categories_to_ignore(cls, categories_to_ignore: Iterable[str]) -> None:
        categories_to_ignore = list(categories_to_ignore)
        for i, cat in enumerate(categories_to_ignore):
            _type_check(cat, f"categories_to_ignore[{i}]", str)
        cls._categories_to_ignore = categories_to_ignore

    # --------------------------------------------------------------------------------
    # Functions
    # --------------------------------------------------------------------------------

    @old_positionals(
        "scanpy",
        "dpi",
        "dpi_save",
        "frameon",
        "vector_friendly",
        "fontsize",
        "figsize",
        "color_map",
        "format",
        "facecolor",
        "transparent",
        "ipython_format",
    )
    def set_figure_params(  # noqa: PLR0913
        cls,
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
            Resolution of saved figures. This should typically be higher to achieve
            publication quality.
        frameon
            Add frames and axes labels to scatter plots.
        vector_friendly
            Plot scatter plots using `png` backend even when exporting as `pdf` or `svg`.
        fontsize
            Set the fontsize for several `rcParams` entries. Ignored if `scanpy=False`.
        figsize
            Set plt.rcParams['figure.figsize'].
        color_map
            Convenience method for setting the default color map. Ignored if `scanpy=False`.
        format
            This sets the default format for saving figures: `file_format_figs`.
        facecolor
            Sets backgrounds via `rcParams['figure.facecolor'] = facecolor` and
            `rcParams['axes.facecolor'] = facecolor`.
        transparent
            Save figures with transparent back ground. Sets
            `rcParams['savefig.transparent']`.
        ipython_format
            Only concerns the notebook/IPython environment; see
            `matplotlib_inline.backend_inline.set_matplotlib_formats
            <https://github.com/ipython/matplotlib-inline/blob/b93777db35267acefe6e37d14214360362d2e8b2/matplotlib_inline/backend_inline.py#L280-L281>`_
            for details.

        """
        if cls._is_run_from_ipython():
            # No docs yet: https://github.com/ipython/matplotlib-inline/issues/12
            from matplotlib_inline.backend_inline import set_matplotlib_formats

            if isinstance(ipython_format, str):
                ipython_format = [ipython_format]

            set_matplotlib_formats(*ipython_format)

        from matplotlib import rcParams

        cls._vector_friendly = vector_friendly
        cls.file_format_figs = format
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
            from .plotting._rcmod import set_rcParams_scanpy

            set_rcParams_scanpy(fontsize=fontsize, color_map=color_map)
        if figsize is not None:
            rcParams["figure.figsize"] = figsize
        cls._frameon = frameon

    @staticmethod
    def _is_run_from_ipython() -> bool:
        """Determine whether we're currently in IPython."""
        import builtins

        return getattr(builtins, "__IPYTHON__", False)

    def __str__(cls) -> str:
        return "\n".join(
            f"{k} = {v!r}"
            for k, v in inspect.getmembers(cls)
            if not k.startswith("_") and k != "getdoc"
        )


class settings(metaclass=SettingsMeta):
    """Settings for scanpy."""

    def __new__(cls) -> type[Self]:
        return cls

    _preset = Preset.ScanpyV1
    # logging
    _root_logger: ClassVar = _RootLogger(logging.INFO)
    _logfile: ClassVar = None
    _logpath: ClassVar = None
    _verbosity: ClassVar = Verbosity.warning
    # rest
    N_PCS: ClassVar = 50
    _plot_suffix: ClassVar = ""
    _file_format_data: ClassVar = "h5ad"
    _file_format_figs: ClassVar = "pdf"
    _autosave: ClassVar = False
    _autoshow: ClassVar = True
    _writedir: ClassVar = Path("./write")
    _cachedir: ClassVar = Path("./cache")
    _datasetdir: ClassVar = Path("./data")
    _figdir: ClassVar = Path("./figures")
    _cache_compression: ClassVar = "lzf"
    _max_memory: ClassVar = 15
    _n_jobs: ClassVar = 1
    _categories_to_ignore: ClassVar = ["N/A", "dontknow", "no_gate", "?"]
    _frameon: ClassVar = True
    _vector_friendly: ClassVar = False
    _low_resolution_warning: ClassVar = True
    _start: ClassVar = time()
    _previous_time: ClassVar = _start
    _previous_memory_usage: ClassVar = -1
