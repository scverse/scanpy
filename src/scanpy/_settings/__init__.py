from __future__ import annotations

import sys
from collections.abc import Container
from itertools import chain
from pathlib import Path
from time import time
from typing import TYPE_CHECKING, Annotated, Literal, Protocol, runtime_checkable

import scverse_misc
from pydantic import (
    AfterValidator,
    PrivateAttr,
    computed_field,
    field_validator,
    model_validator,
)
from scverse_misc import Deprecation, deprecated

from .._compat import set_module
from ..logging import _RootLogger, _set_log_file, _set_log_level
from .presets import Default, Preset
from .verbosity import Verbosity

if sys.version_info >= (3, 14):
    from io import Writer
else:

    @set_module("io")
    @runtime_checkable
    class Writer[T: str | bytes](Protocol):
        def write(self, data: T, /) -> int: ...


if TYPE_CHECKING:
    from typing import Self, TextIO

    from pydantic import ValidationInfo

    from .verbosity import _VerbosityName


__all__ = ["AnnDataFileFormat", "Default", "Preset", "Verbosity"]

AnnDataFileFormat = Literal["h5ad", "zarr"]


def _default_logfile() -> TextIO:
    return sys.stdout if _is_run_from_ipython() else sys.stderr


def _is_run_from_ipython() -> bool:
    """Determine whether we're currently in IPython."""
    import builtins

    return getattr(builtins, "__IPYTHON__", False)


class Settings(scverse_misc.Settings):
    def model_post_init(self, context: object) -> None:
        _set_log_level(self)
        _set_log_file(self)

    preset: Annotated[Preset, AfterValidator(Preset.check)] = Preset.ScanpyV1
    """Preset to use."""

    # logging
    verbosity: Verbosity = Verbosity.warning
    """Verbosity level (default :attr:`Verbosity.warning`)."""

    _root_logger: Annotated[
        _RootLogger, PrivateAttr(default=_RootLogger(Verbosity.warning.level))
    ]

    logfile: Writer[str] = _default_logfile()
    """The open file to write logs to.

    Set it to a :class:`~pathlib.Path` or :class:`str` to open a new one.
    By default uses :obj:`sys.stdout` in jupyter notebooks and :obj:`sys.stderr` otherwise.
    """

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

    n_jobs: int = 4
    """Default number of jobs/ CPUs to use for parallel computing.

    Set to `-1` in order to use all available cores.
    Not all algorithms support special behavior for numbers < `-1`,
    so make sure to leave this setting as >= `-1`.

    .. version-changed:: 1.13.0
        The default changed from `1` to `4`.
    """

    categories_to_ignore: Container[str] = frozenset({
        "N/A",
        "dontknow",
        "no_gate",
        "?",
    })
    """Categories that are omitted in plotting etc."""

    _start: Annotated[float, PrivateAttr(default_factory=time)]
    """Time when the settings module is first imported."""
    _previous_time: Annotated[float, PrivateAttr(default_factory=time)]
    """Variable for timing program parts."""
    _previous_memory_usage: Annotated[int, PrivateAttr(default=-1)]
    """Stores the previous memory usage."""

    @computed_field
    @property
    def logpath(self) -> Path | None:
        """The file path `logfile` was set to."""
        if self.logfile is _default_logfile():
            return None
        if (name := getattr(self.logfile, "name", None)) is None:
            return None
        return Path(name)

    @logpath.setter
    def logpath(self, path: Path | None) -> None:
        self.logfile = _default_logfile() if path is None else path.open("a")

    @field_validator("verbosity", mode="before")
    @classmethod
    def _check_verbosity(
        cls, verbosity: Verbosity | _VerbosityName | int, /
    ) -> Verbosity:
        """Lenient conversion of verbosity from `int` or level name."""
        try:
            return (
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

    @model_validator(mode="after")
    def _logging_side_effects(self, info: ValidationInfo) -> Self:
        """Side effect of setting the verbosity."""
        match info.field_name:
            case "verbosity":
                _set_log_level(self)
            case "logfile":
                _set_log_file(self)
        return self

    # --------------------------------------------------------------------------------
    # Functions
    # --------------------------------------------------------------------------------

    @deprecated(Deprecation("1.11.3", "Use :func:`scanpy.set_figure_params` instead"))
    def set_figure_params(self, *args, **kwargs) -> None:
        from ..plotting.legacy.mpl_settings import set_figure_params

        set_figure_params(*args, **kwargs)

    def __str__(self) -> str:
        return "\n".join(
            f"{k} = {getattr(self, k)!r}"
            for k in chain(type(self).model_fields, type(self).model_computed_fields)
        )

    def __hash__(self) -> int:
        return hash((id(self),))


settings = Settings()
