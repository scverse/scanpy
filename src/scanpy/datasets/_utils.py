from __future__ import annotations

import warnings
from functools import wraps
from typing import TYPE_CHECKING

import anndata as ad
from packaging.version import Version

from .._settings import settings

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import ParamSpec, TypeVar

    P = ParamSpec("P")
    R = TypeVar("R")


def check_datasetdir_exists(f: Callable[P, R]) -> Callable[P, R]:
    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        settings.datasetdir.mkdir(exist_ok=True)
        return f(*args, **kwargs)

    return wrapper


def filter_oldformatwarning(f: Callable[P, R]) -> Callable[P, R]:
    """Filter anndata.OldFormatWarning from being thrown by the wrapped function."""

    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        with warnings.catch_warnings():
            if Version(ad.__version__).release >= (0, 8):
                warnings.filterwarnings(
                    "ignore", category=ad.OldFormatWarning, module="anndata"
                )
            return f(*args, **kwargs)

    return wrapper
