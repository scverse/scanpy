from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING

from .._settings import settings

if TYPE_CHECKING:
    from collections.abc import Callable


def check_datasetdir_exists[**P, R](f: Callable[P, R]) -> Callable[P, R]:
    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        settings.datasetdir.mkdir(exist_ok=True)
        return f(*args, **kwargs)

    return wrapper
