from __future__ import annotations

import warnings
from functools import wraps

import anndata as ad
from packaging import version

from .._settings import settings


def check_datasetdir_exists(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        settings.datasetdir.mkdir(exist_ok=True)
        return f(*args, **kwargs)

    return wrapper


def filter_oldformatwarning(f):
    """
    Filters anndata.OldFormatWarning from being thrown by the wrapped function.
    """

    @wraps(f)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            if version.parse(ad.__version__).release >= (0, 8):
                warnings.filterwarnings(
                    "ignore", category=ad.OldFormatWarning, module="anndata"
                )
            return f(*args, **kwargs)

    return wrapper
