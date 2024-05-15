"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pytest

from .data import (
    _pbmc3ks_parametrized_session,
    backed_adata,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

__all__ = [
    "float_dtype",
    "_doctest_env",
    "_pbmc3ks_parametrized_session",
    "pbmc3k_parametrized",
    "pbmc3k_parametrized_small",
    "backed_adata",
]


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param


@pytest.fixture()
def _doctest_env(cache: pytest.Cache, tmp_path: Path) -> Generator[None, None, None]:
    from scanpy import settings
    from scanpy._compat import chdir

    showwarning_orig = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None):  # noqa: PLR0917
        if file is None:
            if line is None:
                import linecache

                line = linecache.getline(filename, lineno)
            line = line.strip()
            print(f"{category.__name__}: {message}\n    {line}")
        else:
            showwarning_orig(message, category, filename, lineno, file, line)

    # make errors visible and the rest ignored
    warnings.filters = [
        ("default", *rest) for action, *rest in warnings.filters if action == "error"
    ] + [("ignore", None, Warning, None, 0)]

    warnings.showwarning = showwarning
    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
    warnings.showwarning = showwarning_orig
    settings.datasetdir = old_dd
