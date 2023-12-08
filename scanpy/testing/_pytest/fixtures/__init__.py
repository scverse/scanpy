"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
from __future__ import annotations

import sys
import warnings
from typing import TYPE_CHECKING

import numpy as np
import pytest

from .data import (
    _pbmc3ks_parametrized_session,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

__all__ = [
    "float_dtype",
    "doctest_env",
    "_pbmc3ks_parametrized_session",
    "pbmc3k_parametrized",
    "pbmc3k_parametrized_small",
]


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param


@pytest.fixture()
def doctest_env(cache: pytest.Cache, tmp_path: Path) -> Generator[None, None, None]:
    from scanpy import settings
    from scanpy._compat import chdir

    showwarning_orig = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None):
        if file is None:
            file = sys.stdout
        showwarning_orig(message, category, filename, lineno, file, line)

    warnings.filterwarnings("default", category=UserWarning)

    warnings.showwarning = showwarning
    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
    warnings.showwarning = showwarning_orig
    settings.datasetdir = old_dd
