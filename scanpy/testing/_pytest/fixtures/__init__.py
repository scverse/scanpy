"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import numpy as np
from .data import (
    _pbmc3ks_parametrized_session,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)


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
def doctest_env(cache: pytest.Cache, tmp_path: Path) -> None:
    from scanpy import settings
    from scanpy._compat import chdir

    old_dd, settings.datasetdir = settings.datasetdir, cache.mkdir("scanpy-data")
    with chdir(tmp_path):
        yield
    settings.datasetdir = old_dd
