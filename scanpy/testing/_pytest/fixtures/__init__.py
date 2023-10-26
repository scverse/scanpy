"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
from __future__ import annotations

import pytest
import numpy as np
from .data import (
    _pbmc3ks_parametrized_session,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)


__all__ = [
    "float_dtype",
    "_pbmc3ks_parametrized_session",
    "pbmc3k_parametrized",
    "pbmc3k_parametrized_small",
]


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param
