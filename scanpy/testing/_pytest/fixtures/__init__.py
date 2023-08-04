"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
from __future__ import annotations

from typing import TYPE_CHECKING
from collections.abc import Callable

import pytest
import numpy as np
from numpy.typing import ArrayLike
from scipy import sparse
from anndata.tests.helpers import asarray

if TYPE_CHECKING:
    from dask.array import Array as DaskArray

from ..._pytest.marks import needs
from .data import (
    _pbmc3ks_parametrized_session,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)


__all__ = [
    'array_type',
    'float_dtype',
    '_pbmc3ks_parametrized_session',
    'pbmc3k_parametrized',
    'pbmc3k_parametrized_small',
]


def _as_dask_array(x: ArrayLike) -> DaskArray:
    import dask.array as da

    return da.asarray(x)


@pytest.fixture(
    params=[
        pytest.param(asarray, id='numpy-ndarray'),
        pytest.param(sparse.csr_matrix, id='scipy-csr'),
        pytest.param(sparse.csc_matrix, id='scipy-csc'),
        pytest.param(_as_dask_array, marks=[needs('dask')], id='dask-array'),
    ]
)
def array_type(
    request,
) -> Callable[[ArrayLike], DaskArray | np.ndarray | sparse.spmatrix]:
    """Function which converts passed array to one of the common array types."""
    return request.param


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param
