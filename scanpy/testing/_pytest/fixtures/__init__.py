"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
from __future__ import annotations

from collections.abc import Callable

import pytest
import numpy as np
from numpy.typing import ArrayLike
from scipy import sparse
from anndata.tests.helpers import asarray

from ...._compat import DaskArray
from ..._pytest.marks import needs
from .data import (
    _pbmc3ks_parametrized_session,
    pbmc3k_parametrized,
    pbmc3k_parametrized_small,
)


__all__ = [
    "array_type",
    "float_dtype",
    "_pbmc3ks_parametrized_session",
    "pbmc3k_parametrized",
    "pbmc3k_parametrized_small",
]


def _as_dense_da(x: ArrayLike) -> DaskArray:
    if isinstance(x, DaskArray):
        return x

    import dask.array as da

    x = asarray(x)
    quarter_shape = tuple(np.ceil(np.array(x.shape) / 2))
    return da.from_array(x, chunks=quarter_shape)


def _as_csr_da(x: ArrayLike) -> DaskArray:
    return _as_dense_da(x).map_blocks(sparse.csr_matrix)


@pytest.fixture(
    params=[
        pytest.param(asarray, id="numpy_ndarray"),
        pytest.param(sparse.csr_matrix, id="scipy_csr"),
        pytest.param(sparse.csc_matrix, id="scipy_csc"),
        pytest.param(_as_dense_da, marks=[needs("dask")], id="dask_array_dense"),
        pytest.param(_as_csr_da, marks=[needs("dask")], id="dask_array_sparse"),
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
