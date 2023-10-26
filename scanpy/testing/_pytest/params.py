"""Like fixtures, but more flexible"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pytest
from scipy import sparse
from anndata.tests.helpers import asarray, as_dense_dask_array, as_sparse_dask_array

from .._pytest.marks import needs

if TYPE_CHECKING:
    from _pytest.mark.structures import ParameterSet


_AT_MAP: dict[
    tuple[Literal["mem", "dask"], Literal["dense", "sparse"]],
    tuple[ParameterSet, ...],
] = {
    ("mem", "dense"): (pytest.param(asarray, id="numpy_ndarray"),),
    ("mem", "sparse"): (
        pytest.param(sparse.csr_matrix, id="scipy_csr"),
        pytest.param(sparse.csc_matrix, id="scipy_csc"),
    ),
    ("dask", "dense"): (
        pytest.param(as_dense_dask_array, marks=[needs("dask")], id="dask_array_dense"),
    ),
    ("dask", "sparse"): (
        pytest.param(
            as_sparse_dask_array, marks=[needs("dask")], id="dask_array_sparse"
        ),
        # probably not necessary to also do csc
    ),
}

ARRAY_TYPES_MEM = tuple(
    at for (strg, _), ats in _AT_MAP.items() if strg == "mem" for at in ats
)
ARRAY_TYPES_DASK = tuple(
    at for (strg, _), ats in _AT_MAP.items() if strg == "dask" for at in ats
)

ARRAY_TYPES_DENSE = tuple(
    at for (_, spsty), ats in _AT_MAP.items() if spsty == "dense" for at in ats
)
ARRAY_TYPES_SPARSE = tuple(
    at for (_, spsty), ats in _AT_MAP.items() if spsty == "dense" for at in ats
)

ARRAY_TYPES_SUPPORTED = tuple(
    (
        pytest.param(
            *at.values,
            marks=[pytest.mark.xfail(reason="sparse-in-dask not supported"), *at.marks],
            id=at.id,
        )
        if attrs == ("dask", "sparse")
        else at
    )
    for attrs, ats in _AT_MAP.items()
    for at in ats
)
"""
Sparse matrices in dask arrays aren’t officially supported upstream,
so add xfail to them.
"""

ARRAY_TYPES = tuple(at for ats in _AT_MAP.values() for at in ats)
