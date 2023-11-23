"""Like fixtures, but more flexible"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pytest
from anndata.tests.helpers import as_dense_dask_array, as_sparse_dask_array, asarray
from scipy import sparse

from .._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Iterable

    from _pytest.mark.structures import ParameterSet


def param_with(
    at: ParameterSet,
    *,
    marks: Iterable[pytest.Mark | pytest.MarkDecorator] = (),
    id: str | None = None,
) -> ParameterSet:
    return pytest.param(*at.values, marks=[*at.marks, *marks], id=id or at.id)


MAP_ARRAY_TYPES: dict[
    tuple[Literal["mem", "dask"], Literal["dense", "sparse"]],
    tuple[ParameterSet, ...],
] = {
    ("mem", "dense"): (pytest.param(asarray, id="numpy_ndarray"),),
    ("mem", "sparse"): (
        pytest.param(sparse.csr_matrix, id="scipy_csr"),
        pytest.param(sparse.csc_matrix, id="scipy_csc"),
    ),
    ("dask", "dense"): (
        pytest.param(as_dense_dask_array, marks=[needs.dask], id="dask_array_dense"),
    ),
    ("dask", "sparse"): (
        pytest.param(as_sparse_dask_array, marks=[needs.dask], id="dask_array_sparse"),
        # probably not necessary to also do csc
    ),
}

ARRAY_TYPES_MEM = tuple(
    at for (strg, _), ats in MAP_ARRAY_TYPES.items() if strg == "mem" for at in ats
)
ARRAY_TYPES_DASK = tuple(
    at for (strg, _), ats in MAP_ARRAY_TYPES.items() if strg == "dask" for at in ats
)

ARRAY_TYPES_DENSE = tuple(
    at for (_, spsty), ats in MAP_ARRAY_TYPES.items() if spsty == "dense" for at in ats
)
ARRAY_TYPES_SPARSE = tuple(
    at for (_, spsty), ats in MAP_ARRAY_TYPES.items() if spsty == "dense" for at in ats
)

ARRAY_TYPES_SUPPORTED = tuple(
    (
        param_with(at, marks=[pytest.mark.xfail(reason="sparse-in-dask not supported")])
        if attrs == ("dask", "sparse")
        else at
    )
    for attrs, ats in MAP_ARRAY_TYPES.items()
    for at in ats
)
"""
Sparse matrices in dask arrays aren’t officially supported upstream,
so add xfail to them.
"""

ARRAY_TYPES = tuple(at for ats in MAP_ARRAY_TYPES.values() for at in ats)
