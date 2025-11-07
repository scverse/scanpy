"""Like fixtures, but more flexible."""

from __future__ import annotations

from functools import partial, wraps
from importlib.metadata import version
from typing import TYPE_CHECKING

import pytest
from anndata.tests.helpers import asarray
from packaging.version import Version
from scipy import sparse

from .._helpers import as_dense_dask_array, as_sparse_dask_array
from .._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Any, Literal

    import numpy as np
    from _pytest.mark.structures import ParameterSet

    from ....scanpy._compat import DaskArray


skipif_no_sparray = pytest.mark.skipif(
    Version(version("anndata")) < Version("0.11"),
    reason="scipy cs{rc}_array not supported in anndata<0.11",
)


def param_with(
    at: ParameterSet,
    transform: Callable[..., Iterable[Any]] = lambda x: (x,),
    *,
    marks: Iterable[pytest.Mark | pytest.MarkDecorator] = (),
    id: str | None = None,
) -> ParameterSet:
    return pytest.param(
        *transform(*at.values), marks=[*at.marks, *marks], id=id or at.id
    )


def _chunked_1d(
    f: Callable[[np.ndarray], DaskArray],
) -> Callable[[np.ndarray], DaskArray]:
    @wraps(f)
    def wrapper(a: np.ndarray) -> DaskArray:
        da = f(a)
        return da.rechunk(
            (da.chunksize[0], -1)
            if not hasattr(da._meta, "format") or da._meta.format == "csr"
            else (-1, da.chunksize[1])
        )

    wrapper.__name__ = f"{wrapper.__name__}-1d_chunked"
    return wrapper


MAP_ARRAY_TYPES: dict[
    tuple[Literal["mem", "dask"], Literal["dense", "sparse"]],
    tuple[ParameterSet, ...],
] = {
    ("mem", "dense"): (pytest.param(asarray, id="numpy_ndarray"),),
    ("mem", "sparse"): (
        pytest.param(sparse.csr_matrix, id="scipy_csr_mat"),  # noqa: TID251
        pytest.param(sparse.csc_matrix, id="scipy_csc_mat"),  # noqa: TID251
        pytest.param(sparse.csr_array, id="scipy_csr_arr", marks=[skipif_no_sparray]),  # noqa: TID251
    ),
    ("dask", "dense"): tuple(
        pytest.param(
            wrapper(as_dense_dask_array),
            marks=[needs.dask],
            id=f"dask_array_dense{suffix}",
        )
        for wrapper, suffix in [(lambda x: x, ""), (_chunked_1d, "-1d_chunked")]
    ),
    ("dask", "sparse"): tuple(
        pytest.param(
            wrapper(as_sparse_dask_array),
            marks=[needs.dask],
            id=f"dask_array_sparse{suffix}",
        )
        for wrapper, suffix in [
            (lambda x: x, ""),
            *(
                ((_chunked_1d, "-1d_chunked"),)
                if Version(version("anndata")) < Version("0.12.5")
                else (
                    (
                        lambda func,
                        format=format,
                        matrix_or_array=matrix_or_array: _chunked_1d(
                            partial(
                                func, typ=getattr(sparse, f"{format}_{matrix_or_array}")
                            )
                        ),
                        f"-1d_chunked-{format}_{matrix_or_array}",
                    )
                    for format in ["csr", "csc"]
                    # TODO: use `array` as well once anndata 0.13 drops
                    for matrix_or_array in ["matrix"]
                )
            ),
        ]
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
    at for (_, spsty), ats in MAP_ARRAY_TYPES.items() if "sparse" in spsty for at in ats
)
ARRAY_TYPES_SPARSE_DASK_UNSUPPORTED = tuple(
    (
        param_with(at, marks=[pytest.mark.xfail(reason="sparse-in-dask not supported")])
        if attrs[0] == "dask" and "sparse" in attrs[1]
        else at
    )
    for attrs, ats in MAP_ARRAY_TYPES.items()
    for at in ats
)

ARRAY_TYPES = tuple(at for ats in MAP_ARRAY_TYPES.values() for at in ats)
