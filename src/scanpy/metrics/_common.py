from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import InitVar, dataclass, field
from functools import singledispatch
from typing import TYPE_CHECKING, ClassVar, overload

import numpy as np
import pandas as pd
from scipy import sparse

from .._compat import DaskArray, fullname

if TYPE_CHECKING:
    from typing import NoReturn, TypeVar

    from anndata import AnnData
    from numpy.typing import NDArray

    T_NonSparse = TypeVar("T_NonSparse", bound=NDArray | DaskArray)
    V = TypeVar("V", bound=NDArray | sparse.csr_matrix)

    _Vals = NDArray | sparse.spmatrix | DaskArray | pd.DataFrame | pd.Series


__all__ = ["_get_graph", "_SparseMetric"]


@dataclass
class _SparseMetric(ABC):
    graph: sparse.csr_matrix
    vals: InitVar[_Vals]
    _vals: NDArray | sparse.csr_matrix | DaskArray = field(init=False)

    name: ClassVar[str]

    def __post_init__(self, vals: _Vals) -> None:
        assert isinstance(type(self).name, str)
        assert self.graph.shape[0] == self.graph.shape[1], (
            "`g` should be a square adjacency matrix"
        )
        self.graph = self.graph.astype(np.float64, copy=False)
        self._vals = _resolve_vals(vals)

    @abstractmethod
    def mtx(self, vals_het: NDArray | sparse.csr_matrix, /) -> NDArray:
        """Calculate metric when ``._vals`` is a 2D matrix (on an easier to handle version of ``._vals``)."""

    @abstractmethod
    def vec(self) -> np.float64:
        """Calculate metric when ``._vals`` is a 1D vector."""

    def __call__(self) -> np.ndarray:
        match self._vals, self._vals.ndim:
            case sparse.csr_matrix() | np.ndarray(), 2:
                assert self.graph.shape[0] == self._vals.shape[1]
                vals_het, idxer, full_result = _vals_heterogeneous(self._vals)
                result = self.mtx(vals_het.astype(np.float64, copy=False))
                full_result[idxer] = result
                return full_result
            case np.ndarray(), 1:
                assert self.graph.shape[0] == self._vals.shape[0]
                return self.vec()
            case _, _:
                msg = (
                    f"{self.name} metric not implemented for vals of type "
                    f"{fullname(type(self._vals))} and ndim {self._vals.ndim}."
                )
                raise NotImplementedError(msg)


def _get_graph(adata: AnnData, *, use_graph: str | None = None) -> sparse.csr_matrix:
    if use_graph is not None:
        raise NotImplementedError()
    # Fix for anndata<0.7
    if hasattr(adata, "obsp") and "connectivities" in adata.obsp:
        return adata.obsp["connectivities"]
    elif "neighbors" in adata.uns:
        return adata.uns["neighbors"]["connectivities"]
    else:
        msg = "Must run neighbors first."
        raise ValueError(msg)


@overload
def _resolve_vals(val: T_NonSparse) -> T_NonSparse: ...
@overload
def _resolve_vals(val: sparse.spmatrix) -> sparse.csr_matrix: ...
@overload
def _resolve_vals(val: pd.DataFrame | pd.Series) -> NDArray: ...


@singledispatch
def _resolve_vals(val: object) -> NoReturn:
    msg = f"Unsupported type {type(val)}"
    raise TypeError(msg)


@_resolve_vals.register(np.ndarray)
@_resolve_vals.register(sparse.csr_matrix)
@_resolve_vals.register(DaskArray)
def _(
    val: np.ndarray | sparse.csr_matrix | DaskArray,
) -> np.ndarray | sparse.csr_matrix | DaskArray:
    return val


@_resolve_vals.register(sparse.spmatrix)
def _(val: sparse.spmatrix) -> sparse.csr_matrix:
    return sparse.csr_matrix(val)


@_resolve_vals.register(pd.DataFrame)
@_resolve_vals.register(pd.Series)
def _(val: pd.DataFrame | pd.Series) -> NDArray:
    return val.to_numpy()


def _vals_heterogeneous(
    vals: V,
) -> tuple[V, NDArray[np.bool_] | slice, NDArray[np.float64]]:
    """Check that values wont cause issues in computation.

    Returns new set of vals, and indexer to put values back into result.

    For details on why this is neccesary, see:
    https://github.com/scverse/scanpy/issues/1806
    """
    from scanpy._utils import is_constant

    full_result = np.empty(vals.shape[0], dtype=np.float64)
    full_result.fill(np.nan)
    idxer = ~is_constant(vals, axis=1)
    if idxer.all():
        idxer = slice(None)
    else:
        warnings.warn(
            UserWarning(
                f"{len(idxer) - idxer.sum()} variables were constant, will return nan for these.",
            )
        )
    return vals[idxer], idxer, full_result
