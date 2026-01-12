from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import InitVar, dataclass, field
from functools import singledispatch
from typing import TYPE_CHECKING, ClassVar, overload

import numpy as np
import pandas as pd

from .._compat import CSRBase, DaskArray, SpBase, fullname, warn
from .._utils import NeighborsView

if TYPE_CHECKING:
    from typing import NoReturn

    from anndata import AnnData
    from numpy.typing import NDArray

    type _Vals = NDArray | SpBase | DaskArray | pd.DataFrame | pd.Series


__all__ = ["_SparseMetric", "_get_graph"]


@dataclass
class _SparseMetric(ABC):
    graph: CSRBase
    vals: InitVar[_Vals]
    _vals: NDArray | CSRBase | DaskArray = field(init=False)

    name: ClassVar[str]

    def __post_init__(self, vals: _Vals) -> None:
        assert isinstance(type(self).name, str)
        assert self.graph.shape[0] == self.graph.shape[1], (
            "`g` should be a square adjacency matrix"
        )
        self.graph = self.graph.astype(np.float64, copy=False)
        self._vals = _resolve_vals(vals)

    @abstractmethod
    def mtx(self, vals_het: NDArray | CSRBase, /) -> NDArray:
        """Calculate metric when ``._vals`` is a 2D matrix (on an easier to handle version of ``._vals``)."""

    @abstractmethod
    def vec(self) -> np.float64:
        """Calculate metric when ``._vals`` is a 1D vector."""

    def __call__(self) -> np.ndarray:
        match self._vals, self._vals.ndim:
            case _, 2 if isinstance(self._vals, CSRBase | np.ndarray):
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


def _get_graph(
    adata: AnnData,
    *,
    use_graph: str | None = None,
    neighbors_key: str | None = None,
) -> CSRBase:
    if use_graph is not None:
        if neighbors_key is not None:
            msg = "Cannot specify both `use_graph` and `neighbors_key`."
            raise TypeError(msg)
        return adata.obsp[use_graph]
    nv = NeighborsView(adata, neighbors_key)
    if "connectivities" not in nv:
        msg = "Must run neighbors first."
        raise ValueError(msg)
    return nv["connectivities"]


@overload
def _resolve_vals[T: NDArray | DaskArray](val: T) -> T: ...
@overload
def _resolve_vals(val: SpBase) -> CSRBase: ...
@overload
def _resolve_vals(val: pd.DataFrame | pd.Series) -> NDArray: ...


@singledispatch
def _resolve_vals(val: object) -> NoReturn:
    msg = f"Unsupported type {type(val)}"
    raise TypeError(msg)


@_resolve_vals.register(np.ndarray)
@_resolve_vals.register(CSRBase)
@_resolve_vals.register(DaskArray)
def _(
    val: np.ndarray | CSRBase | DaskArray,
) -> np.ndarray | CSRBase | DaskArray:
    return val


@_resolve_vals.register(SpBase)
def _(val: SpBase) -> CSRBase:
    if TYPE_CHECKING:
        from scipy.sparse._base import _spbase

        assert isinstance(val, _spbase)
    return val.tocsr()


@_resolve_vals.register(pd.DataFrame)
@_resolve_vals.register(pd.Series)
def _(val: pd.DataFrame | pd.Series) -> NDArray:
    return val.to_numpy()


def _vals_heterogeneous[V: NDArray | CSRBase](
    vals: V,
) -> tuple[V, NDArray[np.bool_] | slice, NDArray[np.float64]]:
    """Check that values wont cause issues in computation.

    Returns new set of vals, and indexer to put values back into result.

    For details on why this is neccesary, see:
    https://github.com/scverse/scanpy/issues/1806
    """
    from fast_array_utils.stats import is_constant

    full_result = np.empty(vals.shape[0], dtype=np.float64)
    full_result.fill(np.nan)
    idxer = ~is_constant(vals, axis=1)
    if idxer.all():
        idxer = slice(None)
    else:
        msg = f"{len(idxer) - idxer.sum()} variables were constant, will return nan for these."
        warn(msg, UserWarning)
    return vals[idxer], idxer, full_result
