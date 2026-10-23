from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping
from dataclasses import InitVar, dataclass, field
from functools import singledispatch
from typing import TYPE_CHECKING, ClassVar, overload

import numpy as np
import pandas as pd
from scverse_misc import Deprecation

from .._compat import CSRBase, DaskArray, SpBase, fullname, warn
from .._docs import doc_use
from .._utils import NeighborsView
from ..get._types import AdRef, GraphAcc
from ..get._use import _get_arr, _get_vec, _resolve_obs

if TYPE_CHECKING:
    from typing import Literal, NoReturn

    from anndata import AnnData
    from numpy.typing import NDArray

    from ..get._types import ArrAcc, VecRef

    type _Vals = NDArray | SpBase | DaskArray | pd.DataFrame | pd.Series
    type _MetricUse = (
        ArrAcc | VecRef | Mapping[Literal["graph", "vals"], ArrAcc | VecRef]
    )


__all__ = ["DEPR_GRAPH", "DOC_USE", "_SparseMetric", "_graph_and_vals"]


DOC_USE = doc_use(
    "Which values to compute on, e.g. `A.obsm['X_pca']` or the 1D `A.obs['n_genes']`.\n"
    "    Can also be a mapping with the keys `'vals'` (the values) and `'graph'` (an `.obsp` graph,\n"
    "    by default the connectivities :func:`~scanpy.pp.neighbors` wrote),\n"
    "    e.g. `use=dict(graph=A.obsp['foo_connectivities'], vals=A.obsm['X_pca'])`.",
    legacy=("layer", "obsm", "obsp"),
)
DEPR_GRAPH = Deprecation("1.13.0", "Use `use=dict(graph=A.obsp[...])` instead.")


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


def _graph_and_vals(
    adata: AnnData,
    use: _MetricUse | None,
    *,
    vals: _Vals | None,
    use_graph: str | None,
    neighbors_key: str | None,
    use_raw: bool,
    layer: str | None,
    obsm: str | None,
    obsp: str | None,
) -> tuple[CSRBase, _Vals]:
    """Resolve the graph and values of `gearys_c`/`morans_i` called on an `AnnData`."""
    from anndata.acc import A

    if not isinstance(use, Mapping):
        use = {} if use is None else dict(vals=use)
    if unknown := use.keys() - {"graph", "vals"}:
        msg = f"Unknown key(s) {sorted(unknown)} in `use`, valid are ['graph', 'vals']"
        raise TypeError(msg)

    if (graph := use.get("graph")) is None:
        g = _get_graph(adata, use_graph=use_graph, neighbors_key=neighbors_key)
    elif use_graph is not None or neighbors_key is not None:
        msg = "Cannot specify both `use['graph']` and `use_graph`/`neighbors_key`."
        raise TypeError(msg)
    elif not isinstance(graph := _resolve_obs(graph), GraphAcc):
        msg = f"`use['graph']` must point to `.obsp`, was {graph!r}"
        raise TypeError(msg)
    else:
        g = adata[graph]

    # not `@deprecated_arg`, which can’t compare positionally passed arrays to `None`
    if vals is not None:
        warn(
            "The argument vals is deprecated and will be removed in the future. "
            "Use `use` instead.",
            FutureWarning,
        )
        if "vals" in use:
            msg = "Cannot specify both `vals` and `use`."
            raise TypeError(msg)
        return g, vals
    if isinstance(src := use.get("vals"), str):
        src = A.resolve(src)
    if isinstance(src, AdRef):
        if use_raw or layer is not None or obsm is not None or obsp is not None:
            msg = "`use` cannot be combined with `layer`/`use_raw`/`obsm`/`obsp`"
            raise TypeError(msg)
        # e.g. `.obs` columns are pandas arrays
        return g, np.asarray(_get_vec(adata, src, dim="obs"))
    return g, _get_arr(
        adata, _resolve_obs(src), use_raw=use_raw, layer=layer, obsm=obsm, obsp=obsp
    ).T


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
) -> tuple[V, NDArray[np.bool] | slice, NDArray[np.float64]]:
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
