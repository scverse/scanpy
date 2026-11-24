"""Accessor types (with fallbacks for anndata without `anndata.acc`) and type aliases."""

from __future__ import annotations

from collections.abc import Collection
from importlib.util import find_spec
from typing import TYPE_CHECKING, TypedDict

import numpy as np
from anndata import AnnData
from numpy.typing import NDArray

if TYPE_CHECKING:
    from typing import TypeIs

    from anndata.acc import Idx2D


if TYPE_CHECKING or find_spec("anndata.acc"):
    from anndata.acc import AdRef, GraphAcc, LayerAcc, MultiAcc, RefAcc
else:
    AdRef = type("AdRef", (), dict(__module__="anndata.acc"))
    RefAcc = type("RefAcc", (), dict(__module__="anndata.acc"))
    GraphAcc = type("GraphAcc", (), dict(__module__="anndata.acc"))
    LayerAcc = type("LayerAcc", (), dict(__module__="anndata.acc"))
    MultiAcc = type("MultiAcc", (), dict(__module__="anndata.acc"))


type VecRef = AdRef[Idx2D | int, AnnData] | str
"""A reference to a 1D array, e.g. `A.obs['x']` or `'obs.x'`."""


type Mask = NDArray[np.bool] | VecRef
"""A boolean array, or a reference to one (see `_check_mask`)."""


def _collection_of[T](
    thing: object, typ: type[T] | tuple[type[T], ...]
) -> TypeIs[Collection[T]]:
    return (
        isinstance(thing, Collection)
        and not isinstance(thing, typ)
        and len(thing) > 0
        and all(isinstance(e, typ) for e in thing)
    )


class _Rep(TypedDict, total=False):  # noqa: PYI049
    use_raw: bool
    layer: str | None
    obsm: str | None
    obsp: str | None
    varm: str | None
    varp: str | None


class _SetRep(TypedDict, total=False):  # noqa: PYI049
    """Like `_Rep`, but for writing: no `use_raw`, plus 1D destinations."""

    layer: str | None
    obsm: str | None
    obsp: str | None
    varm: str | None
    varp: str | None
    obs: str | None
    var: str | None


type ArrAcc = GraphAcc | LayerAcc | MultiAcc


type RepAcc = LayerAcc | MultiAcc
"""Accessor usable as a representation (`use_rep`), i.e. a non-graph 2D array."""
