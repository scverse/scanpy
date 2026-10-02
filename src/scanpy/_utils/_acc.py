"""Helpers to turn string specs into :mod:`anndata.acc` dimensions/accessors."""

from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING, cast, overload

if TYPE_CHECKING:
    from typing import Literal

    from anndata.acc import GraphAcc, LayerAcc, MultiAcc
    from hv_anndata import AdDim


__all__ = ["VecAcc", "_resolve", "_resolve_all", "_resolve_some"]


type VecAcc = LayerAcc | MultiAcc | GraphAcc
"""Accessor for a whole container, e.g. ``A.X``, ``A.layers["counts"]``, or ``A.obsp["distances"]``."""


@overload
def _resolve(spec: str | AdDim, /, *, vec: Literal[True] = True) -> AdDim: ...
@overload
def _resolve(spec: str | VecAcc, /, *, vec: Literal[False]) -> VecAcc: ...
def _resolve(spec: str | AdDim | VecAcc, /, *, vec: bool = True) -> AdDim | VecAcc:
    """Resolve a string spec (see :meth:`~anndata.acc.AdAcc.resolve`), pass refs through.

    ``vec=True`` yields a dimension (e.g. from ``"obs.a"`` or ``"X[:,g]"``),
    ``vec=False`` a container accessor (e.g. from ``"layers.counts"``).
    """
    from anndata.acc import GraphAcc, LayerAcc, MultiAcc
    from hv_anndata import A, AdDim

    if isinstance(spec, str):
        return A.resolve(spec, vec=vec)
    if not isinstance(spec, AdDim if vec else LayerAcc | MultiAcc | GraphAcc):
        what = "an AdDim" if vec else "a container accessor"
        msg = f"Expected a string or {what}, got {spec!r}."
        raise TypeError(msg)
    return spec


@overload
def _resolve_all(
    specs: Collection[str | AdDim], /, *, vec: Literal[True] = True
) -> list[AdDim]: ...
@overload
def _resolve_all(
    specs: Collection[str | VecAcc], /, *, vec: Literal[False]
) -> list[VecAcc]: ...
def _resolve_all(
    specs: Collection[str | AdDim | VecAcc], /, *, vec: bool = True
) -> list[AdDim] | list[VecAcc]:
    """Resolve a mixed collection of string specs and refs/accessors. See :func:`_resolve`."""
    if isinstance(specs, str) or not isinstance(specs, Collection):
        msg = (
            f"Expected a collection of specs, got {specs!r}. Did you mean [{specs!r}]?"
        )
        raise TypeError(msg)
    if vec:
        return [_resolve(cast("str | AdDim", spec)) for spec in specs]
    return [_resolve(cast("str | VecAcc", spec), vec=False) for spec in specs]


@overload
def _resolve_some(spec: str | AdDim, /) -> AdDim: ...
@overload
def _resolve_some(spec: Collection[str | AdDim], /) -> list[AdDim]: ...
def _resolve_some(
    spec: str | AdDim | Collection[str | AdDim], /
) -> AdDim | list[AdDim]:
    """Resolve a single spec or a collection of them. See :func:`_resolve`."""
    from hv_anndata import AdDim

    return _resolve(spec) if isinstance(spec, str | AdDim) else _resolve_all(spec)
