"""Storing parameters in `.uns`, the way the parameters `use` replaced were stored."""

from __future__ import annotations

from typing import TYPE_CHECKING

from numpy.typing import NDArray

from ._types import AdRef, LayerAcc, MultiAcc
from ._use import _resolve_obs

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from ._types import RepAcc


def _ref_to_uns[M: NDArray | None](ref: AdRef | str | M, resolved: M) -> str | M:
    """Store a vector reference in `.uns` the way its deprecated predecessor was stored.

    Accessors aren't stored in `.uns` (see https://github.com/scverse/anndata/issues/1979),
    so a column of `.obs`/`.var` becomes its name, and any other reference its `resolved` array.
    Plain strings (column names, see `dim_acc`) and arrays are stored unchanged.
    """
    if not isinstance(ref, AdRef):
        return ref
    from anndata.acc import MetaAcc

    return ref.idx if isinstance(ref.acc, MetaAcc) else resolved


def _use_to_uns(
    use: RepAcc | str | None, use_rep: str | None = None, *, rep: bool
) -> dict[str, str]:
    """Store `use` in `.uns` params the way the parameters it replaces were stored.

    Accessors aren't stored in `.uns` (see https://github.com/scverse/anndata/issues/1979),
    so e.g. `A.obsm['pca']` becomes `use_rep='pca'` (or `obsm='pca'` if not `rep`)
    and `A.layers['x']` becomes `layer='x'`.
    The deprecated `use_rep` is stored unchanged.
    """
    if use is None:
        return {} if use_rep is None else dict(use_rep=use_rep)
    match _resolve_obs(use, rep=True):
        case LayerAcc(k=None):
            return dict(use_rep="X") if rep else {}
        case LayerAcc(k=key):
            return dict(layer=key)
        case MultiAcc(k=key):
            return {"use_rep" if rep else "obsm": key}
        case _:  # pragma: no cover
            raise AssertionError


def _use_from_uns(params: Mapping[str, Any]) -> RepAcc | str | None:
    """Read back a `rep`resentation stored by `_use_to_uns(..., rep=True)`.

    Returns a deprecated `use_rep` string (`'X'` or an `.obsm` key),
    an accessor (for representations `use_rep` can’t express), or `None`.
    """
    if (layer := params.get("layer")) is not None:
        from anndata.acc import A

        return A.layers[layer]
    return params.get("use_rep")
