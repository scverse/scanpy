"""Resolving and reading inputs: `use` accessors, legacy choices, masks, and vector references."""

from __future__ import annotations

from collections.abc import Collection
from importlib.util import find_spec
from typing import TYPE_CHECKING, overload

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.typing import NDArray

from .._settings import Default, Preset
from .._utils import dim_acc
from ._types import AdRef, ArrAcc, GraphAcc, LayerAcc, MultiAcc, _collection_of

if TYPE_CHECKING:
    from typing import Any, Literal, Unpack

    from anndata.acc import Idx2D

    from .._compat import CSBase, DaskArray
    from ._types import RefAcc, RepAcc, _Rep


@overload
def _get_arr(
    adata: AnnData,
    use: Collection[ArrAcc | str],
    *,
    dim: Literal["obs", "var"] | None = None,
) -> list[Any]: ...


@overload
def _get_arr(
    adata: AnnData,
    use: ArrAcc | str | None = None,
    *,
    dim: Literal["obs", "var"] | None = None,
    **choices: Unpack[_Rep],
) -> Any: ...


def _get_arr(  # noqa: PLR0911, PLR0912
    adata: AnnData,
    use: ArrAcc | str | Collection[ArrAcc | str] | None = None,
    *,
    dim: Literal["obs", "var"] | None = None,
    **choices: Unpack[_Rep],
) -> Any:
    """Get a 2D array aligned with `dim`, via an `anndata.acc` accessor or old-style choices."""
    if _collection_of(use, (GraphAcc, LayerAcc, MultiAcc, str)):
        return [_get_arr(adata, u, dim=dim, **choices) for u in use]

    if use is not None:
        if isinstance(use, str):
            from anndata.acc import A

            use = A.resolve(use, vec=False)

        if any(v not in (None, False) for v in choices.values()):
            msg = "`use` cannot be combined with `layer`/`use_raw`/`obsm`/`obsp`/`varm`/`varp`"
            raise TypeError(msg)
        if not isinstance(use, GraphAcc | LayerAcc | MultiAcc):
            msg = (
                "`use` must be a `LayerAcc` (e.g. `A.X`, `A.layers[...]`), "
                "`GraphAcc` (e.g. `A.obsp[...]`, `A.varp[...]`), or "
                f"`MultiAcc` (e.g. `A.obsm[...]`, `A.varm[...]`), was {use!r}"
            )
            raise TypeError(msg)
        if isinstance(use, MultiAcc | GraphAcc) and dim is not None and dim != use.dim:
            msg = f"`dim` ({dim!r}) does not match `use`'s ({use.dim!r})"
            raise ValueError(msg)
        data = adata[use]
        if isinstance(use, LayerAcc) and dim == "var":
            data = data.T
        return data

    # https://github.com/scverse/scanpy/issues/1546
    if not isinstance(use_raw := choices.get("use_raw", False), bool):
        msg = f"use_raw expected to be bool, was {type(use_raw)}."
        raise TypeError(msg)
    assert choices.keys() <= {"layer", "use_raw", "obsm", "obsp", "varm", "varp"}
    if dim is None:
        dim = "var" if (choices.get("varm") or choices.get("varp")) else "obs"

    match [(k, v) for k, v in choices.items() if v not in {None, False}]:
        case []:
            return adata.X.T if dim == "var" else adata.X
        # can’t use {"key": v} as match expression, since they allow additional entries
        case [("layer", layer)]:
            return adata.layers[layer].T if dim == "var" else adata.layers[layer]
        case [("use_raw", True)]:
            return adata.raw.X
        case [(("obsm" | "obsp") as k, v)]:
            if dim == "var":
                msg = f"`{k}` cannot be used when `dim` is `var`"
                raise ValueError(msg)
            return adata.obsm[v] if k == "obsm" else adata.obsp[v]
        case [(("varm" | "varp") as k, v)]:
            if dim == "obs":
                msg = f"`{k}` cannot be used when `dim` is `obs`"
                raise ValueError(msg)
            return adata.varm[v] if k == "varm" else adata.varp[v]
        case picked:
            valid = [f"`{k}`" for k, _ in picked]
            valid[-1] = f"or {valid[-1]}"
            msg = f"Only one of {', '.join(valid)} can be specified."
            raise ValueError(msg)


def _is_multi(use: ArrAcc | str | None) -> bool:
    """Whether `use` points at a `.{obs,var}m` array, whose columns aren’t variables."""
    if use is None:
        return False
    if isinstance(use, str):
        from anndata.acc import A

        use = A.resolve(use, vec=False)
    return isinstance(use, MultiAcc)


def _mask_arg[M](
    mask: M | Default, legacy: M | None, *, dim: Literal["obs", "var"]
) -> M | Default:
    """Merge the `mask` argument with its deprecated `mask_{dim}` predecessor."""
    if legacy is not None:
        if mask is not None and not isinstance(mask, Default):
            msg = f"Pass either `mask` or `mask_{dim}`, not both."
            raise TypeError(msg)
        return dim_acc(legacy, dim=dim) if isinstance(legacy, str) else legacy
    if isinstance(mask, str):
        if not find_spec("anndata.acc"):
            msg = (
                f"`mask={mask!r}` requires `anndata>=0.13.3`. "
                f"Pass a boolean array instead, or a column name via the deprecated `mask_{dim}`."
            )
            raise ImportError(msg)
        from anndata.acc import A

        return A.resolve(mask, vec=True)
    return mask


def _mask_hvg[M](adata: AnnData, mask: M | Default) -> M | None:
    """Resolve a `Default` `mask` to `.var['highly_variable']` if there is one."""
    if not isinstance(mask, Default):
        return mask
    if "highly_variable" not in adata.var:
        return None
    return dim_acc("highly_variable", dim="var")


def _check_mask[M: NDArray[np.bool] | NDArray[np.floating] | pd.Series | None](
    data: AnnData | np.ndarray | CSBase | DaskArray,
    mask: str | AdRef[Idx2D | int, AnnData] | M,
    dim: Literal["obs", "var"],
    *,
    allow_probabilities: bool = False,
) -> M:  # Could also be a series, but should be one or the other
    """Validate mask argument.

    Params
    ------
    data
        Annotated data matrix or numpy array.
    mask
        Mask (or probabilities if `allow_probabilities=True`).
        Either an appropriatley sized array, or a reference to one.
    dim
        The dimension being masked.
    allow_probabilities
        Whether to allow probabilities as `mask`
    """
    if mask is None:
        return mask
    desc = "mask/probabilities" if allow_probabilities else "mask"

    if isinstance(mask, str | AdRef):
        mask = _resolve_ref(mask)
        if not isinstance(data, AnnData):
            msg = f"Cannot use refererence for {desc} without providing anndata object as argument"
            raise ValueError(msg)
        try:
            mask_array = np.asarray(_get_vec_compat(data, mask, dim=dim))
        except KeyError:
            if isinstance(mask, AdRef):
                msg = (
                    f"Did not find `{mask}` in `adata`. "
                    f"Either add the {desc} first to `adata.{dim}`"
                    f"or consider using the {desc} argument with an array."
                )
            else:
                msg = f"Did not find `adata.{dim}[{mask!r}]`. "
            raise ValueError(msg) from None
    else:
        if len(mask) != data.shape[0 if dim == "obs" else 1]:
            msg = f"The shape of the {desc} do not match the data."
            raise ValueError(msg)
        mask_array = mask

    is_bool = pd.api.types.is_bool_dtype(mask_array.dtype)
    if not allow_probabilities and not is_bool:
        msg = "Mask array must be boolean."
        raise ValueError(msg)
    elif allow_probabilities and not (
        is_bool or pd.api.types.is_float_dtype(mask_array.dtype)
    ):
        msg = f"{desc} array must be boolean or floating point."
        raise ValueError(msg)

    return mask_array


@overload
def _resolve_ref[R: AdRef](ref: R | str) -> R | str: ...


@overload
def _resolve_ref[R: AdRef](ref: Collection[R | str]) -> list[R] | list[str]: ...


@overload
def _resolve_ref(
    ref: str, *, vec: Literal[False]
) -> LayerAcc | MultiAcc | GraphAcc | str: ...


def _resolve_ref(
    ref: AdRef | str | Collection[AdRef | str],
    *,
    vec: bool = True,
) -> AdRef | LayerAcc | MultiAcc | GraphAcc | str | list[AdRef] | list[str]:
    """Resolve a plain string `ref` into an `AdRef` under the `v2` preset, else leave it as-is.

    With `vec=False`, resolve into a whole array’s accessor (e.g. `A.obsm['x']`) instead,
    for parameters whose strings used to be keys, e.g. of `.obsm`.
    """
    if not vec:
        from scanpy import settings

        if settings.preset is not Preset.ScanpyV2Preview:
            return ref
        from anndata.acc import A

        return A.resolve(ref, vec=False)
    if not isinstance(ref, AdRef | str):
        refs = [_resolve_ref(r) for r in ref]
        if sum(isinstance(r, str) for r in refs) not in {0, len(refs)}:
            msg = "All elements of `refs` must be either `AdRef` or strings, not a mix, if `preset` is not `ScanpyV2Preview`."
            raise TypeError(msg)
        return refs

    from scanpy import settings

    if isinstance(ref, AdRef) or settings.preset is not Preset.ScanpyV2Preview:
        return ref
    from anndata.acc import A

    return A.resolve(ref, vec=True)


def _ref_arg[R: AdRef | str | Collection[AdRef | str] | None](
    ref: R, legacy: R, *, name: str, legacy_name: str, compat: bool = False
) -> R:
    """Merge a reference parameter `name` with the deprecated one it replaces.

    The deprecated parameter’s strings keep naming `.obs` columns (except with the `v2` preset).
    The new one’s are always resolved, unless it was renamed from a positional parameter
    that took column names (`compat`), so positional calls keep working.
    """
    if legacy is not None:
        if ref is not None:
            msg = f"Pass either `{name}` or `{legacy_name}`, not both."
            raise TypeError(msg)
        return _resolve_ref(legacy)
    if ref is None or compat:
        return ref if ref is None else _resolve_ref(ref)
    if isinstance(ref, str) or (
        isinstance(ref, Collection) and any(isinstance(r, str) for r in ref)
    ):
        from anndata.acc import A

        if isinstance(ref, str):
            return A.resolve(ref, vec=True)
        return [A.resolve(r, vec=True) if isinstance(r, str) else r for r in ref]
    return ref


def _ref_dim(
    ref: AdRef | str, *, dim: Literal["obs", "var"] | None
) -> Literal["obs", "var"]:
    """Derive the dim a single `ref` refers to, validating it against `dim` if given."""
    ref = _resolve_ref(ref)
    if isinstance(ref, str):
        return dim or "obs"
    ref_dim = next(iter(ref.dims))
    if dim is not None and dim != ref_dim:
        msg = f"Dimension of `{ref}` ({ref_dim}) does not match `{dim}`."
        raise ValueError(msg)
    return ref_dim


def _refs_dim(
    refs: Collection[AdRef | str], *, dim: Literal["obs", "var"] | None = None
) -> Literal["obs", "var"]:
    """Derive the single dim a collection of `refs` all refer to."""
    dims = {_ref_dim(r, dim=dim) for r in refs}
    if len(dims) != 1:
        msg = (
            f"All refs must refer to the same single axis (`obs` or `var`), got {dims}"
        )
        raise ValueError(msg)
    return next(iter(dims))


def _fetch_vec(adata: AnnData, ref: AdRef | str, *, dim: Literal["obs", "var"]) -> Any:
    """Retrieve the 1D array a single, already dim-resolved `ref` points to."""
    ref = _resolve_ref(ref)
    if isinstance(ref, str):
        return getattr(adata, dim)[ref]
    return adata[ref]


@overload
def _get_vec_compat(
    adata: AnnData,
    ref: Collection[AdRef] | Collection[str],
    *,
    dim: Literal["obs", "var"] | None = None,
) -> list[Any]: ...


@overload
def _get_vec_compat(
    adata: AnnData, ref: AdRef | str, *, dim: Literal["obs", "var"] | None = None
) -> Any: ...


def _get_vec_compat(
    adata: AnnData,
    ref: AdRef | str | Collection[AdRef] | Collection[str],
    *,
    dim: Literal["obs", "var"] | None = None,
) -> Any:
    """Get the 1D array(s) one or more `ref`erences point to.

    Treats strings as `obs` columns instead of `anndata.acc` specs when the preset is v1.
    """
    if _collection_of(ref, (AdRef, str)):
        dim = _refs_dim(ref, dim=dim)
        return [_fetch_vec(adata, r, dim=dim) for r in ref]

    if isinstance(ref, Collection) and not isinstance(ref, str):
        msg = f"Expected a single ref or collection of refs, got {ref!r}"
        raise TypeError(msg)

    dim = _ref_dim(ref, dim=dim)
    return _fetch_vec(adata, ref, dim=dim)


def _obs_df(
    adata: AnnData, refs: AdRef | str | Collection[AdRef | str]
) -> pd.DataFrame:
    """Gather `.obs`-aligned vectors (column names or refs) into a frame.

    Its columns are labelled for display (e.g. `A.obs['batch']`), so only use them as labels.
    """
    refs = [refs] if isinstance(refs, AdRef | str) else list(refs)
    vals = _get_vec_compat(adata, refs, dim="obs")
    return pd.DataFrame(
        {str(_resolve_ref(r)): v for r, v in zip(refs, vals, strict=True)},
        index=adata.obs_names,
    )


@overload
def _get_vec(
    adata: AnnData,
    ref: Collection[AdRef] | Collection[str],
    *,
    dim: Literal["obs", "var"] | None = None,
) -> list[Any]: ...


@overload
def _get_vec(
    adata: AnnData, ref: AdRef | str, *, dim: Literal["obs", "var"] | None = None
) -> Any: ...


def _get_vec(
    adata: AnnData,
    ref: AdRef | str | Collection[AdRef] | Collection[str],
    *,
    dim: Literal["obs", "var"] | None = None,
) -> Any:
    """Get the 1D array(s) one or more `ref`erences point to, using `anndata.acc`."""
    from anndata.acc import A, AdRef

    if _collection_of(ref, (AdRef, str)):
        refs = [A.resolve(r, vec=True) if isinstance(r, str) else r for r in ref]
        _refs_dim(refs, dim=dim)
        return [adata[r] for r in refs]

    if isinstance(ref, Collection) and not isinstance(ref, str):
        msg = f"Expected a single ref or collection of refs, got {ref!r}"
        raise TypeError(msg)

    if isinstance(ref, str):
        ref = A.resolve(ref, vec=True)
    _ref_dim(ref, dim=dim)
    return adata[ref]


@overload
def _resolve_obs[U: ArrAcc | str | None](
    use: U, /, *, rep: Literal[False] = False
) -> U: ...


@overload
def _resolve_obs(use: RefAcc | str, /, *, rep: Literal[True]) -> RepAcc: ...


def _resolve_obs(use: RefAcc | str | None, /, *, rep: bool = False) -> RefAcc | None:
    """Resolve a `use`/`out` string and require the result to be `obs`-aligned.

    Every function that can only handle `obs`-aligned input funnels through here,
    so `grep _resolve_obs` lists the ones that could be extended to work on `var`.

    With `rep=True`, additionally require a `rep`resentation (`LayerAcc`/`MultiAcc`).
    """
    if use is None:
        return use
    if isinstance(use, str):
        from anndata.acc import A

        use = A.resolve(use, vec=False)
    if isinstance(use, MultiAcc | GraphAcc) and use.dim != "obs":
        msg = f"Input must be aligned to `obs`, but {use!r} is aligned to `var`"
        raise ValueError(msg)
    if rep and not isinstance(use, LayerAcc | MultiAcc):
        msg = (
            "Representation must be a `LayerAcc` (e.g. `A.X`, `A.layers[...]`) or a "
            f"`MultiAcc` (e.g. `A.obsm[...]`), was {use!r}"
        )
        raise TypeError(msg)
    return use
