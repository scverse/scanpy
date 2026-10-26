"""Resolving and writing outputs: the `out` parameter, and deprecated `key_added`/`copy`."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from functools import wraps
from inspect import signature
from typing import TYPE_CHECKING, cast

from anndata import AnnData
from scverse_misc import deprecated_arg

from .._docs import DEPR_COPY
from .._settings import Default
from ._types import AdRef, GraphAcc, LayerAcc, MultiAcc, RefAcc

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal, Unpack

    from ._types import ArrAcc, VecRef, _SetRep


def _set_arr(
    adata: AnnData,
    val: Any,
    out: ArrAcc | VecRef | None = None,
    **choices: Unpack[_SetRep],
) -> None:
    """Write `val`, via an `anndata.acc` accessor or old-style choices.

    The inverse of `_get_arr`; see there for the accessor/`choices` split.
    """
    if out is not None:
        if any(v not in (None, False) for v in choices.values()):
            msg = "`out` cannot be combined with `layer`/`obsm`/`obsp`/`varm`/`varp`"
            raise TypeError(msg)
        choices = _dest(out)

    match [(k, v) for k, v in choices.items() if v not in {None, False}]:
        case []:
            adata.X = val
        case [("layer", layer)]:
            adata.layers[layer] = val
        case [(("obsm" | "obsp" | "varm" | "varp" | "obs" | "var") as k, v)]:
            getattr(adata, k)[v] = val
        case picked:
            valid = [f"`{k}`" for k, _ in picked]
            valid[-1] = f"or {valid[-1]}"
            msg = f"Only one of {', '.join(valid)} can be specified."
            raise ValueError(msg)


def _write_out[V](
    adata: AnnData,
    val: V,
    *,
    out: ArrAcc | VecRef | bool,
    use: ArrAcc | VecRef | None = None,
    **choices: Unpack[_SetRep],
) -> V | None:
    """Honor an `out` argument: write `val` and return `None`, or return `val` unwritten.

    `out` is either an accessor to write to,
    `True` (return `val` instead of writing it),
    or `False` (the place `use`/`choices` read from).
    """
    if isinstance(out, AdRef | RefAcc | str):
        _set_arr(adata, val, out)
        return None
    if out is True:
        return val
    _set_arr(adata, val, use, **choices)
    return None


def _use_dest(use: ArrAcc | VecRef | None, **choices: Unpack[_SetRep]) -> _SetRep:
    """Where `use`, or the legacy `layer`/`obsm`/… `choices`, read from, as a destination."""
    if use is not None:
        return _dest(use)
    return cast("_SetRep", {k: v for k, v in choices.items() if v not in {None, False}})


def _dest(acc: ArrAcc | VecRef) -> _SetRep:
    """Translate an `anndata.acc` accessor into `_set_arr`’s old-style choices.

    E.g. `A.obsm['x']` becomes `{'obsm': 'x'}`, and `A.X` becomes `{}`.
    """
    from anndata.acc import A, MetaAcc

    if isinstance(acc, str):
        acc = A.resolve(acc)
    match acc:
        case LayerAcc(k=None):
            return {}
        case LayerAcc(k=key):
            return dict(layer=key)
        case MultiAcc(dim=dim, k=key):
            return cast("_SetRep", {f"{dim}m": key})
        case GraphAcc(dim=dim, k=key):
            return cast("_SetRep", {f"{dim}p": key})
        # assign by column. Assigning to an 1D array slice is done elsewhere
        case AdRef(acc=MetaAcc(dim=dim), idx=str(col)):
            return cast("_SetRep", {dim: col})
        case _:
            msg = f"Cannot write to {acc!r}"
            raise TypeError(msg)


type OutKind = Literal["layer", "obsm", "varm", "obsp", "varp", "obs", "var", "uns"]
"""Kind of place an output can be written to; `'layer'` includes `.X`.

There are no `.uns` accessors yet (https://github.com/scverse/anndata/issues/2657),
so outputs of kind `'uns'` can only be named by `out=False` or the function’s `key_added`.
"""


def _dest_kind(dest: _SetRep) -> OutKind:
    return next(iter(dest), "layer")


def _dest_key(dest: _SetRep | None) -> str | None:
    return None if dest is None else next(iter(dest.values()), None)


def _legacy_out[O](out: O, name: str, value: O | None, *, default: O = False) -> O:
    """Merge a deprecated parameter into `out`, e.g. `inplace=False` into `out=True`.

    `value` is what the deprecated parameter `name` means as `out`, or `None` if it wasn’t passed.
    `default` is `out`’s default, i.e. what it has to be when `name` is passed.
    """
    if value is None:
        return out
    _check_legacy(out, name, default=default)
    return value


def _check_legacy(out: object, name: str, *, default: object = False) -> None:
    """Check that a deprecated parameter `name` isn’t passed together with `out`."""
    if out is not default:
        msg = f"Pass either `out` or `{name}`, not both."
        raise TypeError(msg)


def deprecated_copy[**P, R](func: Callable[P, R]) -> Callable[P, R | AnnData]:
    """Deprecate `func`’s `copy` parameter, and implement it for `AnnData` input.

    With `copy=True`, `func` runs with `copy=False` on a copy of its first argument,
    which is then returned (unless `func` returns something else).
    So functions don’t need to handle `copy` themselves,
    except for non-`AnnData` input (e.g. arrays), which is passed through unchanged.
    """
    sig = signature(func)
    first = next(iter(sig.parameters))

    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R | AnnData:
        bound = sig.bind(*args, **kwargs)
        if not bound.arguments.get("copy") or not isinstance(
            data := bound.arguments[first], AnnData
        ):
            return func(*args, **kwargs)
        if bound.arguments.get("out") is True:
            msg = "`copy=True` cannot be used with `out=True`."
            raise ValueError(msg)
        bound.arguments[first] = data = data.copy()
        bound.arguments["copy"] = False
        res = func(*bound.args, **bound.kwargs)
        return data if res is None else res

    return deprecated_arg("copy", DEPR_COPY)(wrapper)


def _out_key[D: str | Default | None](
    out: ArrAcc | VecRef | bool,  # noqa: FBT001
    key_added: str | None,
    *,
    kind: OutKind,
    default: D = None,
) -> str | D | Literal[True]:
    """Translate the `out` of a function with one output into a key.

    The output goes to a place of `kind`.
    Returns `True` if the output should be returned instead of written,
    `default` if `out=False`, the deprecated `key_added` if given,
    and otherwise the key of the accessor.
    Keys from `out` are literal, i.e. shouldn’t be processed further like `key_added` is.
    `key_added` is ignored for `out=True` (e.g. from `inplace=False`), as nothing is written.
    """
    if key_added is not None and out is not True:
        _check_legacy(out, "key_added")
        return key_added
    match out:
        case True:
            return True
        case False:
            return default
    if kind == "uns":
        msg = f"`.uns` accessors aren’t supported yet, `out` must be a `bool`, was {out!r}"
        raise TypeError(msg)
    dest = _dest(out)
    if _dest_kind(dest) != kind or (key := _dest_key(dest)) is None:
        msg = f"`out` must point to `.{kind}`, was {out!r}"
        raise TypeError(msg)
    return key


def _out_dests[D: tuple, K: tuple](
    out: ArrAcc | VecRef | Mapping[str, ArrAcc | VecRef] | bool,  # noqa: FBT001
    key_added: str | None,
    *,
    dests: Callable[[str | Default | None], D],
    kinds: K,
    uns: bool = False,
    default: str | Default | None = None,
    primary: bool = True,
) -> tuple[str | Default | None, D] | Literal[True]:
    """Translate the `out` of a function with several outputs into one destination each.

    The outputs are described by instances of the function’s output `NamedTuple`,
    whose first field is the primary output:
    `kinds` holds the `OutKind` values each output may be written to,
    and `dests(key_added)` where the deprecated `key_added` put each output
    (`_SetRep` dicts, or `None` to not write one).
    `uns` says if there’s an accompanying `.uns` entry, named after the primary output.

    Returns `True` if the outputs should be returned instead of written.
    Otherwise returns a key for naming `.uns` entries and the destinations:
    `default` and `dests(default)` for `out=False`,
    `key_added` and `dests(key_added)` for the deprecated `key_added`,
    and the primary output’s key for an accessor or mapping
    (or `default` if that is where `out=False` would write it).

    A single accessor (if `primary`) names the primary output `k`, the others follow it:
    an output in a kind of place no other output goes to is named `k` too,
    else `{k}_{name}`, and outputs `out=False` doesn’t write stay unwritten.
    If the primary output is `.X`, the others go to their defaults.
    A mapping writes only the outputs it names, and must name the primary output if there is a `uns` entry.
    `key_added` is ignored for `out=True` (e.g. from `inplace=False`), as nothing is written.
    """
    if key_added is not None and out is not True:
        _check_legacy(out, "key_added")
        return key_added, dests(key_added)
    names: tuple[str, ...] = kinds._fields  # type: ignore[attr-defined]
    defaults = dests(default)
    match out:
        case True:
            return True
        case False:
            return default, defaults
        case AdRef() | RefAcc() | str() if primary:
            first = _check_kind(kinds, names[0], out)
            if first == defaults[0]:  # naming the default place is `out=False`
                return default, defaults
            result = _follow_primary(first, defaults)
        case Mapping():
            result = _resolve_mapping(out, kinds, names)
        case _:
            what = (
                "a `bool`, an accessor for the primary output, or"
                if primary
                else "a `bool` or"
            )
            msg = f"`out` must be {what} a mapping with keys {names}, was {out!r}"
            raise TypeError(msg)
    # falling back to `default` would overwrite the `.uns` entry of outputs this call doesn’t write
    if uns and result[0] is None:
        msg = (
            f"`out` must include {names[0]!r}, after which the `.uns` entry is named, "
            f"until there are `.uns` accessors (https://github.com/scverse/anndata/issues/2657), "
            f"was {out!r}"
        )
        raise TypeError(msg)
    if result[0] == defaults[0] or (key := _dest_key(result[0])) is None:
        return default, result
    return key, result


def _check_kind(kinds: tuple, name: str, acc: ArrAcc | VecRef) -> _SetRep:
    dest = _dest(acc)
    if _dest_kind(dest) not in (valid := getattr(kinds, name)):
        msg = f"Output {name!r} must be written to one of {sorted(valid)}, was {acc!r}"
        raise TypeError(msg)
    return dest


def _resolve_mapping[D: tuple](
    out: Mapping[str, ArrAcc | VecRef], kinds: tuple, names: tuple[str, ...]
) -> D:
    if unknown := out.keys() - set(names):
        msg = f"Unknown output(s) {sorted(unknown)}, valid are {names}"
        raise TypeError(msg)
    resolved = {name: _check_kind(kinds, name, acc) for name, acc in out.items()}
    return type(kinds)(*(resolved.get(name) for name in names))


def _follow_primary[D: tuple](first: _SetRep, defaults: D) -> D:
    """Name the outputs after the primary one, `first`, see `_out_dests`."""
    if (key := _dest_key(first)) is None:  # `.X`
        return type(defaults)(first, *defaults[1:])
    kind_counts = Counter(
        _dest_kind(d) for d in (first, *defaults[1:]) if d is not None
    )
    return type(defaults)(
        first,
        *(
            None
            if d is None
            else {
                (kind := _dest_kind(d)): key
                if kind_counts[kind] == 1
                else f"{key}_{name}"
            }
            for name, d in zip(defaults._fields[1:], defaults[1:], strict=True)
        ),
    )


def _write_dests(
    adata: AnnData, vals: tuple, dests: tuple[_SetRep | None, ...]
) -> None:
    """Write each of `vals` to its destination from `_out_dests`, skipping `None` ones."""
    for val, dest in zip(vals, dests, strict=True):
        if dest is not None:
            _set_arr(adata, val, **dest)
