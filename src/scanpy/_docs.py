"""Shared docstrings for general parameters."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scverse_misc import Deprecation

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Final, Literal

__all__ = [
    "DEPR_COPY",
    "DEPR_KEY_ADDED",
    "DEPR_RAW",
    "doc_mask",
    "doc_out",
    "doc_out_key",
    "doc_out_multi",
    "doc_ref_compat",
    "doc_rng",
    "doc_use",
]

DEPR_COPY: Final = Deprecation(
    "1.13.0", "Copy `adata` before calling this function instead."
)
"""For functions that gained an `out` parameter, making `copy` redundant."""

DEPR_KEY_ADDED: Final = Deprecation("1.13.0", "Use `out` instead.")
"""For `key_added` parameters, replaced by `out` (e.g. `out=sc.Prefix(key_added)`)."""

DEPR_RAW: Final = Deprecation(
    "1.13.0",
    ":attr:`~anndata.AnnData.raw` is deprecated. Keep the data in a layer instead.",
)
"""For `use_raw` parameters, which have no accessor equivalent."""

doc_ref_compat = (
    "If :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`, "
    ":class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d "
    "to :class:`~anndata.acc.AdRef`\\ s (e.g. `'obs.is_control'`), "
    "otherwise interpreted as column names."
)


def doc_mask(desc: str, *, dim: Literal["obs", "var"], extra: str = "") -> str:
    """Docs for a `mask` parameter and the deprecated `mask_{dim}` it replaces."""
    return f"""\
mask
    {desc}
    Given by a boolean array or a reference to one, e.g. `A.{dim}['selected']`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d, e.g. `'{dim}.selected'`.
{f"    {extra}\n" if extra else ""}\
mask_{dim}
    A boolean array, or a column of :attr:`~anndata.AnnData.{dim}` named by a :class:`str`,
    i.e. `mask_{dim}='selected'` is `mask=A.{dim}['selected']`.
"""


def doc_use(desc: str, *, legacy: tuple[str, ...] = ("layer", "obsm")) -> str:
    """Docs for a `use` parameter and the deprecated parameters it replaces."""
    # no trailing newline: the caller’s `{use}` is followed by its own one
    legacy_docs = "\n".join(
        f"{name}\n"
        f"    A key of :attr:`~anndata.AnnData.{attr}`,\n"
        f"    i.e. `{name}='k'` is `use=A.{attr}['k']`."
        for name in legacy
        if (attr := "layers" if name == "layer" else name)
    )
    return f"""\
use
    {desc}
    Given by an accessor, e.g. `A.X`, `A.layers['counts']` or `A.obsm['pca']`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d, e.g. `'layers.counts'`.
{legacy_docs}\
"""


def doc_out(base: str, *, default: str = "the place `use` reads from") -> str:
    """Docs for an `out` parameter of a function with one output, e.g. a matrix.

    `base` is the output’s name, which `sc.Prefix('p')` prefixes:
    a matrix then goes to `A.layers['p_{base}']` (or `.obsm` if `use` points there),
    other outputs to the same kind of place as `default`.
    """
    return f"""\
out
    Where to write the result, e.g. `A.layers['scaled']`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d, e.g. `'layers.scaled'`.
    `sc.Prefix('p')` writes it to the same kind of place as the default, named `'p_{base}'`.
    If :data:`True`, the result is returned instead of written.
    If :data:`False`, it is written to {default}.\
"""


doc_rng = """\
rng
    Random number generation to control stochasticity.

    If a type:`SeedLike` value, it’s used to seed a new random number generator;
    If a :class:`numpy.random.Generator`, `rng`’s state will be directly advanced;
    If :data:`None`, a non-reproducible random number generator is used.
    See :func:`numpy.random.default_rng` for more details.

    The default value matches legacy scanpy behavior and will change to `None` in scanpy 2.0.\
"""


def doc_out_key(
    what: str,
    *,
    default: str,
    prefixed: str,
    key_added: str | None = None,
    example: str | None = None,
) -> str:
    """Docs for the `out` parameter of a function with one output, and the `key_added` it replaces.

    `prefixed` is where `out=sc.Prefix('p')` writes, e.g. "`.obsm['p_umap']` and `.uns['p_umap']`",
    `key_added` says what `key_added='k'` did and how to say it with `out`.
    Without `example`, the output only goes to `.uns`, so `out` can’t be an accessor.
    """
    acc = (
        ""
        if example is None
        else (
            f"    Given by an accessor, e.g. `{example}`.\n"
            "    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d.\n"
            "    Accompanying :attr:`~anndata.AnnData.uns` entries are stored under the same key.\n"
        )
    )
    return f"""\
out
    Where to write {what}.
{acc}\
    `sc.Prefix('p')` writes it to {prefixed}.
    If :data:`True`, {what} is returned instead of written.
    If :data:`False`, it is written to {default}.\
{_doc_key_added(key_added)}"""


def doc_out_multi(
    outputs: Mapping[str, str],
    *,
    example: str,
    returns: str,
    key_added: str | None = None,
    uns: str | None = None,
) -> str:
    """Docs for the `out` parameter of a function with several outputs, and the `key_added` it replaces.

    `outputs` maps each output’s name to a description including its default location
    and where `out=sc.Prefix('p')` writes it.
    `uns` is the name of the accompanying `.uns` entry, if any, which a `Prefix` prefixes.
    `key_added` says what `key_added='k'` did and how to say it with `out`.
    """
    items = "\n".join(f"    `{name}`\n        {desc}" for name, desc in outputs.items())
    uns_doc = (
        "" if uns is None else f" (and :attr:`~anndata.AnnData.uns`\\ `['p_{uns}']`)"
    )
    return f"""\
out
    Where to write the outputs.
    A mapping from output names to accessors writes only the outputs it names,
    e.g. `{example}`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d.
    `sc.Prefix('p')` writes each output under its name below, prefixed by `p_`{uns_doc}.
    If :data:`True`, they are returned as a {returns} instead of written.
    If :data:`False`, they are written to their default locations.
    The outputs are:

{items}\
{_doc_key_added(key_added)}"""


def _doc_key_added(desc: str | None) -> str:
    return "" if desc is None else f"\nkey_added\n    {desc}"
