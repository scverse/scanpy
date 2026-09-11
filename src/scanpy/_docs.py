"""Shared docstrings for general parameters."""

from __future__ import annotations

from typing import TYPE_CHECKING

from scverse_misc import Deprecation

if TYPE_CHECKING:
    from typing import Final, Literal

__all__ = ["DEPR_COPY", "doc_mask", "doc_out", "doc_ref_compat", "doc_rng", "doc_use"]

DEPR_COPY: Final = Deprecation(
    "1.13.0", "Copy `adata` before calling this function instead."
)
"""For functions that gained an `out` parameter, making `copy` redundant."""

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
    legacy_docs = "".join(
        f"{name}\n"
        f"    A key of :attr:`~anndata.AnnData.{attr}`,\n"
        f"    i.e. `{name}='k'` is `use=A.{attr}['k']`.\n"
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


def doc_out(default: str) -> str:
    """Docs for an `out` parameter."""
    return f"""\
out
    Where to write the result, e.g. `A.layers['scaled']`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d, e.g. `'layers.scaled'`.
    If :data:`None`, the result is returned instead of written.
    If not given, it is written to {default}.
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
