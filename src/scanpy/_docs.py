"""Shared docstrings for general parameters."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

__all__ = ["doc_mask", "doc_ref_compat", "doc_rng"]

doc_ref_compat = (
    "If :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`, "
    ":class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d "
    "to :class:`~anndata.acc.AdRef`\\ s (e.g. `'obs.is_control'`), "
    "otherwise interpreted as column names."
)


def doc_mask(desc: str, *, dim: Literal["obs", "var"], extra: str = "") -> str:
    """Docs for a `mask` parameter and its deprecated `mask_{dim}` alias."""
    return f"""\
mask
    {desc}
    Given by a boolean array or a reference to one, e.g. `A.{dim}['selected']`.
    :class:`str`\\ s are :meth:`anndata.acc.AdAcc.resolve`\\ d, e.g. `'{dim}.selected'`.
{f"    {extra}\n" if extra else ""}\
mask_{dim}
    Deprecated alias of `mask`, where a :class:`str` always refers to
    a column of :attr:`~anndata.AnnData.{dim}`.
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
