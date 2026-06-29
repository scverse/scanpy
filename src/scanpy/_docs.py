"""Shared docstrings for general parameters."""

from __future__ import annotations

__all__ = ["doc_rng"]

doc_rng = """\
rng
    Random number generation to control stochasticity.

    If a type:`SeedLike` value, it’s used to seed a new random number generator;
    If a :class:`numpy.random.Generator`, `rng`’s state will be directly advanced;
    If :data:`None`, a non-reproducible random number generator is used.
    See :func:`numpy.random.default_rng` for more details.

    The default value matches legacy scanpy behavior and will change to `None` in scanpy 2.0.\
"""
