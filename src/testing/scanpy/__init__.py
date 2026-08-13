"""Scanpy testing utilities."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._backends import validate_backend

__all__ = ["validate_backend"]


def __getattr__(name: str):
    if name == "validate_backend":
        from ._backends import validate_backend

        return validate_backend
    raise AttributeError(name)
