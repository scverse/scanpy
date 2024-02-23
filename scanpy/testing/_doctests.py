from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

F = TypeVar("F", bound=Callable)


def doctest_needs(mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    def decorator(func: F) -> F:
        func._doctest_needs = mod
        return func

    return decorator


def doctest_skip(reason: str) -> Callable[[F], F]:
    """Mark function so doctest is skipped."""
    if not reason:
        raise ValueError("reason must not be empty")

    def decorator(func: F) -> F:
        func._doctest_skip_reason = reason
        return func

    return decorator
