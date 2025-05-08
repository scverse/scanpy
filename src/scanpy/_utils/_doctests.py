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
        msg = "reason must not be empty"
        raise ValueError(msg)

    def decorator(func: F) -> F:
        func._doctest_skip_reason = reason
        return func

    return decorator


def doctest_internet(func: F) -> F:
    """Mark function so doctest gets the internet mark."""
    func._doctest_internet = True
    return func
