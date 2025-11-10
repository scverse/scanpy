from __future__ import annotations

from collections.abc import Callable


def doctest_needs[F: Callable](mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    def decorator(func: F) -> F:
        func._doctest_needs = mod
        return func

    return decorator


def doctest_skip[F: Callable](reason: str) -> Callable[[F], F]:
    """Mark function so doctest is skipped."""
    if not reason:
        msg = "reason must not be empty"
        raise ValueError(msg)

    def decorator(func: F) -> F:
        func._doctest_skip_reason = reason
        return func

    return decorator


def doctest_internet[F: Callable](func: F) -> F:
    """Mark function so doctest gets the internet mark."""
    func._doctest_internet = True
    return func
