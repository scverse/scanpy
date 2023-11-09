from __future__ import annotations

from types import FunctionType
from typing import TypeVar
from collections.abc import Callable


F = TypeVar("F", bound=FunctionType)


def doctest_needs(mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    def decorator(func: F) -> F:
        try:
            from ._pytest.marks import needs
        except ImportError:
            return func
        try:
            func._doctest_mark = needs[mod]
        except KeyError:
            raise KeyError(
                f"Unknown dependency {mod}. If it isn’t a typo, "
                "please add it to `needs` enum in `scanpy.testing._pytests.marks`."
            ) from None
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
