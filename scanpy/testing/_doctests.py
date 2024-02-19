from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

F = TypeVar("F", bound=Callable)


def doctest_needs(mod: str) -> Callable[[F], F]:
    """Mark function with doctest dependency."""

    try:
        from ._pytest.marks import needs
    except ImportError:
        mark = None
    else:
        try:
            mark = needs[mod]
        except KeyError:
            raise KeyError(
                f"Unknown dependency {mod}. If it isn’t a typo, "
                "please add it to `needs` enum in `scanpy.testing._pytests.marks`."
            ) from None

    def decorator(func: F) -> F:
        if mark is not None:
            func._doctest_needs = mark
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
