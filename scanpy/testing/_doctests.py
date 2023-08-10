from __future__ import annotations

from types import FunctionType
from typing import TypeVar
from collections.abc import Callable


F = TypeVar('F', bound=FunctionType)


def doctest_needs(mod: str) -> Callable[[F], F]:
    def decorator(func: F) -> F:
        try:
            from ._pytest.marks import needs
        except ImportError:
            return func
        func._doctest_mark = needs(mod)
        return func

    return decorator
