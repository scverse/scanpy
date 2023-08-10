from __future__ import annotations

from types import FunctionType
from typing import TYPE_CHECKING, TypeVar
from collections.abc import Callable


if TYPE_CHECKING:
    import pytest


F = TypeVar('F', bound=FunctionType)


def doctest_mark(mark: pytest.MarkDecorator) -> Callable[[F], F]:
    def decorator(func: F) -> F:
        func._doctest_mark = mark
        return func

    return decorator
