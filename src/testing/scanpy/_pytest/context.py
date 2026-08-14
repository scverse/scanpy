from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Generator


# The following is necessary only for subtests: https://github.com/pytest-dev/pytest/issues/14101
@contextmanager
def xfail(
    condition: bool = True,  # noqa: FBT001, FBT002
    /,
    *,
    raises: type[Exception] = Exception,
    reason: str = "",
) -> Generator[None]:
    """If `condition` is true expect `raises` to be raised."""
    if not reason and raises is Exception:
        msg = "Either `reason` or `exc_cls` must be provided"
        raise RuntimeError(msg)

    if not condition:
        yield
        return

    try:
        yield
    except raises:
        pytest.xfail(reason=reason)
    else:
        pytest.fail(reason=f"[XPASS] {reason}")
