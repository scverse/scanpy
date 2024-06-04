from __future__ import annotations

import sys
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING, TypedDict, Union, cast

import pytest

# just import for the IMPORTED check
import scanpy as _sc  # noqa: F401

if TYPE_CHECKING:  # So editors understand that we’re using those fixtures
    import os
    from collections.abc import Generator

    from testing.scanpy._pytest.fixtures import *  # noqa: F403

# define this after importing scanpy but before running tests
IMPORTED = frozenset(sys.modules.keys())


@pytest.fixture(scope="session", autouse=True)
def _manage_log_handlers() -> Generator[None, None, None]:
    """Remove handlers from all loggers on session teardown.

    Fixes <https://github.com/scverse/scanpy/issues/1736>.
    See also <https://github.com/pytest-dev/pytest/issues/5502>.
    """
    import logging

    import scanpy as sc

    yield

    loggers = [
        sc.settings._root_logger,
        logging.getLogger(),
        *logging.Logger.manager.loggerDict.values(),
    ]
    for logger in loggers:
        if not isinstance(logger, logging.Logger):
            continue  # loggerDict can contain `logging.Placeholder`s
        for handler in logger.handlers[:]:
            if isinstance(handler, logging.StreamHandler):
                logger.removeHandler(handler)


@pytest.fixture(autouse=True)
def _caplog_adapter(caplog: pytest.LogCaptureFixture) -> Generator[None, None, None]:
    """Allow use of scanpy’s logger with caplog"""
    import scanpy as sc

    sc.settings._root_logger.addHandler(caplog.handler)
    yield
    sc.settings._root_logger.removeHandler(caplog.handler)


@pytest.fixture
def imported_modules():
    return IMPORTED


class CompareResult(TypedDict):
    rms: float
    expected: str
    actual: str
    diff: str
    tol: int


@pytest.fixture
def check_same_image(add_nunit_attachment):
    from urllib.parse import quote

    from matplotlib.testing.compare import compare_images

    def check_same_image(
        expected: Path | os.PathLike,
        actual: Path | os.PathLike,
        *,
        tol: int,
        basename: str = "",
    ) -> None:
        __tracebackhide__ = True

        def fmt_descr(descr):
            return f"{descr} ({basename})" if basename else descr

        result = cast(
            Union[CompareResult, None],
            compare_images(str(expected), str(actual), tol=tol, in_decorator=True),
        )
        if result is None:
            return

        add_nunit_attachment(result["expected"], fmt_descr("Expected"))
        add_nunit_attachment(result["actual"], fmt_descr("Result"))
        add_nunit_attachment(result["diff"], fmt_descr("Difference"))

        result_urls = {
            k: f"file://{quote(v)}" if isinstance(v, str) else v
            for k, v in result.items()
        }
        msg = dedent(
            """\
            Image files did not match.
            RMS Value:  {rms}
            Expected:   {expected}
            Actual:     {actual}
            Difference: {diff}
            Tolerance:  {tol}
            """
        ).format_map(result_urls)
        raise AssertionError(msg)

    return check_same_image


@pytest.fixture
def image_comparer(check_same_image):
    from matplotlib import pyplot as plt

    def save_and_compare(*path_parts: Path | os.PathLike, tol: int):
        __tracebackhide__ = True

        base_pth = Path(*path_parts)

        if not base_pth.is_dir():
            base_pth.mkdir()
        expected_pth = base_pth / "expected.png"
        actual_pth = base_pth / "actual.png"
        plt.savefig(actual_pth, dpi=40)
        plt.close()
        if not expected_pth.is_file():
            raise OSError(f"No expected output found at {expected_pth}.")
        check_same_image(expected_pth, actual_pth, tol=tol)

    return save_and_compare


@pytest.fixture
def plt():
    from matplotlib import pyplot as plt

    return plt
