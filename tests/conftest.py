from __future__ import annotations

import shutil
import sys
from contextlib import ExitStack
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING, TypedDict, cast

import anndata as ad
import pytest
from packaging.version import Version

# just import for the IMPORTED check
from scanpy._compat import pkg_version

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
        sc.settings.root_logger,
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
    """Allow use of scanpy’s logger with caplog."""
    import scanpy as sc

    sc.settings.root_logger.addHandler(caplog.handler)
    yield
    sc.settings.root_logger.removeHandler(caplog.handler)


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
def check_same_image(cache: pytest.Cache):
    from urllib.parse import quote

    from matplotlib.testing.compare import compare_images

    def check_same_image(
        expected: Path | os.PathLike,
        actual: Path | os.PathLike,
        *,
        tol: int,
        root: Path,
        save: bool = True,
    ) -> None:
        __tracebackhide__ = True

        result = cast(
            "CompareResult | None",
            compare_images(str(expected), str(actual), tol=tol, in_decorator=True),
        )
        if result is None:
            return

        if save:
            d = cache.mkdir("debug")
            for image in ("expected", "actual", "diff"):
                src = Path(result[image])
                dst = d / src.relative_to(root)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

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

    def save_and_compare(root: Path, path_str: Path | os.PathLike, *, tol: int):
        __tracebackhide__ = True

        base_pth = root / path_str

        if not base_pth.is_dir():
            base_pth.mkdir()
        expected_pth = base_pth / "expected.png"
        actual_pth = base_pth / "actual.png"
        plt.savefig(actual_pth, dpi=40)
        plt.close()
        if not expected_pth.is_file():
            msg = f"No expected output found at {expected_pth}."
            raise OSError(msg)
        check_same_image(expected_pth, actual_pth, tol=tol, root=root)

    return save_and_compare


@pytest.fixture
def plt():
    from matplotlib import pyplot as plt

    return plt


@pytest.fixture
def exit_stack() -> Generator[ExitStack]:
    with ExitStack() as stack:
        yield stack


@pytest.fixture(autouse=True)
def anndata_settings():
    if pkg_version("anndata") >= Version("0.12"):
        ad.settings.auto_shard_zarr_v3 = True
        ad.settings.zarr_write_format = 3
