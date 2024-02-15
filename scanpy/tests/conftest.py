from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

# just import for the IMPORTED check
import scanpy as _sc  # noqa: F401

if TYPE_CHECKING:  # So editors understand that we’re using those fixtures
    import os
    from collections.abc import Generator

    from scanpy.testing._pytest.fixtures import *  # noqa: F403

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


@pytest.fixture(autouse=True)
def _caplog_adapter(caplog: pytest.LogCaptureFixture):
    import scanpy as sc

    sc.settings._root_logger.addHandler(caplog.handler)
    yield
    sc.settings._root_logger.removeHandler(caplog.handler)


@pytest.fixture
def imported_modules():
    return IMPORTED


@pytest.fixture
def check_same_image(add_nunit_attachment):
    from matplotlib.testing.compare import compare_images, make_test_filename

    def _(pth1, pth2, *, tol: int, basename: str = ""):
        def fmt_descr(descr):
            if basename != "":
                return f"{descr} ({basename})"
            else:
                return descr

        pth1, pth2 = Path(pth1), Path(pth2)
        try:
            result = compare_images(str(pth1), str(pth2), tol=tol)
            assert result is None, result
        except Exception as e:
            diff_pth = make_test_filename(pth2, "failed-diff")
            add_nunit_attachment(str(pth1), fmt_descr("Expected"))
            add_nunit_attachment(str(pth2), fmt_descr("Result"))
            if Path(diff_pth).is_file():
                add_nunit_attachment(str(diff_pth), fmt_descr("Difference"))
            raise e

    return _


@pytest.fixture
def image_comparer(check_same_image):
    from matplotlib import pyplot as plt

    def save_and_compare(*path_parts: Path | os.PathLike, tol: int):
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


@pytest.fixture
def tmp_dataset_dir(tmp_path_factory):
    import scanpy

    new_dir = tmp_path_factory.mktemp("scanpy_data")
    old_dir = scanpy.settings.datasetdir
    scanpy.settings.datasetdir = new_dir  # Set up
    yield scanpy.settings.datasetdir
    scanpy.settings.datasetdir = old_dir  # Tear down


@pytest.fixture
def tmp_write_dir(tmp_path_factory):
    import scanpy

    new_dir = tmp_path_factory.mktemp("scanpy_write")
    old_dir = scanpy.settings.writedir
    scanpy.settings.writedir = new_dir  # Set up
    yield scanpy.settings.writedir
    scanpy.settings.writedir = old_dir  # Tear down
