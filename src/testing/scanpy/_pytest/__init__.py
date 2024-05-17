"""A private pytest plugin"""

from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

import pytest

from .fixtures import *  # noqa: F403
from .marks import needs

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable


# Defining it here because it’s autouse.
@pytest.fixture(autouse=True)
def _global_test_context(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
) -> Generator[None, None, None]:
    """Switch to agg backend, reset settings, and close all figures at teardown."""
    # make sure seaborn is imported and did its thing
    import seaborn as sns  # noqa: F401
    from matplotlib import pyplot as plt
    from matplotlib.testing import setup

    import scanpy as sc

    setup()
    sc.settings.logfile = sys.stderr
    sc.settings.verbosity = "hint"
    sc.settings.autoshow = True
    sc.settings.datasetdir = tmp_path_factory.mktemp("scanpy_data")
    sc.settings.writedir = tmp_path_factory.mktemp("scanpy_write")

    if isinstance(request.node, pytest.DoctestItem):
        _modify_doctests(request)

    yield

    plt.close("all")


@pytest.fixture(autouse=True, scope="session")
def max_threads() -> Generator[int, None, None]:
    """Limit number of threads used per worker when using pytest-xdist.

    Prevents oversubscription of the CPU when multiple tests with parallel code are
    running at once.
    """
    if (n_workers := os.environ.get("PYTEST_XDIST_WORKER_COUNT")) is not None:
        import threadpoolctl

        n_cpus = os.cpu_count() or 1
        n_workers = int(n_workers)
        max_threads = max(n_cpus // n_workers, 1)

        with threadpoolctl.threadpool_limits(limits=max_threads):
            yield max_threads
    else:
        yield 0


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help=(
            "Run tests that retrieve stuff from the internet. "
            "This increases test time."
        ),
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: Iterable[pytest.Item]
) -> None:
    import pytest

    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--run-internet` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)


def _modify_doctests(request: pytest.FixtureRequest) -> None:
    from scanpy._utils import _import_name

    assert isinstance(request.node, pytest.DoctestItem)

    request.getfixturevalue("_doctest_env")

    func = _import_name(request.node.name)
    needs_mod: str | None
    if needs_mod := getattr(func, "_doctest_needs", None):
        needs_marker = needs[needs_mod]
        if needs_marker.mark.args[0]:
            pytest.skip(reason=needs_marker.mark.kwargs["reason"])
    skip_reason: str | None
    if skip_reason := getattr(func, "_doctest_skip_reason", None):
        pytest.skip(reason=skip_reason)
    if getattr(func, "_doctest_internet", False) and not request.config.getoption(
        "--internet-tests"
    ):
        pytest.skip(reason="need --internet-tests option to run")


def pytest_itemcollected(item: pytest.Item) -> None:
    # Dask AnnData tests require anndata > 0.10
    import anndata
    from packaging.version import Version

    requires_anndata_dask_support = (
        len([mark for mark in item.iter_markers(name="anndata_dask_support")]) > 0
    )

    if requires_anndata_dask_support and Version(anndata.__version__) < Version("0.10"):
        item.add_marker(
            pytest.mark.skip(reason="dask support requires anndata version > 0.10")
        )


assert (
    "scanpy" not in sys.modules
), "scanpy is already imported, this will mess up test coverage"
