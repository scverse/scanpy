"""A private pytest plugin"""
from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

import pytest

from ..._utils import _import_name
from .fixtures import *  # noqa: F403

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable

    from .marks import needs


# Defining it here because it’s autouse.
@pytest.fixture(autouse=True)
def _global_test_context(request: pytest.FixtureRequest) -> Generator[None, None, None]:
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

    if isinstance(request.node, pytest.DoctestItem):
        _modify_doctests(request)

    yield

    plt.close("all")


@pytest.fixture(autouse=True, scope="session")
def limit_multithreading():
    """Limit number of threads used per worker when using pytest-xdist.

    Prevents oversubscription of the CPU when multiple tests with parallel code are
    running at once.
    """
    if "PYTEST_XDIST_WORKER_COUNT" in os.environ:
        import threadpoolctl

        n_workers = int(os.environ["PYTEST_XDIST_WORKER_COUNT"])
        max_threads = os.cpu_count() // n_workers

        with threadpoolctl.threadpool_limits(limits=max_threads):
            yield
    else:
        yield


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
    assert isinstance(request.node, pytest.DoctestItem)

    request.getfixturevalue("_doctest_env")

    func = _import_name(request.node.name)
    needs_marker: needs | None
    if needs_marker := getattr(func, "_doctest_needs", None):
        assert needs_marker.mark.name == "skipif"
        if needs_marker.mark.args[0]:
            pytest.skip(reason=needs_marker.mark.kwargs["reason"])
    skip_reason: str | None
    if skip_reason := getattr(func, "_doctest_skip_reason", None):
        pytest.skip(reason=skip_reason)
