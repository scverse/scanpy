"""A private pytest plugin."""

from __future__ import annotations

import os
import sys
from types import MappingProxyType
from typing import TYPE_CHECKING

import pytest

from .fixtures import *  # noqa: F403
from .marks import needs

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Mapping

_original_settings: Mapping[str, object] | None = None


# Defining it here because it’s autouse.
@pytest.fixture(autouse=True)
def original_settings(
    request: pytest.FixtureRequest,
    cache: pytest.Cache,
    tmp_path_factory: pytest.TempPathFactory,
) -> Generator[Mapping[str, object], None, None]:
    """Switch to agg backend, reset settings, and close all figures at teardown."""
    # make sure seaborn is imported and did its thing
    import seaborn as sns  # noqa: F401
    from matplotlib import pyplot as plt
    from matplotlib.testing import setup

    import scanpy as sc

    global _original_settings  # noqa: PLW0603
    if _original_settings is None:
        _original_settings = MappingProxyType(sc.settings.__dict__.copy())

    setup()
    sc.settings.logfile = sys.stderr
    sc.settings.verbosity = "hint"
    sc.settings.autoshow = True
    # create directory for debug data
    cache.mkdir("debug")
    # reuse data files between test runs (unless overwritten in the test)
    sc.settings.datasetdir = cache.mkdir("scanpy-data")
    # create new writedir for each test run
    sc.settings.writedir = tmp_path_factory.mktemp("scanpy_write")

    if isinstance(request.node, pytest.DoctestItem):
        _modify_doctests(request)

    yield _original_settings

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
            "Run tests that retrieve stuff from the internet. This increases test time."
        ),
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: Iterable[pytest.Item]
) -> None:
    import pytest

    skipif_not_run_internet = pytest.mark.skipif(
        not config.getoption("--internet-tests"),
        reason="need --internet-tests option to run",
    )
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--run-internet` passed
        if "internet" in item.keywords:
            item.add_marker(skipif_not_run_internet)
            item.add_marker(pytest.mark.flaky(reruns=5, reruns_delay=2))


def _modify_doctests(request: pytest.FixtureRequest) -> None:
    from scanpy._utils import _import_name

    assert isinstance(request.node, pytest.DoctestItem)

    request.getfixturevalue("_doctest_env")

    func = _import_name(request.node.name)
    needs_mod: str | None
    skip_reason: str | None
    if (
        (needs_mod := getattr(func, "_doctest_needs", None))
        and (skip_reason := needs[needs_mod].skip_reason)
    ) or (skip_reason := getattr(func, "_doctest_skip_reason", None)):
        pytest.skip(reason=skip_reason)
    if getattr(func, "_doctest_internet", False):
        if not request.config.getoption("--internet-tests"):
            pytest.skip(reason="need --internet-tests option to run")
        request.applymarker(pytest.mark.flaky(reruns=5, reruns_delay=2))


def pytest_itemcollected(item: pytest.Item) -> None:
    # Dask AnnData tests require anndata > 0.10
    import anndata
    from packaging.version import Version

    requires_anndata_dask_support = (
        len(list(item.iter_markers(name="anndata_dask_support"))) > 0
    )

    if requires_anndata_dask_support and Version(anndata.__version__) < Version("0.10"):
        item.add_marker(
            pytest.mark.skip(reason="dask support requires anndata version > 0.10")
        )


assert "scanpy" not in sys.modules, (
    "scanpy is already imported, this will mess up test coverage"
)
