"""A private pytest plugin"""
from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pytest

from .fixtures import *  # noqa: F403


doctest_env_marker = pytest.mark.usefixtures('doctest_env')


# Defining it here because it’s autouse.
@pytest.fixture(autouse=True)
def test_context():
    """Switch to agg backend and close all figures at teardown."""
    from matplotlib import pyplot
    import scanpy as sc

    old_backend = pyplot.rcParams['backend']

    pyplot.switch_backend('agg')
    with sc.settings.verbosity.override('hint'):
        yield
    pyplot.close('all')
    pyplot.switch_backend(old_backend)


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


def pytest_itemcollected(item: pytest.Item) -> None:
    import pytest

    if not isinstance(item, pytest.DoctestItem):
        return

    item.add_marker(doctest_env_marker)

    func = _import_name(item.name)
    if marker := getattr(func, '_doctest_mark', None):
        item.add_marker(marker)
    if skip_reason := getattr(func, '_doctest_skip_reason', False):
        item.add_marker(pytest.mark.skip(reason=skip_reason))


def _import_name(name: str) -> Any:
    from importlib import import_module

    parts = name.split('.')
    obj = import_module(parts[0])
    for i, name in enumerate(parts[1:]):
        try:
            obj = import_module(f'{obj.__name__}.{name}')
        except ModuleNotFoundError:
            break
    for name in parts[i + 1 :]:
        try:
            obj = getattr(obj, name)
        except AttributeError:
            raise RuntimeError(f'{parts[:i]}, {parts[i+1:]}, {obj} {name}')
    return obj
