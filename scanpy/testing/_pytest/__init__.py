"""A private pytest plugin"""
from __future__ import annotations

from collections.abc import Iterable
import re

import pytest

from scanpy.testing._pytest.marks import needs

from .fixtures import *  # noqa: F403


external_tool_to_mod = dict(
    harmony_integrate='harmonypy',
    harmony_timeseries='harmony',
    sam='samalg',
    sandbag='pypairs',
    scanorama_integrate='scanorama',
)


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
    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--run-internet` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)


def pytest_itemcollected(item: pytest.Item) -> None:
    if not isinstance(item, pytest.DoctestItem):
        return
    if not (
        match := re.fullmatch(r'scanpy\.external\..+\.(?P<name>[^\.]+)', item.name)
    ):
        return
    tool = match['name']
    if tool in {'hashsolo'}:
        return  # no deps
    module = external_tool_to_mod.get(tool, tool)
    item.add_marker(needs(module))
