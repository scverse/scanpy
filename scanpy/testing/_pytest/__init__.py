"""A private pytest plugin"""
import pytest

from .fixtures import *  # noqa: F403


# In case pytest-nunit is not installed, defines a dummy fixture
try:
    import pytest_nunit
except ModuleNotFoundError:

    @pytest.fixture
    def add_nunit_attachment(request):
        def noop(file, description):
            pass

        return noop

    def pytest_addoption(parser):
        add_internet_tests_option(parser)
        parser.addini("nunit_attach_on", "Dummy nunit replacement", default="any")

else:

    def pytest_addoption(parser):
        add_internet_tests_option(parser)


def add_internet_tests_option(parser):
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help=(
            "Run tests that retrieve stuff from the internet. "
            "This increases test time."
        ),
    )


def pytest_collection_modifyitems(config, items):
    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--run-internet` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)
