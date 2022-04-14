import sys
from pathlib import Path

import matplotlib as mpl

mpl.use('agg')
from matplotlib import pyplot
from matplotlib.testing.compare import compare_images, make_test_filename
import pytest

import scanpy


scanpy.settings.verbosity = "hint"

# define this after importing scanpy but before running tests
IMPORTED = frozenset(sys.modules.keys())


@pytest.fixture(autouse=True)
def close_figures_on_teardown():
    yield
    pyplot.close("all")


def clear_loggers():
    """Remove handlers from all loggers

    Fixes: https://github.com/scverse/scanpy/issues/1736

    Code from: https://github.com/pytest-dev/pytest/issues/5502#issuecomment-647157873
    """
    import logging

    loggers = [logging.getLogger()] + list(logging.Logger.manager.loggerDict.values())
    for logger in loggers:
        handlers = getattr(logger, 'handlers', [])
        for handler in handlers:
            logger.removeHandler(handler)


@pytest.fixture(scope="session", autouse=True)
def close_logs_on_teardown(request):
    request.addfinalizer(clear_loggers)


@pytest.fixture
def imported_modules():
    return IMPORTED


@pytest.fixture
def check_same_image(add_nunit_attachment):
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
            diff_pth = make_test_filename(pth2, 'failed-diff')
            add_nunit_attachment(str(pth1), fmt_descr("Expected"))
            add_nunit_attachment(str(pth2), fmt_descr("Result"))
            if Path(diff_pth).is_file():
                add_nunit_attachment(str(diff_pth), fmt_descr("Difference"))
            raise e

    return _


@pytest.fixture
def image_comparer(check_same_image):
    def make_comparer(path_expected: Path, path_actual: Path, *, tol: int):
        def save_and_compare(basename, tol=tol):
            path_actual.mkdir(parents=True, exist_ok=True)
            out_path = path_actual / f'{basename}.png'
            pyplot.savefig(out_path, dpi=40)
            pyplot.close()
            check_same_image(path_expected / f'{basename}.png', out_path, tol=tol)

        return save_and_compare

    return make_comparer


@pytest.fixture
def plt():
    return pyplot
