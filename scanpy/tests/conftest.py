from pathlib import Path

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot
from matplotlib.testing.compare import compare_images
import pytest

import scanpy

scanpy.settings.verbosity = "hint"


def make_comparer(path_expected: Path, path_actual: Path, *, tol: int):
    def save_and_compare(basename, tolerance=None):
        path_actual.mkdir(parents=True, exist_ok=True)
        out_path = path_actual / f'{basename}.png'
        pyplot.savefig(out_path, dpi=40)
        pyplot.close()
        if tolerance is None:
            tolerance = tol
        res = compare_images(str(path_expected / f'{basename}.png'), str(out_path), tolerance)
        assert res is None, res
    return save_and_compare


@pytest.fixture
def image_comparer():
    return make_comparer


@pytest.fixture
def plt():
    return pyplot
