"""Helper functions for the scanpy test suite."""

from __future__ import annotations

import warnings
from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass
from importlib.util import find_spec
from itertools import permutations
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
from anndata.tests.helpers import asarray, assert_equal

import scanpy as sc

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableSequence

    from numpy.typing import NDArray

    from scanpy._compat import DaskArray

# TODO: Report more context on the fields being compared on error
# TODO: Allow specifying paths to ignore on comparison

###########################
# Representation choice
###########################
# These functions can be used to check that functions are correctly using arugments like `layers`, `obsm`, etc.


def anndata_v0_8_constructor_compat(X, *args, **kwargs):
    """Construct AnnData that uses dtype of X for test compatibility with older AnnData versions.

    Once the minimum version of AnnData is 0.9, this function can be replaced with the default constructor.
    """
    import anndata as ad
    from packaging.version import Version

    if Version(ad.__version__) < Version("0.9"):
        return ad.AnnData(X, *args, **kwargs, dtype=X.dtype)
    else:
        return ad.AnnData(X, *args, **kwargs)


def check_rep_mutation(func, X, *, fields=("layer", "obsm"), **kwargs):
    """Check that only the array meant to be modified is modified."""
    adata = anndata_v0_8_constructor_compat(X.copy())

    for field in fields:
        sc.get._set_obs_rep(adata, X, **{field: field})
    X_array = asarray(X)

    adata_X = func(adata, copy=True, **kwargs)
    adatas_proc = {
        field: func(adata, copy=True, **{field: field}, **kwargs) for field in fields
    }

    # Modified fields
    for field in fields:
        result_array = asarray(
            sc.get._get_obs_rep(adatas_proc[field], **{field: field})
        )
        np.testing.assert_array_equal(asarray(adata_X.X), result_array)

    # Unmodified fields
    for field in fields:
        np.testing.assert_array_equal(X_array, asarray(adatas_proc[field].X))
        np.testing.assert_array_equal(
            X_array, asarray(sc.get._get_obs_rep(adata_X, **{field: field}))
        )
    for field_a, field_b in permutations(fields, 2):
        result_array = asarray(
            sc.get._get_obs_rep(adatas_proc[field_a], **{field_b: field_b})
        )
        np.testing.assert_array_equal(X_array, result_array)


def check_rep_results(func, X, *, fields: Iterable[str] = ("layer", "obsm"), **kwargs):
    """Check that the results of a computation add values/ mutate the anndata object in a consistent way."""
    # Gen data
    empty_X = np.zeros(shape=X.shape, dtype=X.dtype)
    adata = sc.AnnData(
        X=empty_X.copy(),
        layers={"layer": empty_X.copy()},
        obsm={"obsm": empty_X.copy()},
    )

    adata_X = adata.copy()
    adata_X.X = X.copy()

    adatas_proc = {}
    for field in fields:
        cur = adata.copy()
        sc.get._set_obs_rep(cur, X.copy(), **{field: field})
        adatas_proc[field] = cur

    # Apply function
    func(adata_X, **kwargs)
    for field in fields:
        func(adatas_proc[field], **{field: field}, **kwargs)

    # Reset X
    adata_X.X = empty_X.copy()
    for field in fields:
        sc.get._set_obs_rep(adatas_proc[field], empty_X.copy(), **{field: field})

    for field_a, field_b in permutations(fields, 2):
        assert_equal(adatas_proc[field_a], adatas_proc[field_b])
    for field in fields:
        assert_equal(adata_X, adatas_proc[field])


def _check_check_values_warnings(
    function, adata, expected_warning, kwargs=MappingProxyType({})
):
    """Run `function` on `adata` with provided arguments `kwargs` twice.

    Once with `check_values=True` and once with `check_values=False`.
    Checks that the `expected_warning` is only raised whtn `check_values=True`.
    """
    # expecting 0 no-int warnings
    with warnings.catch_warnings(record=True) as record:
        function(adata.copy(), **kwargs, check_values=False)
    warning_msgs = [w.message.args[0] for w in record]
    assert expected_warning not in warning_msgs

    # expecting 1 no-int warning
    with warnings.catch_warnings(record=True) as record:
        function(adata.copy(), **kwargs, check_values=True)
    warning_msgs = [w.message.args[0] for w in record]
    assert expected_warning in warning_msgs


# Delayed imports for case where we aren't using dask
def as_dense_dask_array(*args, **kwargs) -> DaskArray:
    from anndata.tests.helpers import as_dense_dask_array

    return as_dense_dask_array(*args, **kwargs)


def as_sparse_dask_array(*args, **kwargs) -> DaskArray:
    from anndata.tests.helpers import as_sparse_dask_array

    return as_sparse_dask_array(*args, **kwargs)


@dataclass(init=False)
class MultiContext(AbstractContextManager):
    contexts: MutableSequence[AbstractContextManager]

    def __init__(self, *contexts: AbstractContextManager):
        self.contexts = list(contexts)

    def __enter__(self):
        for ctx in self.contexts:
            ctx.__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        for ctx in reversed(self.contexts):
            ctx.__exit__(exc_type, exc_value, traceback)


@contextmanager
def maybe_dask_process_context():
    """Switch to a single-threaded scheduler for tests that use numba.

    Running numba with dask's threaded scheduler causes crashes,
    so we need to switch to single-threaded (or processes, which is slower).
    """
    if not find_spec("dask"):
        yield
        return

    import dask.config

    prev_scheduler = dask.config.get("scheduler", "threads")
    dask.config.set(scheduler="single-threaded")
    try:
        yield
    finally:
        dask.config.set(scheduler=prev_scheduler)


def random_mask(n: int, *, rng: np.random.Generator | None = None) -> NDArray[np.bool_]:
    """Generate a random mask.

    Makes sure that at least 2 mask entries are True and at least 2 are False.
    This avoids off-by-1 errors even in e.g. neighbors (which already cuts 1 off).
    """
    assert n >= 4, "n must be at least 4"
    rng = np.random.default_rng(rng)
    mask = rng.choice([True, False], n)
    if (n_false := (~mask).sum()) < 2:
        mask[rng.choice(np.flatnonzero(mask), 2 - n_false, replace=False)] = False
    if (n_true := mask.sum()) < 2:
        mask[rng.choice(np.flatnonzero(~mask), 2 - n_true, replace=False)] = True
    return mask
