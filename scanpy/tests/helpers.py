"""
This file contains helper functions for the scanpy test suite.
"""

from itertools import permutations

import scanpy as sc
import numpy as np
import warnings
import pytest
from anndata.tests.helpers import asarray, assert_equal
from scanpy.tests._data._cached_datasets import pbmc3k

# TODO: Report more context on the fields being compared on error
# TODO: Allow specifying paths to ignore on comparison

###########################
# Representation choice
###########################
# These functions can be used to check that functions are correctly using arugments like `layers`, `obsm`, etc.


def check_rep_mutation(func, X, *, fields=["layer", "obsm"], **kwargs):
    """Check that only the array meant to be modified is modified."""
    adata = sc.AnnData(X=X.copy(), dtype=X.dtype)
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


def check_rep_results(func, X, *, fields=["layer", "obsm"], **kwargs):
    """Checks that the results of a computation add values/ mutate the anndata object in a consistent way."""
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


def _prepare_pbmc_testdata(sparsity_func, dtype, small=False):
    """Prepares 3k PBMC dataset with batch key `batch` and defined datatype/sparsity.

    Params
    ------
    sparsity_func
        sparsity function applied to adata.X (e.g. csr_matrix.toarray for dense or csr_matrix for sparse)
    dtype
        numpy dtype applied to adata.X (e.g. 'float32' or 'int64')
    small
        False (default) returns full data, True returns small subset of the data."""

    adata = pbmc3k().copy()

    if small:
        adata = adata[:1000, :500]
        sc.pp.filter_cells(adata, min_genes=1)
    np.random.seed(42)
    adata.obs['batch'] = np.random.randint(0, 3, size=adata.shape[0])
    sc.pp.filter_genes(adata, min_cells=1)
    adata.X = sparsity_func(adata.X.astype(dtype))
    return adata


def _check_check_values_warnings(function, adata, expected_warning, kwargs={}):
    '''Runs `function` on `adata` with provided arguments `kwargs` twice: once with `check_values=True` and once with `check_values=False`. Checks that the `expected_warning` is only raised whtn `check_values=True`.'''

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
