"""
This file contains helper functions for the scanpy test suite.
"""

from itertools import permutations

import scanpy as sc
import numpy as np

from anndata.tests.helpers import asarray, assert_equal

###########################
# Representation choice
###########################
# These functions can be used to check that functions are correctly using arugments like `layers`, `obsm`, etc.


def check_rep_mutation(func, X, fields=["layer", "obsm"], **kwargs):
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


def check_rep_results(func, X, **kwargs):
    """Checks that the results of a computation add values/ mutate the anndata object in a consistent way."""
    # Gen data
    adata_X = sc.AnnData(
        X=X.copy(),
        layers={"layer": np.zeros(shape=X.shape, dtype=X.dtype)},
        obsm={"obsm": np.zeros(shape=X.shape, dtype=X.dtype)},
    )
    adata_layer = sc.AnnData(
        X=np.zeros(shape=X.shape, dtype=X.dtype),
        layers={"layer": X.copy()},
        obsm={"obsm": np.zeros(shape=X.shape, dtype=X.dtype)},
    )
    adata_obsm = sc.AnnData(
        X=np.zeros(shape=X.shape, dtype=X.dtype),
        layers={"layer": np.zeros(shape=X.shape, dtype=X.dtype)},
        obsm={"obsm": X.copy()},
    )

    # Apply function
    func(adata_X, **kwargs)
    func(adata_layer, layer="layer", **kwargs)
    func(adata_obsm, obsm="obsm", **kwargs)

    # Reset X
    adata_X.X = np.zeros(shape=X.shape, dtype=X.dtype)
    adata_layer.layers["layer"] = np.zeros(shape=X.shape, dtype=X.dtype)
    adata_obsm.obsm["obsm"] = np.zeros(shape=X.shape, dtype=X.dtype)

    # Check equality
    assert_equal(adata_X, adata_layer)
    assert_equal(adata_X, adata_obsm)
