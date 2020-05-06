"""
This file contains helper functions for the scanpy test suite.
"""

import scanpy as sc
import numpy as np

from anndata.tests.helpers import asarray, assert_equal

###########################
# Representation choice
###########################
# These functions can be used to check that functions are correctly using arugments like `layers`, `obsm`, etc.


def check_rep_mutation(func, X, **kwargs):
    """Check that only the array meant to be modified is modified."""
    adata = sc.AnnData(
        X=X.copy(),
        layers={"layer": X.copy()},
        obsm={"obsm": X.copy()},
        dtype=np.float64,
    )
    adata_X = func(adata, copy=True, **kwargs)
    adata_layer = func(adata, layer="layer", copy=True, **kwargs)
    adata_obsm = func(adata, obsm="obsm", copy=True, **kwargs)

    assert np.array_equal(asarray(adata_X.X), asarray(adata_layer.layers["layer"]))
    assert np.array_equal(asarray(adata_X.X), asarray(adata_obsm.obsm["obsm"]))

    assert np.array_equal(asarray(adata_layer.X), asarray(adata_layer.obsm["obsm"]))
    assert np.array_equal(asarray(adata_obsm.X), asarray(adata_obsm.layers["layer"]))
    assert np.array_equal(
        asarray(adata_X.layers["layer"]), asarray(adata_X.obsm["obsm"])
    )


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
