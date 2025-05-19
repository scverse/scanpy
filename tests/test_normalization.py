from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from anndata.tests.helpers import assert_equal
from scipy import sparse

import scanpy as sc
from scanpy._compat import CSBase
from scanpy._utils import axis_sum
from scanpy.preprocessing._normalization import _compute_nnz_median
from testing.scanpy._helpers import (
    _check_check_values_warnings,
    check_rep_mutation,
    check_rep_results,
)

# TODO: Add support for sparse-in-dask
from testing.scanpy._pytest.params import ARRAY_TYPES, ARRAY_TYPES_DENSE

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

X_total = np.array([[1, 0], [3, 0], [5, 6]])
X_frac = np.array([[1, 0, 1], [3, 0, 1], [5, 6, 1]])


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
@pytest.mark.parametrize("target_sum", [None, 1.0])
@pytest.mark.parametrize("exclude_highly_expressed", [True, False])
def test_normalize_matrix_types(
    array_type, dtype, target_sum, exclude_highly_expressed
):
    adata = sc.datasets.pbmc68k_reduced()
    adata.X = (adata.raw.X).astype(dtype)
    adata_casted = adata.copy()
    adata_casted.X = array_type(adata_casted.raw.X).astype(dtype)
    sc.pp.normalize_total(
        adata, target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed
    )
    sc.pp.normalize_total(
        adata_casted,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
    )
    X = adata_casted.X
    if "dask" in array_type.__name__:
        X = X.compute()
    if isinstance(X, CSBase):
        X = X.todense()
    if isinstance(adata.X, CSBase):
        adata.X = adata.X.todense()
    np.testing.assert_allclose(X, adata.X, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total(array_type, dtype):
    adata = AnnData(array_type(X_total).astype(dtype))
    sc.pp.normalize_total(adata, key_added="n_counts")
    assert np.allclose(np.ravel(axis_sum(adata.X, axis=1)), [3.0, 3.0, 3.0])
    sc.pp.normalize_total(adata, target_sum=1, key_added="n_counts2")
    assert np.allclose(np.ravel(axis_sum(adata.X, axis=1)), [1.0, 1.0, 1.0])

    adata = AnnData(array_type(X_frac).astype(dtype))
    sc.pp.normalize_total(adata, exclude_highly_expressed=True, max_fraction=0.7)
    assert np.allclose(np.ravel(axis_sum(adata.X[:, 1:3], axis=1)), [1.0, 1.0, 1.0])


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total_rep(array_type, dtype):
    # Test that layer kwarg works
    X = array_type(sparse.random(100, 50, format="csr", density=0.2, dtype=dtype))
    check_rep_mutation(sc.pp.normalize_total, X, fields=["layer"])
    check_rep_results(sc.pp.normalize_total, X, fields=["layer"])


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total_view(array_type, dtype):
    adata = AnnData(array_type(X_total).astype(dtype))
    v = adata[:, :]

    sc.pp.normalize_total(v)
    sc.pp.normalize_total(adata)

    assert not v.is_view
    assert_equal(adata, v)


def test_normalize_pearson_residuals_warnings(pbmc3k_parametrized):
    adata = pbmc3k_parametrized()

    if np.issubdtype(adata.X.dtype, np.integer):
        pytest.skip("Can’t store non-integral data with int dtype")

    # depending on check_values, warnings should be raised for non-integer data
    adata_noninteger = adata.copy()
    x, y = np.nonzero(adata_noninteger.X)
    adata_noninteger.X[x[0], y[0]] = 0.5

    _check_check_values_warnings(
        function=sc.experimental.pp.normalize_pearson_residuals,
        adata=adata_noninteger,
        expected_warning="`normalize_pearson_residuals()` expects raw count data, but non-integers were found.",
    )


@pytest.mark.parametrize(
    ("params", "match"),
    [
        pytest.param(dict(theta=0), r"Pearson residuals require theta > 0", id="theta"),
        pytest.param(
            dict(theta=-1), r"Pearson residuals require theta > 0", id="theta"
        ),
        pytest.param(
            dict(clip=-1), r"Pearson residuals require `clip>=0` or `clip=None`."
        ),
    ],
)
def test_normalize_pearson_residuals_errors(pbmc3k_parametrized, params, match):
    adata = pbmc3k_parametrized()

    with pytest.raises(ValueError, match=match):
        sc.experimental.pp.normalize_pearson_residuals(adata, **params)


@pytest.mark.parametrize(
    "sparsity_func",
    [np.array, sparse.csr_matrix],  # noqa: TID251
    ids=lambda x: x.__name__,
)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
@pytest.mark.parametrize("theta", [0.01, 1, 100, np.inf])
@pytest.mark.parametrize("clip", [None, 1, np.inf])
def test_normalize_pearson_residuals_values(sparsity_func, dtype, theta, clip):
    # toy data
    X = np.array([[3, 6], [2, 4], [1, 0]])
    ns = np.sum(X, axis=1)
    ps = np.sum(X, axis=0) / np.sum(X)
    mu = np.outer(ns, ps)

    # compute reference residuals
    if np.isinf(theta):
        # Poisson case
        residuals_reference = (X - mu) / np.sqrt(mu)
    else:
        # NB case
        residuals_reference = (X - mu) / np.sqrt(mu + mu**2 / theta)

    # compute output to test
    adata = AnnData(sparsity_func(X).astype(dtype))
    output = sc.experimental.pp.normalize_pearson_residuals(
        adata, theta=theta, clip=clip, inplace=False
    )
    output_X = output["X"]
    sc.experimental.pp.normalize_pearson_residuals(
        adata, theta=theta, clip=clip, inplace=True
    )

    # check for correct new `adata.uns` keys
    assert {"pearson_residuals_normalization"} <= adata.uns.keys()
    assert {"theta", "clip", "computed_on"} <= adata.uns[
        "pearson_residuals_normalization"
    ].keys()
    # test against inplace
    np.testing.assert_array_equal(adata.X, output_X)

    if clip is None:
        # default clipping: compare to sqrt(n) threshold
        clipping_threshold = np.sqrt(adata.shape[0]).astype(np.float32)
        assert np.max(output_X) <= clipping_threshold
        assert np.min(output_X) >= -clipping_threshold
    elif np.isinf(clip):
        # no clipping: compare to raw residuals
        assert np.allclose(output_X, residuals_reference)
    else:
        # custom clipping: compare to custom threshold
        assert np.max(output_X) <= clip
        assert np.min(output_X) >= -clip


def _check_pearson_pca_fields(ad, n_cells, n_comps):
    assert {"pearson_residuals_normalization", "pca"} <= ad.uns.keys(), (
        "Missing `.uns` keys. Expected `['pearson_residuals_normalization', 'pca']`, "
        f"but only {list(ad.uns.keys())} were found"
    )
    assert "X_pca" in ad.obsm, (
        f"Missing `obsm` key `'X_pca'`, only {list(ad.obsm.keys())} were found"
    )
    assert "PCs" in ad.varm, (
        f"Missing `varm` key `'PCs'`, only {list(ad.varm.keys())} were found"
    )
    assert ad.obsm["X_pca"].shape == (
        n_cells,
        n_comps,
    ), "Wrong shape of PCA output in `X_pca`"


@pytest.mark.parametrize("n_hvgs", [100, 200])
@pytest.mark.parametrize("n_comps", [30, 50])
@pytest.mark.parametrize(
    ("do_hvg", "params", "n_var_copy_name"),
    [
        pytest.param(False, dict(), "n_genes", id="no_hvg"),
        pytest.param(True, dict(), "n_hvgs", id="hvg_default"),
        pytest.param(
            True, dict(use_highly_variable=False), "n_genes", id="hvg_opt_out"
        ),
        pytest.param(False, dict(mask_var="test_mask"), "n_unmasked", id="mask"),
    ],
)
def test_normalize_pearson_residuals_pca(
    *,
    pbmc3k_parametrized_small: Callable[[], AnnData],
    n_hvgs: int,
    n_comps: int,
    do_hvg: bool,
    params: dict[str, Any],
    n_var_copy_name: str,  # number of variables in output if inplace=False
):
    adata = pbmc3k_parametrized_small()
    n_cells, n_genes = adata.shape
    n_unmasked = n_genes - 5
    adata.var["test_mask"] = np.r_[
        np.repeat(True, n_unmasked), np.repeat(False, n_genes - n_unmasked)  # noqa: FBT003
    ]
    n_var_copy = locals()[n_var_copy_name]
    assert isinstance(n_var_copy, int | np.integer)

    if do_hvg:
        sc.experimental.pp.highly_variable_genes(
            adata, flavor="pearson_residuals", n_top_genes=n_hvgs
        )

    # inplace=False
    adata_pca = sc.experimental.pp.normalize_pearson_residuals_pca(
        adata.copy(), inplace=False, n_comps=n_comps, **params
    )
    # inplace=True modifies the input adata object
    sc.experimental.pp.normalize_pearson_residuals_pca(
        adata, inplace=True, n_comps=n_comps, **params
    )

    for ad, n_var_ret in (
        (adata_pca, n_var_copy),
        # inplace adatas should always retains original shape
        (adata, n_genes),
    ):
        _check_pearson_pca_fields(ad, n_cells, n_comps)

        # check adata shape to see if all genes or only HVGs are in the returned adata
        assert ad.shape == (n_cells, n_var_ret)

        # check PC shapes to see whether or not HVGs were used for PCA
        assert ad.varm["PCs"].shape == (n_var_ret, n_comps)

    # check if there are columns of all-zeros in the PCs shapes
    # to see whether or not HVGs were used for PCA
    # either no all-zero-colums or all number corresponding to non-hvgs should exist
    assert sum(np.sum(np.abs(adata.varm["PCs"]), axis=1) == 0) == (n_genes - n_var_copy)

    # compare PCA results beteen inplace / copied
    np.testing.assert_array_equal(adata.obsm["X_pca"], adata_pca.obsm["X_pca"])


@pytest.mark.parametrize("n_hvgs", [100, 200])
@pytest.mark.parametrize("n_comps", [30, 50])
def test_normalize_pearson_residuals_recipe(pbmc3k_parametrized_small, n_hvgs, n_comps):
    adata = pbmc3k_parametrized_small()
    n_cells, n_genes = adata.shape

    ### inplace = False ###
    # outputs the (potentially hvg-restricted) adata_pca object
    # PCA on all genes
    adata_pca, hvg = sc.experimental.pp.recipe_pearson_residuals(
        adata.copy(), inplace=False, n_comps=n_comps, n_top_genes=n_hvgs
    )

    # check PCA fields
    _check_pearson_pca_fields(adata_pca, n_cells, n_comps)
    # check adata output shape (only HVGs in output)
    assert adata_pca.shape == (n_cells, n_hvgs)
    # check PC shape (non-hvgs are removed, so only `n_hvgs` genes)
    assert adata_pca.varm["PCs"].shape == (n_hvgs, n_comps)

    # check hvg df
    assert {
        "means",
        "variances",
        "residual_variances",
        "highly_variable_rank",
        "highly_variable",
    } <= set(hvg.columns)
    assert np.sum(hvg["highly_variable"]) == n_hvgs
    assert hvg.shape[0] == n_genes

    ### inplace = True ###
    # modifies the input adata object
    # PCA on all genes
    sc.experimental.pp.recipe_pearson_residuals(
        adata, inplace=True, n_comps=n_comps, n_top_genes=n_hvgs
    )

    # check PCA fields and output shape
    _check_pearson_pca_fields(adata, n_cells, n_comps)
    # check adata shape (no change to input)
    assert adata.shape == (n_cells, n_genes)
    # check PC shape (non-hvgs are masked with 0s, so original number of genes)
    assert adata.varm["PCs"].shape == (n_genes, n_comps)
    # number of all-zero-colums should be number of non-hvgs
    assert sum(np.sum(np.abs(adata.varm["PCs"]), axis=1) == 0) == n_genes - n_hvgs


@pytest.mark.parametrize("array_type", ARRAY_TYPES_DENSE)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_compute_nnz_median(array_type, dtype):
    data = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=dtype)
    data = array_type(data)
    np.testing.assert_allclose(_compute_nnz_median(data), 5)
