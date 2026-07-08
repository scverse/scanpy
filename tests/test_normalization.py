from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from anndata.tests.helpers import assert_equal
from fast_array_utils import conv, stats
from scipy import sparse

import scanpy as sc
from scanpy._compat import CSBase
from scanpy.preprocessing._normalization import _compute_nnz_median
from testing.scanpy._helpers import (
    _check_check_values_warnings,
    check_rep_mutation,
    check_rep_results,
)
from testing.scanpy._pytest.marks import needs

# TODO: Add support for sparse-in-dask
from testing.scanpy._pytest.params import (
    ARRAY_TYPES,
    ARRAY_TYPES_DENSE,
    ARRAY_TYPES_MEM,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

to_ndarray = partial(conv.to_dense, to_cpu_memory=True)

X_total = np.array([[1, 0], [3, 0], [5, 6]])
X_frac = np.array([[1, 0, 1], [3, 0, 1], [5, 6, 1]])


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
@pytest.mark.parametrize("target_sum", [None, 1.0], ids=["no_target_sum", "target_sum"])
@pytest.mark.parametrize(
    "exclude_highly_expressed", [True, False], ids=["excl_hi", "no_excl_hi"]
)
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
    adata.X = conv.to_dense(adata.X)
    adata_casted.X = conv.to_dense(adata_casted.X, to_cpu_memory=True)
    np.testing.assert_allclose(adata_casted.X, adata.X, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total(array_type, dtype):
    adata = AnnData(array_type(X_total).astype(dtype))
    sc.pp.normalize_total(adata, key_added="n_counts")
    assert np.allclose(to_ndarray(stats.sum(adata.X, axis=1)), [3.0, 3.0, 3.0])
    sc.pp.normalize_total(adata, target_sum=1, key_added="n_counts2")
    assert np.allclose(to_ndarray(stats.sum(adata.X, axis=1)), [1.0, 1.0, 1.0])

    adata = AnnData(array_type(X_frac).astype(dtype))
    sc.pp.normalize_total(adata, exclude_highly_expressed=True, max_fraction=0.7)
    assert np.allclose(to_ndarray(stats.sum(adata.X[:, 1:3], axis=1)), [1.0, 1.0, 1.0])


@pytest.mark.filterwarnings("ignore:Some cells have zero counts:UserWarning")
@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total_rep(array_type, dtype):
    """Test that layer/obsm kwargs work."""
    x = array_type(sparse.random(100, 50, format="csr", density=0.2, dtype=dtype))
    check_rep_mutation(sc.pp.normalize_total, x)
    check_rep_results(sc.pp.normalize_total, x)


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_total_view(array_type, dtype):
    adata = AnnData(array_type(X_total).astype(dtype))
    v = adata[:, :]

    with pytest.warns(UserWarning, match=r"Received a view"):
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
        pytest.param(
            dict(theta=0), r"Pearson residuals require theta > 0", id="theta=0"
        ),
        pytest.param(
            dict(theta=-1), r"Pearson residuals require theta > 0", id="theta=-1"
        ),
        pytest.param(
            dict(clip=-1),
            r"Pearson residuals require `clip>=0` or `clip=None`.",
            id="clip=-1",
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
    x = np.array([[3, 6], [2, 4], [1, 0]])
    ns = np.sum(x, axis=1)
    ps = np.sum(x, axis=0) / np.sum(x)
    mu = np.outer(ns, ps)

    # compute reference residuals
    if np.isinf(theta):
        # Poisson case
        residuals_reference = (x - mu) / np.sqrt(mu)
    else:
        # NB case
        residuals_reference = (x - mu) / np.sqrt(mu + mu**2 / theta)

    # compute output to test
    adata = AnnData(sparsity_func(x).astype(dtype))
    output = sc.experimental.pp.normalize_pearson_residuals(
        adata, theta=theta, clip=clip, inplace=False
    )
    output_x = output["X"]
    sc.experimental.pp.normalize_pearson_residuals(
        adata, theta=theta, clip=clip, inplace=True
    )

    # check for correct new `adata.uns` keys
    assert {"pearson_residuals_normalization"} <= adata.uns.keys()
    assert {"theta", "clip", "computed_on"} <= adata.uns[
        "pearson_residuals_normalization"
    ].keys()
    # test against inplace
    np.testing.assert_array_equal(adata.X, output_x)

    if clip is None:
        # default clipping: compare to sqrt(n) threshold
        clipping_threshold = np.sqrt(adata.shape[0]).astype(np.float32)
        assert np.max(output_x) <= clipping_threshold
        assert np.min(output_x) >= -clipping_threshold
    elif np.isinf(clip):
        # no clipping: compare to raw residuals
        assert np.allclose(output_x, residuals_reference)
    else:
        # custom clipping: compare to custom threshold
        assert np.max(output_x) <= clip
        assert np.min(output_x) >= -clip


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
        pytest.param(True, dict(mask_var=None), "n_genes", id="hvg_opt_out"),
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
def test_normalize_pearson_residuals_recipe(
    pbmc3k_parametrized_small: Callable[[], AnnData], n_hvgs: int, n_comps: int
) -> None:
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


# ------------------------------------------------------------------------------
# normalize_clr (shifted CLR / PFlog)
# ------------------------------------------------------------------------------

# A small count matrix with no empty cells, used for the value/equivalence tests.
X_clr = np.array(
    [[5, 0, 3, 2], [1, 1, 0, 4], [0, 7, 2, 1], [3, 3, 3, 3]], dtype="float32"
)


def _estimate_alpha_reference(x) -> float:
    """Calculate OLS overdispersion for reference."""
    x = np.asarray(to_ndarray(x), dtype=np.float64)
    mu = x.mean(axis=0)
    var = (x**2).mean(axis=0) - mu**2
    mu2 = mu**2
    return float(np.sum((var - mu) * mu2) / np.sum(mu2 * mu2))


def _clr_reference(x, *, target="auto", alpha=None) -> np.ndarray:
    """Calculate shifted CLR densely for reference.

    PFlog uses a constant log1p(4 * alpha * x) scale. Depth targets use the
    fixed-target scale log1p(K * x / depth). Empty cells are left as all-zero rows.
    """
    x = np.asarray(to_ndarray(x), dtype=np.float64)
    depths = x.sum(axis=1)
    if alpha is not None:
        log_u = np.log1p(4.0 * alpha * x)
    elif target == "auto":
        log_u = np.log1p(4.0 * _estimate_alpha_reference(x) * x)
    else:
        if target == "mean":
            target_sum = depths.mean()
        elif target == "median":
            target_sum = np.median(depths)
        else:
            target_sum = float(target)
        safe_depths = np.where(depths == 0, 1.0, depths)
        log_u = np.log1p(x * (target_sum / safe_depths)[:, None])
    return log_u - log_u.mean(axis=1, keepdims=True)


def _reconstruct_clr(adata, key="pflog") -> np.ndarray:
    return (
        to_ndarray(adata.layers[key]) - adata.obs[f"{key}_center"].to_numpy()[:, None]
    )


def _materialize(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return to_ndarray(x)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
@pytest.mark.parametrize("dtype", ["float32", "int64"])
def test_normalize_clr_values(array_type, dtype):
    """Check values against the reference and zero-sum cells."""
    adata = AnnData(array_type(X_clr).astype(dtype))
    sc.pp.normalize_clr(adata)
    result = _reconstruct_clr(adata)

    np.testing.assert_allclose(result, _clr_reference(X_clr), rtol=1e-5, atol=1e-5)
    # zero-sum (Aitchison) hyperplane
    np.testing.assert_allclose(result.sum(axis=1), 0.0, atol=1e-5)
    assert adata.layers["pflog"].nnz == sparse.csr_matrix(X_clr).nnz  # noqa: TID251
    assert adata.uns["pflog"]["encoding_type"] == "shifted_clr"
    assert adata.uns["pflog"]["row_center_key"] == "pflog_center"


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
@pytest.mark.parametrize(
    "kwargs",
    [{}, {"target": 1e4}, {"target": "mean"}, {"target": "median"}, {"alpha": 0.5}],
    ids=["default", "fixed_target", "mean_target", "median_target", "alpha"],
)
def test_normalize_clr_params(array_type, kwargs):
    adata = AnnData(array_type(X_clr).astype("float32"))
    sc.pp.normalize_clr(adata, **kwargs)
    np.testing.assert_allclose(
        _reconstruct_clr(adata), _clr_reference(X_clr, **kwargs), rtol=1e-5, atol=1e-5
    )


def test_normalize_clr_alpha_is_constant_scale_and_overrides_target():
    """`alpha` stores uncentered log1p(4 * alpha * x), independent of depth."""
    alpha = 0.5

    via_alpha = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    sc.pp.normalize_clr(via_alpha, alpha=alpha)
    np.testing.assert_allclose(
        via_alpha.layers["pflog"].toarray(),
        np.log1p(4.0 * alpha * X_clr),
        rtol=1e-5,
        atol=1e-5,
    )

    both = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    sc.pp.normalize_clr(both, alpha=alpha, target=999.0)
    np.testing.assert_allclose(
        _reconstruct_clr(both), _reconstruct_clr(via_alpha), rtol=1e-5, atol=1e-5
    )


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_normalize_clr_alpha_auto(array_type):
    """Check that `target="auto"` matches explicit alpha."""
    estimated = _estimate_alpha_reference(X_clr)
    assert estimated > 0

    auto = AnnData(array_type(X_clr).astype("float32"))
    sc.pp.normalize_clr(auto)

    explicit = AnnData(array_type(X_clr).astype("float32"))
    sc.pp.normalize_clr(explicit, alpha=estimated)
    np.testing.assert_allclose(
        _reconstruct_clr(auto), _reconstruct_clr(explicit), rtol=1e-5, atol=1e-5
    )


@pytest.mark.parametrize("alpha", [0.0, -0.5], ids=["zero", "negative"])
def test_normalize_clr_nonpositive_alpha_raises(alpha):
    """Raise for non-positive `alpha`."""
    adata = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    with pytest.raises(ValueError, match=r"alpha.*positive"):
        sc.pp.normalize_clr(adata, alpha=alpha)


def test_normalize_clr_alpha_auto_zero_mean_raises():
    """`target="auto"` cannot estimate overdispersion when every gene mean is zero."""
    adata = AnnData(np.zeros((3, 4), dtype="float32"))
    with pytest.raises(ValueError, match="Cannot estimate overdispersion"):
        sc.pp.normalize_clr(adata)


@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_normalize_clr_zero_cell(array_type):
    """Keep an empty cell finite and all-zero."""
    x = X_clr.copy()
    x[1] = 0  # make the second cell empty
    adata = AnnData(array_type(x))
    with pytest.warns(UserWarning, match="Some cells have zero counts"):
        sc.pp.normalize_clr(adata)
    result = _reconstruct_clr(adata)
    assert np.isfinite(result).all()
    np.testing.assert_allclose(result[1], 0.0, atol=1e-6)


def test_normalize_clr_inplace_false():
    adata = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    x_before = to_ndarray(adata.X).copy()
    out = sc.pp.normalize_clr(adata, inplace=False)

    assert isinstance(out, dict)
    np.testing.assert_allclose(
        to_ndarray(out["X"]) - out["row_center"][:, None],
        _clr_reference(X_clr),
        rtol=1e-5,
        atol=1e-5,
    )
    # input is left untouched
    assert isinstance(adata.X, CSBase)
    np.testing.assert_array_equal(to_ndarray(adata.X), x_before)


def test_normalize_clr_copy():
    adata = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    returned = sc.pp.normalize_clr(adata, copy=True)

    assert isinstance(returned, AnnData)
    assert returned is not adata
    np.testing.assert_allclose(
        _reconstruct_clr(returned), _clr_reference(X_clr), rtol=1e-5, atol=1e-5
    )
    # original is left untouched
    assert isinstance(adata.X, CSBase)


def test_normalize_clr_copy_inplace_error():
    adata = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    with pytest.raises(
        ValueError, match="`copy=True` cannot be used with `inplace=False`"
    ):
        sc.pp.normalize_clr(adata, copy=True, inplace=False)


def test_normalize_clr_layer():
    """`layer` selects the input layer and leaves `X` untouched."""
    adata = AnnData(
        sparse.csr_matrix(X_clr),  # noqa: TID251
        layers={"counts": sparse.csr_matrix(X_clr)},  # noqa: TID251
    )
    x_before = to_ndarray(adata.X).copy()
    sc.pp.normalize_clr(adata, layer="counts")

    np.testing.assert_array_equal(to_ndarray(adata.X), x_before)
    np.testing.assert_allclose(
        _reconstruct_clr(adata),
        _clr_reference(X_clr),
        rtol=1e-5,
        atol=1e-5,
    )


def test_normalize_clr_densify():
    adata = AnnData(sparse.csr_matrix(X_clr))  # noqa: TID251
    sc.pp.normalize_clr(adata, densify=True)
    np.testing.assert_allclose(to_ndarray(adata.X), _clr_reference(X_clr), rtol=1e-5)
    assert "pflog" in adata.layers


@needs.dask
@pytest.mark.parametrize("densify", [False, True], ids=["sparse_encoded", "densified"])
@pytest.mark.parametrize(
    "sparse_blocks", [False, True], ids=["dense_dask", "sparse_dask"]
)
def test_normalize_clr_dask(sparse_blocks, densify):
    import dask.array as da

    chunks = (2, X_clr.shape[1])
    x = (
        da.from_array(sparse.csr_matrix(X_clr), chunks=chunks, asarray=False)  # noqa: TID251
        if sparse_blocks
        else da.from_array(X_clr, chunks=chunks)
    )
    adata = AnnData(x.astype("float32"))

    sc.pp.normalize_clr(adata, densify=densify)

    result = (
        _materialize(adata.layers["pflog"])
        - adata.obs["pflog_center"].to_numpy()[:, None]
    )
    np.testing.assert_allclose(result, _clr_reference(X_clr), rtol=1e-5, atol=1e-5)
    if densify:
        np.testing.assert_allclose(
            _materialize(adata.X), _clr_reference(X_clr), rtol=1e-5, atol=1e-5
        )


def test_normalize_clr_view():
    adata = AnnData(X_clr.copy())
    v = adata[:, :]
    with pytest.warns(UserWarning, match=r"Received a view"):
        sc.pp.normalize_clr(v)
    assert not v.is_view
