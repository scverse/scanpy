from __future__ import annotations

import itertools
from pathlib import Path
from string import ascii_letters
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from pandas.testing import assert_frame_equal, assert_index_equal

import scanpy as sc
from scanpy._compat import CSRBase
from scanpy.preprocessing._utils import _get_mean_var
from testing.scanpy._helpers import _check_check_values_warnings
from testing.scanpy._helpers.data import pbmc3k, pbmc68k_reduced
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

FILE = Path(__file__).parent / Path("_scripts/seurat_hvg.csv")
FILE_V3 = Path(__file__).parent / Path("_scripts/seurat_hvg_v3.csv.gz")
FILE_V3_BATCH = Path(__file__).parent / Path("_scripts/seurat_hvg_v3_batch.csv")
FILE_CELL_RANGER = Path(__file__).parent / "_scripts/cell_ranger_hvg.csv"


@pytest.fixture(scope="session")
def adata_sess() -> AnnData:
    adata = sc.datasets.blobs()
    rng = np.random.default_rng(0)
    adata.var_names = rng.choice(list(ascii_letters), adata.n_vars, replace=False)
    return adata


@pytest.fixture
def adata(adata_sess: AnnData) -> AnnData:
    return adata_sess.copy()


def test_runs(adata):
    sc.pp.highly_variable_genes(adata)


def test_supports_batch(adata):
    gen = np.random.default_rng(0)
    adata.obs["batch"] = pd.array(
        gen.binomial(3, 0.5, size=adata.n_obs), dtype="category"
    )
    sc.pp.highly_variable_genes(adata, batch_key="batch")
    assert "highly_variable_nbatches" in adata.var.columns
    assert "highly_variable_intersection" in adata.var.columns


def test_supports_layers(adata_sess):
    def execute(layer: str | None) -> AnnData:
        gen = np.random.default_rng(0)
        adata = adata_sess.copy()
        assert isinstance(adata.X, np.ndarray)
        if layer:
            adata.X, adata.layers[layer] = None, adata.X.copy()
            gen.shuffle(adata.layers[layer])
        adata.obs["batch"] = pd.array(
            gen.binomial(4, 0.5, size=adata.n_obs), dtype="category"
        )
        sc.pp.highly_variable_genes(
            adata, batch_key="batch", n_top_genes=3, layer=layer
        )
        assert "highly_variable_nbatches" in adata.var.columns
        assert adata.var["highly_variable"].sum() == 3
        return adata

    adata1, adata2 = map(execute, [None, "test_layer"])
    assert (adata1.var["highly_variable"] != adata2.var["highly_variable"]).any()


def test_no_batch_matches_batch(adata):
    sc.pp.highly_variable_genes(adata)
    no_batch_hvg = adata.var["highly_variable"].copy()
    assert no_batch_hvg.any()
    adata.obs["batch"] = pd.array(["batch"], dtype="category").repeat(len(adata))
    sc.pp.highly_variable_genes(adata, batch_key="batch")
    assert np.all(no_batch_hvg == adata.var["highly_variable"])
    assert np.all(
        adata.var["highly_variable_intersection"] == adata.var["highly_variable"]
    )


@pytest.mark.parametrize("batch_key", [None, "batch"], ids=["single", "batched"])
@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_no_inplace(adata, array_type, batch_key):
    """Tests that, with `n_top_genes=None` the returned dataframe has the expected columns."""
    adata.X = array_type(adata.X)
    if batch_key:
        adata.obs[batch_key] = np.tile(["a", "b"], adata.shape[0] // 2)
    sc.pp.highly_variable_genes(adata, batch_key=batch_key, n_bins=3)
    assert adata.var["highly_variable"].any()

    colnames = {"means", "dispersions", "dispersions_norm", "highly_variable"} | (
        {"mean_bin"}
        if batch_key is None
        else {"highly_variable_nbatches", "highly_variable_intersection"}
    )
    hvg_df = sc.pp.highly_variable_genes(
        adata, batch_key=batch_key, n_bins=3, inplace=False
    )
    assert isinstance(hvg_df, pd.DataFrame)
    assert colnames == set(hvg_df.columns)


@pytest.mark.parametrize("base", [None, 10])
@pytest.mark.parametrize("flavor", ["seurat", "cell_ranger"])
def test_keep_layer(base, flavor):
    adata = pbmc3k()
    # cell_ranger flavor can raise error if many 0 genes
    sc.pp.filter_genes(adata, min_counts=1)

    sc.pp.log1p(adata, base=base)
    assert isinstance(adata.X, CSRBase)
    X_orig = adata.X.copy()

    if flavor == "seurat":
        sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor=flavor)
    elif flavor == "cell_ranger":
        sc.pp.highly_variable_genes(adata, flavor=flavor)
    else:
        pytest.fail(f"Unknown {flavor=}")

    assert np.allclose(X_orig.toarray(), adata.X.toarray())


@pytest.mark.parametrize(
    "flavor",
    [
        "seurat",
        pytest.param(
            "cell_ranger",
            marks=pytest.mark.xfail(reason="can’t deal with duplicate bin edges"),
        ),
    ],
)
def test_no_filter_genes(flavor):
    """Test that even with columns containing all-zeros in the data, n_top_genes is respected."""
    adata = sc.datasets.pbmc3k()
    means, _ = _get_mean_var(adata.X)
    assert (means == 0).any()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=10000)
    assert adata.var["highly_variable"].sum() == 10000


def _check_pearson_hvg_columns(output_df: pd.DataFrame, n_top_genes: int):
    assert pd.api.types.is_float_dtype(output_df["residual_variances"].dtype)

    assert output_df["highly_variable"].to_numpy().dtype is np.dtype("bool")
    assert np.sum(output_df["highly_variable"]) == n_top_genes

    assert np.nanmax(output_df["highly_variable_rank"].to_numpy()) <= n_top_genes - 1


def test_pearson_residuals_inputchecks(pbmc3k_parametrized_small):
    adata = pbmc3k_parametrized_small()

    # depending on check_values, warnings should be raised for non-integer data
    if adata.X.dtype == "float32":
        adata_noninteger = adata.copy()
        x, y = np.nonzero(adata_noninteger.X)
        adata_noninteger.X[x[0], y[0]] = 0.5

        _check_check_values_warnings(
            function=sc.experimental.pp.highly_variable_genes,
            adata=adata_noninteger,
            expected_warning="`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
            kwargs=dict(
                flavor="pearson_residuals",
                n_top_genes=100,
            ),
        )

    # errors should be raised for invalid theta values
    for theta in [0, -1]:
        with pytest.raises(ValueError, match="Pearson residuals require theta > 0"):
            sc.experimental.pp.highly_variable_genes(
                adata.copy(), theta=theta, flavor="pearson_residuals", n_top_genes=100
            )

    with pytest.raises(
        ValueError, match="Pearson residuals require `clip>=0` or `clip=None`."
    ):
        sc.experimental.pp.highly_variable_genes(
            adata.copy(), clip=-1, flavor="pearson_residuals", n_top_genes=100
        )


@pytest.mark.parametrize("subset", [True, False], ids=["subset", "full"])
@pytest.mark.parametrize(
    "clip", [None, np.inf, 30], ids=["noclip", "infclip", "30clip"]
)
@pytest.mark.parametrize("theta", [100, np.inf], ids=["100theta", "inftheta"])
@pytest.mark.parametrize("n_top_genes", [100, 200], ids=["100n", "200n"])
def test_pearson_residuals_general(
    pbmc3k_parametrized_small, subset, clip, theta, n_top_genes
):
    adata = pbmc3k_parametrized_small()
    # cleanup var
    del adata.var

    # compute reference output
    residuals_res = sc.experimental.pp.normalize_pearson_residuals(
        adata, clip=clip, theta=theta, inplace=False
    )
    assert isinstance(residuals_res, dict)
    residual_variances_reference = np.var(residuals_res["X"], axis=0)

    if subset:
        # lazyly sort by residual variance and take top N
        top_n_idx = np.argsort(-residual_variances_reference)[:n_top_genes]
        # (results in sorted "gene order" in reference)
        residual_variances_reference = residual_variances_reference[top_n_idx]

    # compute output to be tested
    output_df = sc.experimental.pp.highly_variable_genes(
        adata,
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=False,
        clip=clip,
        theta=theta,
    )
    assert output_df is not None

    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        subset=subset,
        inplace=True,
        clip=clip,
        theta=theta,
    )

    # compare inplace=True and inplace=False output
    pd.testing.assert_frame_equal(output_df, adata.var)

    # check output is complete
    for key in [
        "highly_variable",
        "means",
        "variances",
        "residual_variances",
        "highly_variable_rank",
    ]:
        assert key in output_df.columns

    # check consistency with normalization method
    if subset:
        # sort values before comparing as reference is sorted as well for subset case
        sort_output_idx = np.argsort(-output_df["residual_variances"].to_numpy())
        assert np.allclose(
            output_df["residual_variances"].to_numpy()[sort_output_idx],
            residual_variances_reference,
        )
    else:
        assert np.allclose(
            output_df["residual_variances"].to_numpy(), residual_variances_reference
        )

    # check hvg flag
    hvg_idx = np.where(output_df["highly_variable"])[0]
    topn_idx = np.sort(
        np.argsort(-output_df["residual_variances"].to_numpy())[:n_top_genes]
    )
    assert np.all(hvg_idx == topn_idx)

    # check ranks
    assert np.nanmin(output_df["highly_variable_rank"].to_numpy()) == 0

    # more general checks on ranks, hvg flag and residual variance
    _check_pearson_hvg_columns(output_df, n_top_genes)


@pytest.mark.parametrize("subset", [True, False], ids=["subset", "full"])
@pytest.mark.parametrize("n_top_genes", [100, 200], ids=["100n", "200n"])
def test_pearson_residuals_batch(pbmc3k_parametrized_small, subset, n_top_genes):
    adata = pbmc3k_parametrized_small()
    # cleanup var
    del adata.var
    n_genes = adata.shape[1]

    output_df = sc.experimental.pp.highly_variable_genes(
        adata,
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        batch_key="batch",
        subset=subset,
        inplace=False,
    )
    assert output_df is not None

    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor="pearson_residuals",
        n_top_genes=n_top_genes,
        batch_key="batch",
        subset=subset,
        inplace=True,
    )

    # compare inplace=True and inplace=False output
    pd.testing.assert_frame_equal(output_df, adata.var)

    # check output is complete
    for key in [
        "highly_variable",
        "means",
        "variances",
        "residual_variances",
        "highly_variable_rank",
        "highly_variable_nbatches",
        "highly_variable_intersection",
    ]:
        assert key in output_df.columns

    # general checks on ranks, hvg flag and residual variance
    _check_pearson_hvg_columns(output_df, n_top_genes)

    # check intersection flag
    nbatches = len(np.unique(adata.obs["batch"]))
    assert output_df["highly_variable_intersection"].to_numpy().dtype is np.dtype(
        "bool"
    )
    assert np.sum(output_df["highly_variable_intersection"]) <= n_top_genes * nbatches
    assert np.all(output_df["highly_variable"][output_df.highly_variable_intersection])

    # check ranks (with batch_key these are the median of within-batch ranks)
    assert pd.api.types.is_float_dtype(output_df["highly_variable_rank"].dtype)

    # check nbatches
    assert output_df["highly_variable_nbatches"].to_numpy().dtype is np.dtype("int")
    assert np.min(output_df["highly_variable_nbatches"].to_numpy()) >= 0
    assert np.max(output_df["highly_variable_nbatches"].to_numpy()) <= nbatches

    # check subsetting
    if subset:
        assert len(output_df) == n_top_genes
    else:
        assert len(output_df) == n_genes


@pytest.mark.parametrize("func", ["hvg", "fgd"])
@pytest.mark.parametrize(
    ("flavor", "params", "ref_path"),
    [
        pytest.param(
            "seurat", dict(min_mean=0.0125, max_mean=3, min_disp=0.5), FILE, id="seurat"
        ),
        pytest.param(
            "cell_ranger", dict(n_top_genes=100), FILE_CELL_RANGER, id="cell_ranger"
        ),
    ],
)
@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_compare_to_upstream(
    *,
    request: pytest.FixtureRequest,
    func: Literal["hvg", "fgd"],
    flavor: Literal["seurat", "cell_ranger"],
    params: dict[str, float | int],
    ref_path: Path,
    array_type: Callable,
):
    if func == "fgd" and flavor == "cell_ranger":
        reason = "The deprecated filter_genes_dispersion behaves differently with cell_ranger"
        request.applymarker(pytest.mark.xfail(reason=reason))
    hvg_info = pd.read_csv(ref_path)

    pbmc = pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    pbmc.X = array_type(pbmc.X)
    pbmc.var_names_make_unique()
    sc.pp.filter_cells(pbmc, min_counts=1)
    sc.pp.normalize_total(pbmc, target_sum=1e4)

    if func == "hvg":
        sc.pp.log1p(pbmc)
        sc.pp.highly_variable_genes(pbmc, flavor=flavor, **params, inplace=True)
    elif func == "fgd":
        sc.pp.filter_genes_dispersion(
            pbmc, flavor=flavor, **params, log=True, subset=False
        )
    else:
        raise AssertionError()

    np.testing.assert_array_equal(
        hvg_info["highly_variable"], pbmc.var["highly_variable"]
    )

    # (still) Not equal to tolerance rtol=2e-05, atol=2e-05
    # np.testing.assert_allclose(4, 3.9999, rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(
        hvg_info["means"],
        pbmc.var["means"],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        hvg_info["dispersions"],
        pbmc.var["dispersions"],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        hvg_info["dispersions_norm"],
        pbmc.var["dispersions_norm"],
        rtol=2e-05 if "dask" not in array_type.__name__ else 1e-4,
        atol=2e-05 if "dask" not in array_type.__name__ else 1e-4,
    )


@needs.skmisc
def test_compare_to_seurat_v3():
    ### test without batch
    seurat_hvg_info = pd.read_csv(FILE_V3)

    pbmc = pbmc3k()
    sc.pp.filter_cells(pbmc, min_genes=200)  # this doesnt do anything btw
    sc.pp.filter_genes(pbmc, min_cells=3)

    pbmc_dense = pbmc.copy()
    pbmc_dense.X = pbmc_dense.X.toarray()

    sc.pp.highly_variable_genes(pbmc, n_top_genes=1000, flavor="seurat_v3")
    sc.pp.highly_variable_genes(pbmc_dense, n_top_genes=1000, flavor="seurat_v3")

    np.testing.assert_allclose(
        seurat_hvg_info["variance"],
        pbmc.var["variances"],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        seurat_hvg_info["variance.standardized"],
        pbmc.var["variances_norm"],
        rtol=2e-05,
        atol=2e-05,
    )
    np.testing.assert_allclose(
        pbmc_dense.var["variances_norm"],
        pbmc.var["variances_norm"],
        rtol=2e-05,
        atol=2e-05,
    )

    ### test with batch
    # introduce a dummy "technical covariate"; this is used in Seurat's SelectIntegrationFeatures
    pbmc.obs["dummy_tech"] = (
        "source_" + pd.array([*range(1, 6), 5]).repeat(500).astype("string")
    )[: pbmc.n_obs]

    seurat_v3_paper = sc.pp.highly_variable_genes(
        pbmc,
        n_top_genes=2000,
        flavor="seurat_v3_paper",
        batch_key="dummy_tech",
        inplace=False,
    )

    seurat_v3 = sc.pp.highly_variable_genes(
        pbmc,
        n_top_genes=2000,
        flavor="seurat_v3",
        batch_key="dummy_tech",
        inplace=False,
    )

    seurat_hvg_info_batch = pd.read_csv(FILE_V3_BATCH)
    seu = pd.Index(seurat_hvg_info_batch["x"].to_numpy())

    gene_intersection_paper = seu.intersection(
        seurat_v3_paper[seurat_v3_paper["highly_variable"]].index
    )
    gene_intersection_impl = seu.intersection(
        seurat_v3[seurat_v3["highly_variable"]].index
    )
    assert len(gene_intersection_paper) / 2000 > 0.95
    assert len(gene_intersection_impl) / 2000 < 0.95


@needs.skmisc
def test_seurat_v3_warning():
    pbmc = pbmc3k()[:200].copy()
    sc.pp.log1p(pbmc)
    with pytest.warns(
        UserWarning,
        match="`flavor='seurat_v3'` expects raw count data, but non-integers were found.",
    ):
        sc.pp.highly_variable_genes(pbmc, flavor="seurat_v3")


def test_batches():
    adata = pbmc68k_reduced()
    adata[:100, :100].X = np.zeros((100, 100))

    adata.obs["batch"] = ["0" if i < 100 else "1" for i in range(adata.n_obs)]
    adata_1 = adata[adata.obs["batch"] == "0"].copy()
    adata_2 = adata[adata.obs["batch"] == "1"].copy()

    sc.pp.highly_variable_genes(
        adata,
        batch_key="batch",
        flavor="cell_ranger",
        n_top_genes=200,
    )

    sc.pp.filter_genes(adata_1, min_cells=1)
    sc.pp.filter_genes(adata_2, min_cells=1)
    hvg1 = sc.pp.highly_variable_genes(
        adata_1, flavor="cell_ranger", n_top_genes=200, inplace=False
    )
    assert hvg1 is not None
    hvg2 = sc.pp.highly_variable_genes(
        adata_2, flavor="cell_ranger", n_top_genes=200, inplace=False
    )
    assert hvg2 is not None

    np.testing.assert_allclose(
        adata.var["dispersions_norm"].iat[100],
        0.5 * hvg1["dispersions_norm"].iat[0] + 0.5 * hvg2["dispersions_norm"].iat[100],
        rtol=1.0e-7,
        atol=1.0e-7,
    )
    np.testing.assert_allclose(
        adata.var["dispersions_norm"].iat[101],
        0.5 * hvg1["dispersions_norm"].iat[1] + 0.5 * hvg2["dispersions_norm"].iat[101],
        rtol=1.0e-7,
        atol=1.0e-7,
    )
    np.testing.assert_allclose(
        adata.var["dispersions_norm"].iat[0],
        0.5 * hvg2["dispersions_norm"].iat[0],
        rtol=1.0e-7,
        atol=1.0e-7,
    )

    colnames = [
        "means",
        "dispersions",
        "dispersions_norm",
        "highly_variable",
    ]

    assert np.all(np.isin(colnames, hvg1.columns))


def test_degenerate_batches():
    adata = AnnData(
        X=np.random.randn(10, 100),
        obs=dict(batch=pd.Categorical([*([1] * 4), *([2] * 5), 3])),
    )
    sc.pp.highly_variable_genes(adata, batch_key="batch")


@needs.skmisc
def test_seurat_v3_mean_var_output_with_batchkey():
    pbmc = pbmc3k()
    pbmc.var_names_make_unique()
    n_cells = pbmc.shape[0]
    batch = np.zeros((n_cells), dtype=int)
    batch[1500:] = 1
    pbmc.obs["batch"] = batch

    # true_mean, true_var = _get_mean_var(pbmc.X)
    true_mean = np.mean(pbmc.X.toarray(), axis=0)
    true_var = np.var(pbmc.X.toarray(), axis=0, dtype=np.float64, ddof=1)

    result_df = sc.pp.highly_variable_genes(
        pbmc, batch_key="batch", flavor="seurat_v3", n_top_genes=4000, inplace=False
    )
    np.testing.assert_allclose(true_mean, result_df["means"], rtol=2e-05, atol=2e-05)
    np.testing.assert_allclose(true_var, result_df["variances"], rtol=2e-05, atol=2e-05)


def test_cellranger_n_top_genes_warning():
    X = np.random.poisson(2, (100, 30))
    adata = AnnData(X)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    with pytest.warns(
        UserWarning,
        match="`n_top_genes` > number of normalized dispersions, returning all genes with normalized dispersions.",
    ):
        sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="cell_ranger")


def test_cutoff_info():
    adata = pbmc3k()[:200].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    with pytest.warns(UserWarning, match="pass `n_top_genes`, all cutoffs are ignored"):
        sc.pp.highly_variable_genes(adata, n_top_genes=10, max_mean=3.1)


@pytest.mark.parametrize("flavor", ["seurat", "cell_ranger"])
@pytest.mark.parametrize("array_type", ARRAY_TYPES)
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_subset_inplace_consistency(flavor, array_type, batch_key):
    """Tests `n_top_genes=n`.

    - if `inplace` and `subset` interact correctly
    - for both the `seurat` and `cell_ranger` flavors
    - for dask arrays and non-dask arrays
    - for both with and without batch_key
    """
    adata = sc.datasets.blobs(n_observations=20, n_variables=80, random_state=0)
    rng = np.random.default_rng(0)
    adata.obs["batch"] = rng.choice(["a", "b"], adata.shape[0])
    adata.X = array_type(np.abs(adata.X).astype(int))

    if flavor in {"seurat", "cell_ranger"}:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    elif flavor == "seurat_v3":
        pass

    else:
        msg = f"Unknown flavor {flavor}"
        raise ValueError(msg)

    n_genes = adata.shape[1]

    adatas: dict[bool, AnnData] = {}
    dfs: dict[bool, pd.DataFrame] = {}
    # for loops instead of parametrization to compare between settings
    for subset, inplace in itertools.product([True, False], repeat=2):
        adata_copy = adata.copy()

        output_df = sc.pp.highly_variable_genes(
            adata_copy,
            flavor=flavor,
            n_top_genes=15,
            batch_key=batch_key,
            subset=subset,
            inplace=inplace,
        )

        assert (output_df is None) == inplace
        assert len(adata_copy.var if inplace else output_df) == (
            15 if subset else n_genes
        )
        assert sum((adata_copy.var if inplace else output_df)["highly_variable"]) == 15

        if not inplace:
            assert isinstance(output_df, pd.DataFrame)

        if inplace:
            assert subset not in adatas
            adatas[subset] = adata_copy
        else:
            assert subset not in dfs
            dfs[subset] = output_df

    # check that the results are consistent for subset True/False: inplace True
    adata_subset = adatas[False][:, adatas[False].var["highly_variable"]]
    assert adata_subset.var_names.equals(adatas[True].var_names)

    # check that the results are consistent for subset True/False: inplace False
    df_subset = dfs[False][dfs[False]["highly_variable"]]
    assert df_subset.index.equals(dfs[True].index)

    # check that the results are consistent for inplace True/False: subset True
    assert adatas[True].var_names.equals(dfs[True].index)


@pytest.mark.parametrize("flavor", ["seurat", "cell_ranger"])
@pytest.mark.parametrize("batch_key", [None, "batch"], ids=["single", "batched"])
@pytest.mark.parametrize(
    "to_dask", [p for p in ARRAY_TYPES if "dask" in p.values[0].__name__]
)
def test_dask_consistency(adata: AnnData, flavor, batch_key, to_dask):
    adata.X = np.abs(adata.X).astype(int)
    if batch_key is not None:
        adata.obs[batch_key] = np.tile(["a", "b"], adata.shape[0] // 2)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata_dask = adata.copy()
    adata_dask.X = to_dask(adata_dask.X)

    output_mem, output_dask = (
        sc.pp.highly_variable_genes(ad, flavor=flavor, n_top_genes=15, inplace=False)
        for ad in [adata, adata_dask]
    )

    assert isinstance(output_mem, pd.DataFrame)
    assert isinstance(output_dask, pd.DataFrame)

    assert_index_equal(adata.var_names, output_mem.index, check_names=False)
    assert_index_equal(adata.var_names, output_dask.index, check_names=False)

    assert_frame_equal(output_mem, output_dask, atol=1e-4)
