"""Tests for `flavor="poisson_gene_selection"` in `pp.highly_variable_genes`."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import scanpy as sc
from testing.scanpy._helpers import _check_check_values_warnings
from testing.scanpy._pytest.params import ARRAY_TYPES

STAT_COLUMNS = [
    "observed_fraction_zeros",
    "expected_fraction_zeros",
    "prob_zero_enrichment",
    "prob_zero_enrichment_rank",
]
FLAVOR = "poisson_gene_selection"


def _counts(
    seed: int, n_cells: int, n_genes: int, *, zero_inflate: float = 0.0
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    lam = rng.uniform(0.05, 2.0, size=n_genes)
    x = rng.poisson(lam, size=(n_cells, n_genes)).astype(np.float32)
    if zero_inflate:
        x[rng.random(x.shape) < zero_inflate] = 0
    return x


def _adata(x, *, batches: list[str] | None = None) -> AnnData:
    obs = pd.DataFrame(index=[f"c{i}" for i in range(x.shape[0])])
    if batches is not None:
        obs["batch"] = pd.Categorical(batches)
    var = pd.DataFrame(index=[f"g{i}" for i in range(x.shape[1])])
    return AnnData(x.copy(), obs=obs, var=var)


def _reference_stats(x: np.ndarray) -> dict[str, np.ndarray]:
    """Independent dense NumPy implementation of the three formulas."""
    x = x.astype(np.float64)
    lib = x.sum(axis=1)
    p = x.sum(axis=0) / x.sum()
    obs = 1 - (x > 0).sum(axis=0) / x.shape[0]
    exp_ = np.exp(-p[:, None] * lib[None, :]).mean(axis=1)
    return dict(obs=obs, exp=exp_, score=obs * (1 - exp_))


def _reference_combine(
    scores: list[np.ndarray], n_top_genes: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Port of the batch combination in rapids-singlecell / scvi-tools.

    Returns (median rank, n batches, highly_variable).
    """
    s = np.vstack(scores)
    ranks = s.argsort(axis=1).argsort(axis=1)
    df = pd.DataFrame({
        "prob_zero_enriched_nbatches": (ranks >= s.shape[1] - n_top_genes).sum(0),
        "prob_zero_enrichment_rank": np.median(ranks, axis=0),
    })
    top = df.nlargest(n_top_genes, list(df.columns)).index
    hv = np.zeros(s.shape[1], dtype=bool)
    hv[top] = True
    return df.iloc[:, 1].to_numpy(), df.iloc[:, 0].to_numpy(), hv


def test_manual_math():
    # 4 cells x 2 genes, every library size is 10.
    # gene a: [0, 0, 0, 10]  -> p_a = 10/40, p_exp = exp(-2.5), p_obs = 3/4
    # gene b: [10, 10, 10, 0] -> p_b = 30/40, p_exp = exp(-7.5), p_obs = 1/4
    x = np.array([[0, 10], [0, 10], [0, 10], [10, 0]], dtype=np.float64)
    adata = AnnData(x, var=pd.DataFrame(index=["a", "b"]))

    sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=1)

    v = adata.var
    np.testing.assert_allclose(v["observed_fraction_zeros"], [0.75, 0.25])
    np.testing.assert_allclose(
        v["expected_fraction_zeros"], [np.exp(-2.5), np.exp(-7.5)], rtol=1e-12
    )
    np.testing.assert_allclose(
        v["prob_zero_enrichment"],
        [0.75 * (1 - np.exp(-2.5)), 0.25 * (1 - np.exp(-7.5))],
        rtol=1e-12,
    )
    # 0-based, higher = more enriched (scvi-tools / rapids-singlecell)
    assert v["prob_zero_enrichment_rank"].tolist() == [1.0, 0.0]
    assert v["highly_variable"].tolist() == [True, False]
    assert "prob_zero_enriched_nbatches" not in v
    assert adata.uns["hvg"] == {"flavor": FLAVOR}


def test_all_zero_gene_ranks_last():
    x = np.array([[0, 10, 0], [0, 10, 0], [0, 10, 0], [10, 0, 0]], dtype=np.float64)
    adata = AnnData(x, var=pd.DataFrame(index=["a", "b", "zero"]))

    sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=2)

    v = adata.var
    # an all-zero gene leaves the other genes' statistics unchanged
    np.testing.assert_allclose(
        v.loc["a", "expected_fraction_zeros"], np.exp(-2.5), rtol=1e-12
    )
    # p_g = 0 -> expected zeros = observed zeros = 1 -> score 0
    assert v.loc["zero", "expected_fraction_zeros"] == 1.0
    assert v.loc["zero", "prob_zero_enrichment"] == 0.0
    assert v.loc["zero", "prob_zero_enrichment_rank"] == 0.0
    assert v["highly_variable"].tolist() == [True, True, False]


@pytest.mark.parametrize(
    "array_type",
    [
        p
        for p in ARRAY_TYPES
        if "dask" not in p.id or ("1d_chunked" in p.id and "csc" not in p.id)
    ],
)
@pytest.mark.parametrize("batch_key", [None, "batch"])
def test_array_types_match_reference(array_type, batch_key):
    x = _counts(36, 200, 40, zero_inflate=0.5)
    batches = ["a"] * 90 + ["b"] * 110
    adata = _adata(x, batches=batches)
    adata.X = array_type(adata.X)

    sc.pp.highly_variable_genes(
        adata, flavor=FLAVOR, n_top_genes=10, batch_key=batch_key
    )
    v = adata.var

    parts = [x] if batch_key is None else [x[:90], x[90:]]
    refs = [_reference_stats(xi) for xi in parts]
    # The references order exact ties via NumPy's default (unstable) argsort,
    # which is implementation-defined; this seed has no ties in any batch.
    assert all(np.unique(r["score"]).size == x.shape[1] for r in refs)
    for col, key in [
        ("observed_fraction_zeros", "obs"),
        ("expected_fraction_zeros", "exp"),
        ("prob_zero_enrichment", "score"),
    ]:
        expected = np.median(np.vstack([r[key] for r in refs]), axis=0)
        np.testing.assert_allclose(v[col], expected, rtol=1e-10)

    med_rank, nbatches, hv = _reference_combine([r["score"] for r in refs], 10)
    np.testing.assert_array_equal(v["prob_zero_enrichment_rank"], med_rank)
    np.testing.assert_array_equal(v["highly_variable"], hv)
    if batch_key is not None:
        np.testing.assert_array_equal(v["prob_zero_enriched_nbatches"], nbatches)


@pytest.mark.parametrize(
    "array_type",
    [p for p in ARRAY_TYPES if "dask" in p.id and "1d_chunked" not in p.id],
)
def test_dask_feature_axis_chunking_raises(array_type):
    adata = _adata(_counts(0, 50, 20))
    adata.X = array_type(adata.X)
    with pytest.raises(ValueError, match=r"chunking along the first axis"):
        sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=5)


def test_matches_scvi_monte_carlo():
    """The score is the limit of scvi-tools' Monte-Carlo estimate.

    scvi-tools draws `Binomial(probs=p)` with the default `total_count=1`,
    i.e. Bernoulli, and averages `obs_draw > exp_draw` over samples.
    """
    adata = _adata(_counts(7, 150, 30, zero_inflate=0.4))
    df = sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=5, inplace=False)
    p_obs = df["observed_fraction_zeros"].to_numpy()
    p_exp = df["expected_fraction_zeros"].to_numpy()

    n_samples = 20_000
    rng = np.random.default_rng(0)
    mc = (
        (rng.random((n_samples, p_obs.size)) < p_obs)
        > (rng.random((n_samples, p_exp.size)) < p_exp)
    ).mean(axis=0)
    # 5 standard errors of a proportion at worst case p = 0.5
    np.testing.assert_allclose(
        df["prob_zero_enrichment"], mc, atol=5 * 0.5 / np.sqrt(n_samples)
    )


def test_ties_favor_earlier_gene():
    # genes 0 and 2 are identical, so their scores tie exactly
    x = _counts(4, 50, 5, zero_inflate=0.3)
    x[:, 2] = x[:, 0]
    adata = _adata(x)
    sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=2)
    v = adata.var
    assert v["prob_zero_enrichment"].iloc[0] == v["prob_zero_enrichment"].iloc[2]
    assert (
        v["prob_zero_enrichment_rank"].iloc[0]
        == v["prob_zero_enrichment_rank"].iloc[2] + 1
    )


def test_sparse_explicit_zeros_are_not_expressing():
    # same matrix as test_manual_math, plus an explicitly stored zero at (0, 0)
    x = sp.coo_array(
        (
            np.array([0.0, 10.0, 10.0, 10.0, 10.0]),
            (np.array([0, 0, 1, 2, 3]), np.array([0, 1, 1, 1, 0])),
        ),
        shape=(4, 2),
    ).tocsr()
    assert x.nnz == 5
    adata = AnnData(x, var=pd.DataFrame(index=["a", "b"]))

    sc.pp.highly_variable_genes(adata, flavor=FLAVOR, n_top_genes=1)

    np.testing.assert_allclose(adata.var["observed_fraction_zeros"], [0.75, 0.25])
    assert adata.var["prob_zero_enrichment_rank"].tolist() == [1.0, 0.0]


def test_identical_batches_match_single_batch():
    # duplicating cells into two batches leaves every statistic unchanged
    x = _counts(1, 60, 30, zero_inflate=0.4)
    single = _adata(x)
    sc.pp.highly_variable_genes(single, flavor=FLAVOR, n_top_genes=8)

    double = _adata(np.vstack([x, x]), batches=["a"] * 60 + ["b"] * 60)
    sc.pp.highly_variable_genes(double, flavor=FLAVOR, n_top_genes=8, batch_key="batch")

    for col in STAT_COLUMNS:
        np.testing.assert_allclose(double.var[col], single.var[col], rtol=1e-10)
    np.testing.assert_array_equal(
        double.var["highly_variable"], single.var["highly_variable"]
    )
    hv = double.var["highly_variable"]
    assert (double.var.loc[hv, "prob_zero_enriched_nbatches"] == 2).all()


def test_batch_key_validation():
    adata = _adata(_counts(2, 20, 10))
    adata.obs["batch"] = pd.Categorical(["a"] * 10 + [None] * 10)
    with pytest.raises(ValueError, match=r"missing batch labels"):
        sc.pp.highly_variable_genes(
            adata, flavor=FLAVOR, n_top_genes=5, batch_key="batch"
        )


def test_empty_batch_categories_are_skipped():
    x = _counts(3, 40, 15, zero_inflate=0.3)
    plain = _adata(x, batches=["a"] * 20 + ["b"] * 20)
    extra = _adata(x)
    extra.obs["batch"] = pd.Categorical(
        ["a"] * 20 + ["b"] * 20, categories=["a", "unused", "b"]
    )
    kw = dict(flavor=FLAVOR, n_top_genes=5, batch_key="batch", inplace=False)
    pd.testing.assert_frame_equal(
        sc.pp.highly_variable_genes(plain, **kw),
        sc.pp.highly_variable_genes(extra, **kw),
    )


def test_check_values():
    x = np.array([[0.0, 1.0], [2.5, 0.0], [0.0, 3.0], [1.0, 2.0]])
    _check_check_values_warnings(
        function=sc.pp.highly_variable_genes,
        adata=_adata(x),
        expected_warning=(
            "`flavor='poisson_gene_selection'` expects raw count data, "
            "but non-integers were found."
        ),
        kwargs=dict(flavor=FLAVOR, n_top_genes=1),
    )


def test_n_top_genes():
    adata = _adata(_counts(6, 100, 2500, zero_inflate=0.3))
    sc.pp.highly_variable_genes(adata, flavor=FLAVOR)  # default
    assert adata.var["highly_variable"].sum() == 2000

    small = _adata(_counts(6, 30, 4))
    sc.pp.highly_variable_genes(small, flavor=FLAVOR, n_top_genes=10)
    assert small.var["highly_variable"].all()

    with pytest.raises(ValueError, match=r"positive integer"):
        sc.pp.highly_variable_genes(small, flavor=FLAVOR, n_top_genes=0)


def test_layer():
    x = _counts(5, 80, 30, zero_inflate=0.3)
    layered = _adata(np.log1p(x))
    layered.layers["counts"] = x.copy()
    sc.pp.highly_variable_genes(layered, flavor=FLAVOR, n_top_genes=6, layer="counts")
    direct = _adata(x)
    sc.pp.highly_variable_genes(direct, flavor=FLAVOR, n_top_genes=6)

    pd.testing.assert_frame_equal(layered.var, direct.var)
    np.testing.assert_allclose(layered.X, np.log1p(x))  # X untouched
