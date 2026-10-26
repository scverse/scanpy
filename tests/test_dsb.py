from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from scanpy.preprocessing._dsb import _mixture_fits, _threshold
from testing.scanpy._helpers import check_rep_mutation, check_rep_results
from testing.scanpy._pytest.params import ARRAY_TYPES

pytestmark = pytest.mark.filterwarnings("ignore:Proteins with low background")

# Reference values from the dsb R package, see tests/_scripts/dsb_reference.R
DATA = Path(__file__).parent / "_data" / "dsb"
ISOTYPES = [
    "MouseIgG1kappaisotype_PROT",
    "MouseIgG2akappaisotype_PROT",
    "Mouse IgG2bkIsotype_PROT",
    "RatIgG2bkIsotype_PROT",
]
F32_TOL = 1e-5  # float32 output of float64 computations


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / f"{name}.csv.gz")


@pytest.fixture(scope="module")
def counts() -> tuple[AnnData, AnnData]:
    cells, empty = _read("cells"), _read("empty")

    def to_adata(df, prefix):
        obs = pd.DataFrame(index=[f"{prefix}{i}" for i in range(len(df))])
        return AnnData(
            df.to_numpy(dtype=np.float64), obs=obs, var=pd.DataFrame(index=df.columns)
        )

    return to_adata(cells, "c"), to_adata(empty, "e")


@pytest.fixture(scope="module")
def expected() -> np.ndarray:
    return _read("default").to_numpy()


def _tol(dtype):
    return 1e-10 if dtype == np.float64 else F32_TOL


CONFIGS = {
    "default": ({"isotype_controls": ISOTYPES}, True),
    "meansub_clip": (
        {
            "isotype_controls": ISOTYPES,
            "pseudocount": 5,
            "scale_factor": "mean_subtract",
            "quantile_clip": (0.001, 0.9995),
        },
        True,
    ),
    "negative": ({"isotype_controls": ISOTYPES}, False),
    "no_isotypes": ({"use_isotype_controls": False}, True),
    "negative_no_isotypes": ({"use_isotype_controls": False}, False),
}


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
@pytest.mark.parametrize("config", list(CONFIGS))
def test_matches_r(counts, config, dtype):
    if not (DATA / f"{config}.csv.gz").exists():
        pytest.skip("regenerate the references with tests/_scripts/dsb_reference.R")
    cells, empty = counts
    kwargs, with_empty = CONFIGS[config]
    res = sc.pp.dsb(
        cells, empty if with_empty else None, dtype=dtype, copy=True, **kwargs
    )
    assert res.X.dtype == dtype
    np.testing.assert_allclose(
        res.X, _read(config).to_numpy(), rtol=0, atol=_tol(dtype)
    )


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_stats(counts, dtype):
    cells, empty = counts
    res = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES, dtype=dtype, copy=True)
    ref = _read("default_technical_stats")
    tol = _tol(dtype)
    np.testing.assert_allclose(
        res.obs["dsb_cell_background_mean"],
        ref["cellwise_background_mean"],
        rtol=0,
        atol=tol,
    )
    # the sign of a principal component is arbitrary in R
    tech = res.obs["dsb_technical_component"].to_numpy()
    ref_tech = ref["dsb_technical_component"].to_numpy()
    sign = np.sign(tech @ ref_tech)
    np.testing.assert_allclose(tech, sign * ref_tech, rtol=0, atol=tol)
    assert np.corrcoef(tech, res.obs["dsb_cell_background_mean"])[0, 1] > 0

    log_empty = np.log(empty.X + 10)
    np.testing.assert_allclose(
        res.var["dsb_ambient_mean"], log_empty.mean(0), rtol=1e-6
    )
    np.testing.assert_allclose(
        res.var["dsb_ambient_sd"], log_empty.std(0, ddof=1), rtol=1e-6
    )


def test_no_denoising(counts):
    cells, empty = counts
    res = sc.pp.dsb(cells, empty, denoise_counts=False, copy=True)
    log_empty = np.log(empty.X + 10)
    step1 = (np.log(cells.X + 10) - log_empty.mean(0)) / log_empty.std(0, ddof=1)
    np.testing.assert_allclose(res.X, step1, rtol=1e-12, atol=1e-12)
    assert "dsb_technical_component" not in res.obs


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def test_input_types(counts, expected, array_type):
    cells, empty = counts
    a = AnnData(array_type(cells.X), obs=cells.obs, var=cells.var)
    e = AnnData(array_type(empty.X), obs=empty.obs, var=empty.var)
    sc.pp.dsb(a, e, isotype_controls=ISOTYPES)
    assert isinstance(a.X, np.ndarray)
    np.testing.assert_allclose(a.X, expected, rtol=0, atol=1e-10)


@pytest.mark.parametrize("fmt", ["dense", "csr", "csc"])
def test_backed_empty(counts, expected, tmp_path, fmt):
    cells, empty = counts
    x = empty.X if fmt == "dense" else getattr(sparse, f"{fmt}_matrix")(empty.X)
    AnnData(x, obs=empty.obs, var=empty.var).write_h5ad(tmp_path / "empty.h5ad")
    backed = ad.read_h5ad(tmp_path / "empty.h5ad", backed="r")
    res = sc.pp.dsb(cells, backed, isotype_controls=ISOTYPES, copy=True)
    np.testing.assert_allclose(res.X, expected, rtol=0, atol=1e-10)


@pytest.mark.parametrize("negative", [False, True], ids=["empty", "negative"])
def test_layer(counts, negative):
    cells, empty = counts
    iso = list(cells.var_names.get_indexer(ISOTYPES).astype(str))
    e = AnnData(empty.X.copy(), layers={"layer": empty.X.copy()})

    def func(adata, **kwargs):
        return sc.pp.dsb(adata, None if negative else e, isotype_controls=iso, **kwargs)

    check_rep_mutation(func, cells.X, fields=("layer",))
    check_rep_results(func, cells.X, fields=("layer",))


def test_key_added(counts, expected):
    cells, empty = counts
    a = cells.copy()
    sc.pp.dsb(a, empty, isotype_controls=ISOTYPES, key_added="dsb")
    np.testing.assert_allclose(a.layers["dsb"], expected, rtol=0, atol=1e-10)
    np.testing.assert_array_equal(a.X, cells.X)


def test_copy_leaves_input(counts):
    cells, empty = counts
    before = cells.X.copy()
    res = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES, copy=True)
    assert res is not cells
    np.testing.assert_array_equal(cells.X, before)
    assert "dsb_ambient_mean" not in cells.var


def test_view_leaves_parent(counts, expected):
    cells, empty = counts
    parent = cells.copy()
    view = parent[:50]
    with pytest.warns(UserWarning, match="Received a view"):
        sc.pp.dsb(view, empty, isotype_controls=ISOTYPES)
    np.testing.assert_array_equal(parent.X, cells.X)
    assert "dsb_ambient_mean" not in parent.var
    ref = sc.pp.dsb(cells[:50].copy(), empty, isotype_controls=ISOTYPES, copy=True)
    np.testing.assert_array_equal(view.X, ref.X)


def test_stats_are_reset(counts):
    cells, empty = counts
    a = cells.copy()
    sc.pp.dsb(a, empty, isotype_controls=ISOTYPES, key_added="with_empty")
    sc.pp.dsb(a, isotype_controls=ISOTYPES, key_added="negative")
    assert "dsb_ambient_sd" not in a.var
    sc.pp.dsb(a, denoise_counts=False, key_added="step1")
    assert not a.obs.columns.str.startswith("dsb_").any()


def test_empty_protein_order(counts, expected):
    cells, empty = counts
    shuffled = empty[:, empty.var_names[::-1]].copy()
    with pytest.warns(UserWarning, match="Reordered"):
        res = sc.pp.dsb(cells, shuffled, isotype_controls=ISOTYPES, copy=True)
    np.testing.assert_allclose(res.X, expected, rtol=0, atol=1e-10)
    with pytest.raises(ValueError, match="different proteins"):
        sc.pp.dsb(cells, empty[:, 1:].copy(), isotype_controls=ISOTYPES, copy=True)


@pytest.mark.filterwarnings("ignore:Variable names are not unique")
def test_duplicate_proteins(counts):
    cells, empty = counts
    names = cells.var_names.to_numpy().copy()
    names[1] = names[0]
    dup = AnnData(empty.X, var=pd.DataFrame(index=names))
    with pytest.raises(ValueError, match="duplicated"):
        sc.pp.dsb(cells, dup, isotype_controls=ISOTYPES, copy=True)
    with pytest.raises(ValueError, match="not unique"):
        sc.pp.dsb(
            AnnData(cells.X, var=pd.DataFrame(index=names)),
            isotype_controls=ISOTYPES,
            copy=True,
        )


@pytest.mark.parametrize(
    ("kwargs", "with_empty", "match"),
    [
        pytest.param({}, True, "Potential isotype controls", id="no_isotypes"),
        pytest.param(
            {"isotype_controls": []}, True, "Potential isotype", id="empty_isotypes"
        ),
        pytest.param(
            {"isotype_controls": ["nope"]}, True, "not found", id="missing_isotype"
        ),
        pytest.param(
            {"isotype_controls": ["nope"], "denoise_counts": False},
            True,
            "not found",
            id="missing_isotype_no_denoising",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "use_isotype_controls": False},
            True,
            "not both",
            id="isotypes_unused",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "scale_factor": "scale"},
            True,
            "scale_factor",
            id="scale_factor",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "scale_factor": "mean_subtract"},
            False,
            "only applies with `adata_empty`",
            id="scale_factor_negative",
        ),
        pytest.param(
            {"fast_km": True, "denoise_counts": False},
            True,
            "fast_km",
            id="fast_km_no_denoising",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "pseudocount": 0},
            True,
            "pseudocount",
            id="pseudocount",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "quantile_clip": (0.9, 0.1)},
            True,
            "quantile_clip",
            id="quantile_clip",
        ),
        pytest.param(
            {"isotype_controls": ISOTYPES, "dtype": np.int32},
            True,
            "dtype",
            id="dtype",
        ),
    ],
)
def test_invalid_arguments(counts, kwargs, with_empty, match):
    cells, empty = counts
    with pytest.raises(ValueError, match=match):
        sc.pp.dsb(cells, empty if with_empty else None, copy=True, **kwargs)


def test_isotype_inputs(counts):
    cells, empty = counts
    as_list = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES[:1], copy=True).X
    as_str = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES[0], copy=True).X
    np.testing.assert_array_equal(as_str, as_list)
    full = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES, copy=True).X
    gen = (p for p in cells.var_names if p in ISOTYPES)
    np.testing.assert_array_equal(
        sc.pp.dsb(cells, empty, isotype_controls=gen, copy=True).X, full
    )


@pytest.mark.parametrize(
    "case", ["negative", "nan", "x_none", "missing_layer", "one_empty"]
)
def test_invalid_counts(counts, case):
    cells, empty = counts
    a, e, kwargs = cells.copy(), empty, {}
    if case == "negative":
        a.X[0, 0], match = -1, "raw counts"
    elif case == "nan":
        a.X[0, 0], match = np.nan, "raw counts"
    elif case == "x_none":
        a = AnnData(obs=cells.obs, var=cells.var, layers={"counts": cells.X})
        match = r"`adata.X` is None"
    elif case == "missing_layer":
        a.layers["counts"] = a.X
        kwargs, match = {"layer": "counts"}, "`adata_empty` has no layer"
    else:
        e, match = empty[:1], "at least 2 empty droplets"
    with pytest.raises(ValueError, match=match):
        sc.pp.dsb(a, e, isotype_controls=ISOTYPES, copy=True, **kwargs)


def test_constant_ambient_protein(counts):
    cells, empty = counts
    e = empty.copy()
    e.X[:, 0] = 0
    match = f"same count in every empty droplet.*{re.escape(e.var_names[0])}"
    with pytest.raises(ValueError, match=match):
        sc.pp.dsb(cells, e, isotype_controls=ISOTYPES, copy=True)
    # centering alone is well defined
    sc.pp.dsb(
        cells, e, isotype_controls=ISOTYPES, scale_factor="mean_subtract", copy=True
    )


def test_constant_isotype(counts):
    cells, empty = counts
    a = cells.copy()
    a.X[:, a.var_names.get_loc(ISOTYPES[0])] = 0
    with pytest.raises(
        ValueError, match=re.escape(f"constant across cells: {ISOTYPES[:1]}")
    ):
        sc.pp.dsb(a, empty, isotype_controls=ISOTYPES, copy=True)


def test_identical_cells(counts):
    # a constant covariate leaves the data unchanged; with isotypes, constant
    # columns cannot be standardized
    cells, empty = counts
    same = AnnData(np.repeat(cells.X[:1], 20, axis=0), var=cells.var)
    step1 = sc.pp.dsb(same, empty, denoise_counts=False, copy=True).X
    with pytest.warns(UserWarning, match="technical component is constant"):
        res = sc.pp.dsb(same, empty, use_isotype_controls=False, copy=True)
    np.testing.assert_array_equal(res.X, step1)
    with pytest.raises(ValueError, match="constant across cells"):
        sc.pp.dsb(same, empty, isotype_controls=ISOTYPES, copy=True)


def test_fast_km_is_optimal_two_means(counts):
    cells, empty = counts
    res = sc.pp.dsb(cells, empty, isotype_controls=ISOTYPES, fast_km=True, copy=True)
    step1 = sc.pp.dsb(cells, empty, denoise_counts=False, copy=True).X
    for x, low in zip(step1, res.obs["dsb_cell_background_mean"], strict=True):
        xs = np.sort(x)
        sse = [
            ((xs[:s] - xs[:s].mean()) ** 2).sum()
            + ((xs[s:] - xs[s:].mean()) ** 2).sum()
            for s in range(1, len(xs))
        ]
        best = int(np.argmin(sse)) + 1
        assert low == pytest.approx(xs[:best].mean(), abs=1e-12)


def test_cell_fit_fallback():
    # every protein has the same ambient distribution, so a cell with only two
    # distinct counts has two distinct normalized values: no mixture fits
    gen = np.random.default_rng(0)
    n_proteins = 6
    empty = np.repeat(gen.poisson(3, (200, 1)), n_proteins, axis=1).astype(float)
    cells = gen.poisson(gen.uniform(1, 40, n_proteins), (30, n_proteins)).astype(float)
    cells[0] = [0, 0, 0, 20, 20, 20]
    with pytest.warns(UserWarning, match=r"could not be fit for \d+ cells"):
        res = sc.pp.dsb(
            AnnData(cells), AnnData(empty), use_isotype_controls=False, copy=True
        )
    step1 = sc.pp.dsb(AnnData(cells), AnnData(empty), denoise_counts=False, copy=True)
    assert res.obs["dsb_cell_background_mean"].iloc[0] == step1.X[0].min()


def test_protein_fit_fallback(counts):
    cells, _ = counts
    a = cells.copy()
    a.X[:, 0] = np.r_[np.zeros(60), np.full(40, 5.0)]  # two values: no fits
    name = a.var_names[0]
    with pytest.warns(
        UserWarning, match=f"could not be fit for proteins .*{re.escape(name)}"
    ):
        res = sc.pp.dsb(a, denoise_counts=False, copy=True)
    assert res.var["dsb_ambient_mean"].iloc[0] == 0
    np.testing.assert_array_equal(res.X[:, 0], np.log(a.X[:, 0] + 1))


def test_negative_model_rng(counts):
    # more than 2000 cells: the per-protein fits start from random subsets
    cells, _ = counts
    gen = np.random.default_rng(0)
    big = AnnData(cells.X[gen.integers(0, cells.n_obs, 2500)], var=cells.var)
    res = [
        sc.pp.dsb(big, isotype_controls=ISOTYPES, rng=0, copy=True).X for _ in range(2)
    ]
    np.testing.assert_array_equal(*res)


# --- mixture fits ----------------------------------------------------------


def _ref_threshold(x):
    """Type 7 quantile at the simplest j / k with min < q < max."""
    xs, n = np.sort(x), len(x)
    for k in range(2, 2 * n):
        for j in range(1, k):
            lo, rem = divmod((n - 1) * j, k)
            q = xs[lo] + rem / k * (xs[min(lo + 1, n - 1)] - xs[lo])
            if xs[0] < q < xs[-1]:
                return q
    return np.nan


def _ref_mixture(x, sub=None):
    """NumPy reference of the fits (see ``_mixture_fits``).

    `sub` is the random subset that starts long vectors.
    """
    eps = np.finfo(np.float64).eps

    def mstep(v, z, equal):
        nk = z.sum(0)
        mu = (z * v[:, None]).sum(0) / nk
        ss = (z * (v[:, None] - mu) ** 2).sum(0)
        return mu, np.full(2, ss.sum() / len(v)) if equal else ss / nk, nk / len(v)

    def estep(mu, s2, pro):
        lg = (
            np.log(pro)
            - 0.5 * np.log(2 * np.pi * s2)
            - 0.5 * (x[:, None] - mu) ** 2 / s2
        )
        ll = np.logaddexp(lg[:, 0], lg[:, 1])
        return np.exp(lg - ll[:, None]), ll.sum()

    def valid(p):
        return np.all(np.isfinite(p[0])) and np.all(p[1] > eps)

    start = x if sub is None else sub
    c = start < _ref_threshold(start)
    fits = []
    for equal, npar in ((True, 4), (False, 5)):
        p = mstep(start, np.stack([c, ~c], 1).astype(float), equal)
        skip, prev, ok = sub is not None, None, valid(p)
        while ok:
            z, ll = estep(*p)
            pn = mstep(x, z, equal)
            if skip:
                skip = False
            elif prev is not None and abs(ll - prev) < 1e-5 * (1 + abs(ll)):
                break
            else:
                prev = ll
            ok = valid(pn)
            p = pn if ok else p
        if not ok:
            continue
        if sub is not None or ((p[2] - pn[2]) ** 2).sum() > np.sqrt(eps):
            p = pn
        fits.append((2 * ll - npar * np.log(len(x)), p))
    return max(fits, key=lambda f: f[0])[1] if fits else None


def _random_rows(rng, n_rows):
    rows = []
    for i in range(n_rows):
        n = int(rng.integers(20, 150))
        kind = i % 3
        if kind == 0:
            x = np.r_[rng.normal(0, 1, n // 2), rng.normal(rng.uniform(1, 6), 2, n)]
        elif kind == 1:  # log counts with many ties
            x = np.log(rng.poisson(rng.uniform(0.05, 3), 2 * n) + 1.0)
        else:
            x = np.log(rng.poisson(rng.uniform(0.05, 30), 2 * n) + 10.0)
        rows.append(x[: 2 * n])
    return rows


def test_mixture_fits_match_reference():
    rng = np.random.default_rng(0)
    for x in _random_rows(rng, 60):
        assert _threshold(np.sort(x)) == pytest.approx(_ref_threshold(x), abs=1e-12)
        got = _mixture_fits(x[None, :], rng).params[0]
        ref = _ref_mixture(x)
        if ref is None:
            assert np.isnan(got).all()
            continue
        np.testing.assert_allclose(got[:2], ref[0], rtol=1e-9, atol=1e-10)
        np.testing.assert_allclose(got[2:4], ref[1], rtol=1e-9, atol=1e-10)
        np.testing.assert_allclose(got[4:], ref[2], rtol=0, atol=1e-10)


def _two_component(n=3000):
    gen = np.random.default_rng(1)
    return np.r_[gen.normal(0, 1, 2000), gen.normal(4, 1.5, n - 2000)]


def test_mixture_fits_subset_path():
    # vectors longer than 2000 start from a random subset of 2000 values
    n = 3000
    gen = np.random.default_rng(1)
    mat = np.stack([
        np.log(gen.poisson(gen.uniform(0.1, 5), n) + 1.0),
        _two_component(n),
    ])
    fit = _mixture_fits(mat, np.random.default_rng(7))
    again = _mixture_fits(mat, np.random.default_rng(7))
    np.testing.assert_array_equal(fit.params, again.params)

    sub_rng = np.random.default_rng(7)
    for x, got in zip(mat, fit.params, strict=True):
        ref = _ref_mixture(x, x[sub_rng.choice(n, 2000, replace=False)])
        np.testing.assert_allclose(got[:2], ref[0], rtol=1e-9, atol=1e-10)


def test_mixture_fits_terminate_on_ties():
    # log(counts + 10) with mostly zero counts: heavy ties, near-singular fits
    rng = np.random.default_rng(3)
    counts = rng.poisson(0.05, (500, 87))
    counts[:, :3] = rng.poisson(50, (500, 3))
    mat = np.log(counts + 10.0)
    # per-row fits, and whole-column fits of the same data (subset path)
    for m in (mat, np.tile(mat, (5, 1)).T):
        fit = _mixture_fits(m, rng)
        ok = ~np.isnan(fit.params[:, 0])
        assert ok.any()
        assert np.isfinite(fit.params[ok]).all()
        assert not fit.capped.any()


def test_mixture_constant_rows():
    mat = np.r_[np.ones((2, 30)), np.arange(60.0).reshape(2, 30)]
    fit = _mixture_fits(mat, np.random.default_rng(0))
    np.testing.assert_array_equal(
        np.isnan(fit.params[:, 0]), [True, True, False, False]
    )


@pytest.mark.parametrize("n", [100, 3000], ids=["full", "subset"])
def test_iteration_cap_is_reported(n):
    mat = np.stack([np.random.default_rng(0).normal(0, 1, n), _two_component(3000)[:n]])
    with pytest.warns(UserWarning, match="iteration limit"):
        fit = _mixture_fits(mat, np.random.default_rng(0), max_iter=2)
    assert fit.capped.all()
    assert not _mixture_fits(mat, np.random.default_rng(0)).capped.any()
