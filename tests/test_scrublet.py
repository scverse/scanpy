from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, concat
from anndata.tests.helpers import assert_equal
from numpy.testing import assert_allclose, assert_array_equal
from scipy import sparse

import scanpy as sc
from testing.scanpy._pytest.marks import needs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    from scanpy._compat import CSBase

pytestmark = [
    needs.skimage,
    # most tests use the deprecated `batch_key`, so they also run without `anndata.acc`
    pytest.mark.filterwarnings(
        "ignore:.*argument batch_key is deprecated:FutureWarning"
    ),
]


def pbmc200() -> AnnData:
    from testing.scanpy._helpers.data import _pbmc3k

    return _pbmc3k()[200:400].copy()


def paul500() -> AnnData:
    from testing.scanpy._helpers.data import _paul15

    return _paul15()[:500].copy()


@pytest.mark.parametrize(
    ("mk_data", "expected_idx", "expected_scores"),
    [
        pytest.param(pbmc200, [13, 138], [0.149254] * 2, id="sparse"),
        pytest.param(paul500, [180], [0.219178], id="dense"),
    ],
)
@pytest.mark.parametrize("use_approx_neighbors", [True, False, None])
def test_scrublet(
    *,
    mk_data: Callable[[], AnnData],
    expected_idx: list[int],
    expected_scores: list[float],
    use_approx_neighbors: bool | None,
):
    """Check that scrublet runs and detects some doublets."""
    adata = mk_data()
    sc.pp.scrublet(adata, use_approx_neighbors=use_approx_neighbors)

    doublet_idx = np.flatnonzero(adata.obs["predicted_doublet"]).tolist()
    assert doublet_idx == expected_idx
    assert_allclose(
        adata.obs["doublet_score"].iloc[doublet_idx],
        expected_scores,
        atol=1e-5,
        rtol=1e-5,
    )


def test_scrublet_batched():
    """Test that Scrublet run works with batched data."""
    adata = pbmc200()
    adata.obs["batch"] = 100 * ["a"] + 100 * ["b"]
    split = [adata[adata.obs["batch"] == x].copy() for x in ("a", "b")]

    sc.pp.scrublet(adata, use_approx_neighbors=False, batch_key="batch")

    doublet_idx = np.flatnonzero(adata.obs["predicted_doublet"]).tolist()
    # only one in the first batch (<100)
    assert doublet_idx == [0, 2, 8, 15, 43, 88, 108, 113, 115, 132, 135, 175]
    assert_allclose(
        adata.obs["doublet_score"].iloc[doublet_idx],
        np.array([0.109375, 0.164835])[([0] * 4 + [1] + [0] * 3 + [1] + [0] * 3)],
        atol=1e-5,
        rtol=1e-5,
    )
    assert adata.uns["scrublet"]["batches"].keys() == {"a", "b"}

    # Check that results are independent
    for s in split:
        sc.pp.scrublet(s, use_approx_neighbors=False)
    merged = concat(split)

    pd.testing.assert_frame_equal(adata.obs[merged.obs.columns], merged.obs)


@pytest.mark.filterwarnings("ignore:Observation names are not unique:UserWarning")
@pytest.mark.parametrize("batch_key", [None, "batch"], ids=["unbatched", "batched"])
def test_scrublet_duplicate_obs_names(batch_key: str | None):
    """Duplicate `obs_names` must not break mapping results back (see #4098)."""
    adata = pbmc200()
    adata.obs["batch"] = 100 * ["a"] + 100 * ["b"]
    # the batches now consist of different cells sharing the same names
    obs_names = [*adata.obs_names[:100]] * 2
    adata.obs_names = obs_names
    assert not adata.obs_names.is_unique

    sc.pp.scrublet(adata, use_approx_neighbors=False, batch_key=batch_key, rng=0)

    assert_array_equal(adata.obs_names, obs_names)  # left untouched
    assert adata.obs["doublet_score"].notna().all()
    assert adata.obs["predicted_doublet"].dtype == bool
    if batch_key is not None:
        # check that the rows line up, using the sub-generators the batched path spawns
        expected = np.concatenate([
            sc.pp.scrublet(
                adata[adata.obs["batch"] == batch].copy(),
                use_approx_neighbors=False,
                rng=sub_rng,
                out=True,
            ).doublet_score
            for batch, sub_rng in zip(
                "ab", np.random.default_rng(0).spawn(2), strict=True
            )
        ])
        assert_allclose(adata.obs["doublet_score"], expected)


@pytest.mark.parametrize("batch_key", [None, "batch"], ids=["unbatched", "batched"])
def test_scrublet_filtered_cells(batch_key: str | None):
    """Cells dropped by Scrublet’s internal filtering get missing scores, not an error."""
    adata = pbmc200()
    adata.obs["batch"] = 100 * ["a"] + 100 * ["b"]
    # make these too sparse to survive `pp.filter_cells(min_genes=3)`, one per batch
    dropped = [5, 105]
    x = cast("CSBase", adata.X).tolil()
    x[dropped, :] = 0
    adata.X = x.tocsr()

    sc.pp.scrublet(adata, use_approx_neighbors=False, batch_key=batch_key)

    assert_array_equal(
        np.flatnonzero(adata.obs["doublet_score"].isna()),
        dropped,
    )
    assert adata.obs["predicted_doublet"].isna().sum() == len(dropped)


def _preprocess_for_scrublet(adata: AnnData) -> AnnData:
    adata_pp = adata.copy()
    sc.pp.filter_genes(adata_pp, min_cells=3)
    sc.pp.filter_cells(adata_pp, min_genes=3)
    adata_pp.layers["raw"] = adata_pp.X.copy()
    sc.pp.normalize_total(adata_pp)
    logged = adata_pp.copy()
    sc.pp.log1p(logged)
    sc.pp.highly_variable_genes(logged)
    return adata_pp[:, logged.var["highly_variable"]].copy()


def _create_sim_from_parents(adata: AnnData, parents: np.ndarray) -> AnnData:
    """Simulate doublets based on the randomly selected parents used previously."""
    n_sim = parents.shape[0]
    entries = sparse.coo_matrix(
        (
            np.ones(2 * n_sim),
            (np.repeat(np.arange(n_sim), 2), parents.flat),
        ),
        (n_sim, adata.n_obs),
    )
    # maintain data type, just like the real scrublet function.
    x = (entries @ adata.layers["raw"]).astype(adata.X.dtype)
    return AnnData(
        x,
        var=pd.DataFrame(index=adata.var_names),
        obs={"total_counts": np.ravel(x.sum(axis=1))},
        obsm={"doublet_parents": parents.copy()},
    )


@pytest.mark.parametrize("rng_arg", ["rng", "random_state"])
def test_scrublet_data(rng_arg: Literal["rng", "random_state"]) -> None:
    """Test that Scrublet processing is arranged correctly.

    Check that simulations run on raw data.
    """
    seed = 1234

    # Run Scrublet and let the main function run simulations
    adata_scrublet_auto_sim = pbmc200()
    sc.pp.scrublet(
        adata_scrublet_auto_sim,
        use_approx_neighbors=False,
        **{rng_arg: seed},
    )

    # Now make our own simulated data so we can check the result from function
    # is the same, and by inference that the processing steps have not been
    # broken

    # Replicate the preprocessing steps used by the main function
    adata_obs = _preprocess_for_scrublet(pbmc200())
    # Simulate doublets using the same parents
    adata_sim = _create_sim_from_parents(
        adata_obs, adata_scrublet_auto_sim.uns["scrublet"]["doublet_parents"]
    )

    # Apply the same post-normalisation the Scrublet function would
    sc.pp.normalize_total(adata_obs, target_sum=1e6)
    sc.pp.normalize_total(adata_sim, target_sum=1e6)

    adata_scrublet_manual_sim = adata_obs.copy()
    sc.pp.scrublet(
        adata_scrublet_manual_sim,
        adata_sim=adata_sim,
        use_approx_neighbors=False,
        **{rng_arg: seed},
    )

    # Require that the doublet scores are the same whether simulation is via
    # the main function or manually provided
    assert_allclose(
        adata_scrublet_manual_sim.obs["doublet_score"],
        adata_scrublet_auto_sim.obs["doublet_score"],
        atol=1e-15,
        rtol=1e-15,
    )


@pytest.fixture(scope="module")
def scrub_small_sess() -> AnnData:
    # Reduce size of input for faster test
    adata = pbmc200()
    sc.pp.filter_genes(adata, min_counts=100)

    sc.pp.scrublet(adata, use_approx_neighbors=False)
    return adata


@pytest.fixture
def scrub_small(scrub_small_sess: AnnData):
    return scrub_small_sess.copy()


test_params = {
    "expected_doublet_rate": 0.1,
    "synthetic_doublet_umi_subsampling": 0.8,
    "knn_dist_metric": "manhattan",
    "normalize_variance": False,
    "log_transform": True,
    "mean_center": False,
    "n_prin_comps": 10,
    "n_neighbors": 2,
    "threshold": 0.1,
}


@pytest.mark.parametrize(("param", "value"), test_params.items())
def test_scrublet_params(scrub_small: AnnData, param: str, value: Any):
    """Test that Scrublet args are passed.

    Check that changes to parameters change scrublet results.
    """
    curr = scrub_small.copy()
    sc.pp.scrublet(adata=curr, use_approx_neighbors=False, **{param: value})
    with pytest.raises(AssertionError):
        assert_equal(scrub_small, curr)


@needs.anndata_acc
def test_scrublet_simulate_doublets():
    """Check that doublet simulation runs and simulates some doublets."""
    adata_obs = pbmc200()
    sc.pp.filter_genes(adata_obs, min_cells=3)
    sc.pp.filter_cells(adata_obs, min_genes=3)
    adata_obs.layers["raw"] = adata_obs.X
    sc.pp.normalize_total(adata_obs)
    # HVG needs log'd data, but the simulation uses the un-log'd counts
    sc.pp.log1p(adata_obs, out="layers.log1p")
    sc.pp.highly_variable_genes(adata_obs, use="layers.log1p")
    adata_obs = adata_obs[:, adata_obs.var["highly_variable"]]

    adata_sim = sc.pp.scrublet_simulate_doublets(
        adata_obs, sim_doublet_ratio=0.02, use="layers.raw"
    )

    assert_array_equal(
        adata_sim.obsm["doublet_parents"],
        np.array([[13, 132], [106, 43], [152, 3], [160, 103]]),
    )


def test_scrublet_dtypes() -> None:
    """Test that Scrublet does not change dtypes of existing data.obs cols."""
    adata = pbmc200()
    adata.obs["batch"] = pd.Categorical(100 * ["a"] + 100 * ["b"])

    sc.pp.scrublet(adata, use_approx_neighbors=False, batch_key="batch")

    assert adata.obs["batch"].dtype == "category"


@needs.anndata_acc
def test_scrublet_out() -> None:
    """`out=True` returns the columns, a mapping names them and the `.uns` entry."""
    from anndata.acc import A

    adata = pbmc200()
    res = sc.pp.scrublet(adata, use_approx_neighbors=False, out=True, rng=0)
    assert "doublet_score" not in adata.obs
    assert "scrublet" not in adata.uns
    out = dict(doublet_score=A.obs["p"], predicted_doublet=A.obs["p_doublet"])
    sc.pp.scrublet(adata, use_approx_neighbors=False, out=out, rng=0)
    np.testing.assert_array_equal(adata.obs["p"], res.doublet_score)
    np.testing.assert_array_equal(adata.obs["p_doublet"], res.predicted_doublet)
    assert "threshold" in adata.uns["p"]

    with pytest.raises(TypeError, match=r"must be a `bool` or a mapping"):
        sc.pp.scrublet(adata, use_approx_neighbors=False, out=A.obs["x"], rng=0)


@needs.anndata_acc
def test_scrublet_batch() -> None:
    """`batch` takes references (or resolvable strings) instead of `batch_key`’s column names."""
    from anndata.acc import A

    adata = pbmc200()
    adata.obs["batch"] = 100 * ["a"] + 100 * ["b"]
    kw = dict(use_approx_neighbors=False, rng=0, out=True)
    expected = sc.pp.scrublet(adata, batch_key="batch", **kw)
    for batch in [A.obs["batch"], "obs.batch"]:
        res = sc.pp.scrublet(adata, batch=batch, **kw)
        np.testing.assert_allclose(res.doublet_score, expected.doublet_score)
