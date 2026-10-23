"""Tests for the `use`/`out` parameters replacing `layer`/`obsm`/`use_rep`/`inplace`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from anndata.tests.helpers import asarray
from scipy import sparse

import scanpy as sc
from scanpy.neighbors import NeighborsOut
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES_MEM

if TYPE_CHECKING:
    from collections.abc import Callable

pytestmark = needs.anndata_acc


@pytest.fixture
def adata() -> AnnData:
    x = np.random.default_rng(0).poisson(3, (20, 5)).astype("float32")
    adata = AnnData(x.copy())
    adata.layers["counts"] = x.copy()
    adata.obs["batch"] = np.tile([0.0, 1.0], 10)
    return adata


def _orig(adata: AnnData) -> np.ndarray:
    return asarray(adata.layers["counts"])


# `use`/`out` on the same-shape transforms.
# Each entry: id -> (call, name of the `.X`-derived result check)
TRANSFORMS: dict[str, Callable[..., object]] = {
    "scale": lambda adata, **kw: sc.pp.scale(adata, zero_center=True, **kw),
    "log1p": sc.pp.log1p,
    "normalize_total": lambda adata, **kw: sc.pp.normalize_total(
        adata, target_sum=1, **kw
    ),
    "regress_out": lambda adata, **kw: sc.pp.regress_out(adata, "batch", **kw),
}
# transforms with several outputs, whose `out` is a mapping (name of the matrix output)
MULTI = {"scale": "scaled", "normalize_total": "normalized"}


def _out(name: str, acc: object) -> object:
    return {MULTI[name]: acc} if name in MULTI else acc


@pytest.mark.filterwarnings("ignore:zero-centering a sparse:UserWarning")
@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_writes_elsewhere(
    adata: AnnData, name: str, array_type: Callable[[np.ndarray], object]
) -> None:
    """`out` writes to its own destination and leaves `use`’s source alone."""
    from anndata.acc import A

    adata.layers["counts"] = array_type(adata.layers["counts"])
    before = _orig(adata).copy()
    TRANSFORMS[name](adata, use=A.layers["counts"], out=_out(name, A.layers["result"]))
    assert "result" in adata.layers
    np.testing.assert_array_equal(_orig(adata), before)
    assert not np.allclose(asarray(adata.layers["result"]), before)


@pytest.mark.filterwarnings("ignore:zero-centering a sparse:UserWarning")
@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_none_returns(
    adata: AnnData, name: str, array_type: Callable[[np.ndarray], object]
) -> None:
    """`out=True` returns the result without touching `adata`."""
    from anndata.acc import A

    adata.layers["counts"] = array_type(adata.layers["counts"])
    before = _orig(adata).copy()
    res = TRANSFORMS[name](adata, use=A.layers["counts"], out=True)
    if name in MULTI:
        res = getattr(res, MULTI[name])
    assert res is not None
    assert np.shape(res) == adata.shape
    np.testing.assert_array_equal(_orig(adata), before)
    assert "result" not in adata.layers


@pytest.mark.parametrize("name", TRANSFORMS)
def test_out_default_writes_back(adata: AnnData, name: str) -> None:
    """Without `out`, the result goes back where `use` read it from."""
    from anndata.acc import A

    before = _orig(adata).copy()
    x_before = np.asarray(adata.X).copy()
    assert TRANSFORMS[name](adata, use=A.layers["counts"]) is None
    assert not np.allclose(_orig(adata), before)
    np.testing.assert_array_equal(np.asarray(adata.X), x_before)


@pytest.mark.parametrize("name", TRANSFORMS)
def test_use_str_is_resolved(adata: AnnData, name: str) -> None:
    """A `use` string is always an `anndata.acc` spec, regardless of preset."""
    from anndata.acc import A

    expected = adata.copy()
    TRANSFORMS[name](
        expected, use=A.layers["counts"], out=_out(name, A.layers["result"])
    )
    TRANSFORMS[name](adata, use="layers.counts", out=_out(name, "layers.result"))
    np.testing.assert_allclose(adata.layers["result"], expected.layers["result"])


def test_use_bare_str_rejected(adata: AnnData) -> None:
    """A bare (non-period-separated) `use` string is not a legacy key."""
    with pytest.raises(ValueError, match=r"Cannot parse accessor"):
        sc.pp.log1p(adata, use="counts")


@pytest.mark.parametrize(
    ("call", "legacy"),
    [
        pytest.param(sc.pp.log1p, "layer", id="log1p"),
        pytest.param(
            lambda a, **kw: sc.pp.scale(a, zero_center=True, **kw), "layer", id="scale"
        ),
        pytest.param(lambda a, **kw: sc.pp.pca(a, n_comps=2, **kw), "layer", id="pca"),
    ],
)
def test_legacy_layer_deprecated(
    adata: AnnData, call: Callable[..., object], legacy: str
) -> None:
    """The legacy `layer` argument still works, but warns."""
    with pytest.warns(FutureWarning, match=rf"argument {legacy} is deprecated"):
        call(adata, **{legacy: "counts"})


@pytest.mark.parametrize(
    "call",
    [
        pytest.param(
            lambda a: sc.pp.calculate_qc_metrics(a, percent_top=None, use_raw=False),
            id="qc",
        ),
        pytest.param(
            lambda a: sc.tl.score_genes(a, list(a.var_names[:2]), use_raw=False),
            id="score_genes",
        ),
        pytest.param(
            lambda a: sc.tl.rank_genes_groups(a, "group", use_raw=False),
            id="rank_genes_groups",
        ),
        pytest.param(lambda a: sc.metrics.gearys_c(a, use_raw=False), id="gearys_c"),
    ],
)
def test_use_raw_deprecated(adata: AnnData, call: Callable[[AnnData], object]) -> None:
    """`.raw` is going away, so `use_raw` warns wherever it's passed."""
    adata.obs["group"] = np.tile(["a", "b"], 10)
    sc.pp.neighbors(adata, n_neighbors=3)
    with pytest.warns(FutureWarning, match=r"argument use_raw is deprecated"):
        call(adata)


@pytest.mark.parametrize(
    "call",
    [
        pytest.param(sc.pp.log1p, id="log1p"),
        pytest.param(
            lambda d, **kw: sc.pp.scale(d, zero_center=True, **kw), id="scale"
        ),
    ],
)
@pytest.mark.parametrize(
    "data",
    [
        pytest.param(np.ones((4, 3), dtype="float32"), id="ndarray"),
        pytest.param(sparse.csr_matrix(np.ones((4, 3), dtype="float32")), id="csr"),  # noqa: TID251
        pytest.param([[1.0, 2.0], [3.0, 4.0]], id="unregistered"),
    ],
)
@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param(dict(use="layers.x"), r"Arguments \['use'\]", id="use"),
        pytest.param(dict(out="layers.x"), r"Argument `out`", id="out"),
    ],
)
@pytest.mark.filterwarnings("ignore:.*is deprecated:FutureWarning")
def test_anndata_only_args_rejected(
    call: Callable[..., object], data: object, kwargs: dict[str, str], match: str
) -> None:
    """`use`/`out` need an `AnnData`, and every input kind says so the same way.

    The check runs before `singledispatch` routes, so registered types reach it too.
    """
    with pytest.raises(TypeError, match=match):
        call(data, **kwargs)


def test_use_and_layer_conflict(adata: AnnData) -> None:
    from anndata.acc import A

    with (
        pytest.warns(FutureWarning, match=r"argument layer is deprecated"),
        pytest.raises(TypeError, match=r"`use` cannot be combined with"),
    ):
        sc.pp.log1p(adata, use=A.layers["counts"], layer="counts")


def test_score_genes_out(adata: AnnData) -> None:
    """`score_genes`’ `out` references an `.obs` column; `score_name` is its deprecated version."""
    from anndata.acc import A

    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    genes = list(adata.var_names[:2])

    sc.tl.score_genes(adata, genes, out=A.obs["mine"])
    assert "mine" in adata.obs

    scores = sc.tl.score_genes(adata, genes, out=True)
    assert scores.shape == (adata.n_obs,)
    assert "score" not in adata.obs

    with pytest.warns(FutureWarning, match=r"argument score_name is deprecated"):
        sc.tl.score_genes(adata, genes, score_name="legacy")
    assert "legacy" in adata.obs

    # `score_name` names a column, `out` takes an accessor: not interchangeable
    with pytest.raises(ValueError, match=r"Cannot parse accessor"):
        sc.tl.score_genes(adata, genes, out="legacy")


def test_score_genes_out_and_score_name(adata: AnnData) -> None:
    from anndata.acc import A

    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    with (
        pytest.warns(FutureWarning, match=r"argument score_name is deprecated"),
        pytest.raises(TypeError, match=r"Pass either `out` or `score_name`"),
    ):
        sc.tl.score_genes(
            adata, list(adata.var_names[:2]), out=A.obs["a"], score_name="b"
        )


@pytest.mark.filterwarnings("ignore:.*argument key_added is deprecated:FutureWarning")
@pytest.mark.parametrize("kw", [{}, dict(key_added="k")], ids=["plain", "key_added"])
def test_normalize_total_inplace_deprecated(
    adata: AnnData, kw: dict[str, object]
) -> None:
    """`inplace=False` is the deprecated way to say `out=True`, which ignores `key_added`."""
    with pytest.warns(FutureWarning, match=r"argument inplace is deprecated"):
        res = sc.pp.normalize_total(adata, target_sum=1, inplace=False, **kw)
    assert set(res) == {"X", "norm_factor"}
    assert "k" not in adata.obs
    np.testing.assert_allclose(np.asarray(res["X"]).sum(axis=1), 1, rtol=1e-6)


def test_pca_use(adata: AnnData) -> None:
    """`pca` takes `use`, and `sc.Prefix` names all of obsm, varm and uns."""
    from anndata.acc import A

    sc.pp.pca(adata, n_comps=2, use=A.layers["counts"], out=sc.Prefix("mine"))
    assert "mine_pca" in adata.obsm
    assert "mine_pca" in adata.varm
    assert "mine_pca" in adata.uns


@pytest.fixture
def adata_neighbors(adata: AnnData) -> AnnData:
    sc.pp.neighbors(adata, n_neighbors=3)
    return adata


def test_out_single_acc(adata_neighbors: AnnData) -> None:
    """An accessor `out` names the output, and `.uns` follows it."""
    from anndata.acc import A

    sc.tl.umap(adata_neighbors, out=A.obsm["mine"])
    assert "mine" in adata_neighbors.obsm
    assert "params" in adata_neighbors.uns["mine"]
    assert set(adata_neighbors.obsm) == {"mine"}  # nothing at the default key


def test_out_single_prefix(adata_neighbors: AnnData) -> None:
    sc.tl.umap(adata_neighbors, out=sc.Prefix("mine"))
    assert "mine_umap" in adata_neighbors.obsm
    assert "mine_umap" in adata_neighbors.uns


def test_out_single_true(adata_neighbors: AnnData) -> None:
    """`out=True` returns the output and writes nothing, not even `.uns`."""
    before = set(adata_neighbors.uns)
    emb = sc.tl.umap(adata_neighbors, out=True)
    assert emb.shape == (adata_neighbors.n_obs, 2)
    assert not adata_neighbors.obsm
    assert set(adata_neighbors.uns) == before


def test_out_single_wrong_kind(adata_neighbors: AnnData) -> None:
    from anndata.acc import A

    with pytest.raises(TypeError, match=r"`out` must point to `\.obsm`"):
        sc.tl.umap(adata_neighbors, out=A.obs["mine"])


def test_out_multi_mapping(adata: AnnData) -> None:
    """A mapping `out` only writes what it names, and `.uns` follows the primary output."""
    from anndata.acc import A

    sc.pp.neighbors(adata, n_neighbors=3, out=dict(connectivities=A.obsp["c"]))
    assert set(adata.obsp) == {"c"}
    assert adata.uns["c"]["connectivities_key"] == "c"
    assert "distances_key" not in adata.uns["c"]
    sc.tl.umap(adata, use=dict(connectivities=A.obsp["c"]))  # readers find it


def test_out_multi_prefix(adata: AnnData) -> None:
    sc.pp.neighbors(adata, n_neighbors=3, out=sc.Prefix("p"))
    assert set(adata.obsp) == {"p_connectivities", "p_distances"}
    assert "p" not in adata.uns
    assert adata.uns["p_neighbors"]["connectivities_key"] == "p_connectivities"
    assert adata.uns["p_neighbors"]["distances_key"] == "p_distances"


@pytest.mark.parametrize(
    ("call", "expected"),
    [
        pytest.param(sc.pp.log1p, dict(layers={"p_log1p"}), id="log1p"),
        pytest.param(
            lambda a, **kw: sc.pp.scale(a, zero_center=True, **kw),
            dict(layers={"p_scaled"}, var={"p_mean", "p_std"}),
            id="scale",
        ),
        pytest.param(
            lambda a, **kw: sc.pp.normalize_total(a, target_sum=1, **kw),
            dict(layers={"p_normalized"}, obs={"p_norm_factor"}),
            id="normalize_total",
        ),
        pytest.param(
            lambda a, **kw: sc.tl.score_genes(a, list(a.var_names[:2]), **kw),
            dict(obs={"p_score"}),
            id="score_genes",
        ),
    ],
)
@pytest.mark.parametrize("use", ["X", "obsm"])
def test_out_transform_prefix(
    adata: AnnData,
    call: Callable[..., object],
    expected: dict[str, set[str]],
    use: str,
) -> None:
    """`sc.Prefix('p')` names each output `p_{default name}`, matrices follow `use`’s kind."""
    from anndata.acc import A

    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    kw = {}
    if use == "obsm":
        if "layers" not in expected:
            pytest.skip("no matrix output")
        adata.obsm["m"] = np.asarray(adata.X).copy()
        kw["use"] = A.obsm["m"]
        expected = dict(expected, obsm=expected["layers"] | {"m"}, layers=set())
    x_before = np.asarray(adata.X).copy()
    obs_before, var_before = set(adata.obs), set(adata.var)
    call(adata, **kw, out=sc.Prefix("p"))
    assert set(adata.layers) - {None, "counts"} == expected.get("layers", set())
    assert set(adata.obsm) == expected.get("obsm", set())
    assert set(adata.obs) - obs_before == expected.get("obs", set())
    assert set(adata.var) - var_before == expected.get("var", set())
    np.testing.assert_array_equal(adata.X, x_before)


def test_out_multi_true(adata: AnnData) -> None:
    res = sc.pp.neighbors(adata, n_neighbors=3, out=True)
    assert isinstance(res, NeighborsOut)
    assert res.distances.shape == res.connectivities.shape == (adata.n_obs,) * 2
    assert not adata.obsp
    assert "neighbors" not in adata.uns


@pytest.mark.parametrize(
    ("out", "match"),
    [
        pytest.param(
            "obsp.c", r"`out` must be a `bool`, `sc.Prefix`, or a mapping", id="bare"
        ),
        pytest.param(
            dict(nope="obsp.c"), r"Unknown output\(s\) \['nope'\]", id="unknown"
        ),
        pytest.param(
            dict(distances="obsm.c"),
            r"Output 'distances' must be written to one of \['obsp'\]",
            id="kind",
        ),
        pytest.param(
            dict(distances="obsp.d"),
            r"`out` must include 'connectivities'",
            id="no-primary",
        ),
    ],
)
def test_out_multi_invalid(adata: AnnData, out: object, match: str) -> None:
    with pytest.raises(TypeError, match=match):
        sc.pp.neighbors(adata, n_neighbors=3, out=out)


@pytest.mark.parametrize(
    ("call", "equiv"),
    [
        pytest.param(sc.tl.umap, lambda a: a.obsm["k"], id="umap"),
        pytest.param(
            lambda a, **kw: sc.pp.neighbors(a, n_neighbors=3, **kw),
            lambda a: sc.Prefix("k"),
            id="neighbors",
        ),
        pytest.param(
            lambda a, **kw: sc.pp.pca(a, n_comps=2, **kw),
            lambda a: dict(embedding=a.obsm["k"], loadings=a.varm["k"]),
            id="pca",
        ),
        pytest.param(
            sc.pp.normalize_total,
            lambda a: dict(normalized=a.X, norm_factor=a.obs["k"]),
            id="normalize_total",
        ),
    ],
)
def test_key_added_deprecated(
    request: pytest.FixtureRequest,
    adata_neighbors: AnnData,
    call: Callable[..., object],
    equiv: Callable[[object], object],
) -> None:
    """`key_added='k'` still works like its documented `out` equivalent, but warns, and conflicts with `out`."""
    from anndata.acc import A

    expected = adata_neighbors.copy()
    call(expected, out=equiv(A))
    if request.node.callspec.id == "neighbors":
        # the only difference: `key_added` names the `.uns` entry `k`, not `k_neighbors`
        expected.uns["k"] = expected.uns.pop("k_neighbors")
    with pytest.warns(FutureWarning, match=r"argument key_added is deprecated"):
        call(adata_neighbors, key_added="k")
    assert adata_neighbors.uns.keys() == expected.uns.keys()
    assert adata_neighbors.obsm.keys() == expected.obsm.keys()
    assert adata_neighbors.obsp.keys() == expected.obsp.keys()
    assert adata_neighbors.varm.keys() == expected.varm.keys()
    assert adata_neighbors.obs.keys().tolist() == expected.obs.keys().tolist()

    with (
        pytest.warns(FutureWarning, match=r"argument key_added is deprecated"),
        pytest.raises(TypeError, match=r"Pass either `out` or `key_added`"),
    ):
        call(adata_neighbors, key_added="k", out=equiv(A))


def test_pca_out_true(adata: AnnData) -> None:
    """`out=True` returns embedding and (mask-padded) loadings, writing nothing."""
    adata.var["highly_variable"] = [True, True, True, False, True]
    res = sc.pp.pca(adata, n_comps=2, out=True)
    assert isinstance(res, sc.pp.PcaOut)
    assert res.embedding.shape == (adata.n_obs, 2)
    assert res.loadings.shape == (adata.n_vars, 2)
    np.testing.assert_array_equal(res.loadings[3], 0)
    assert not adata.obsm
    assert not adata.varm
    assert not adata.uns


def test_pca_out_mapping(adata: AnnData) -> None:
    """A mapping writes only what it names, `.uns` follows the embedding."""
    from anndata.acc import A

    sc.pp.pca(adata, n_comps=2, out=dict(embedding=A.obsm["e"]))
    assert set(adata.obsm) == {"e"}
    assert not adata.varm
    assert "variance" in adata.uns["e"]


def test_pca_copy_out_true(adata: AnnData) -> None:
    with (
        pytest.warns(FutureWarning, match=r"argument copy is deprecated"),
        pytest.raises(ValueError, match=r"`copy=True` cannot be used with `out=True`"),
    ):
        sc.pp.pca(adata, n_comps=2, out=True, copy=True)


def test_normalize_total_out(adata: AnnData) -> None:
    from anndata.acc import A

    res = sc.pp.normalize_total(adata, target_sum=1, out=True)
    assert isinstance(res, sc.pp.NormalizeTotalOut)
    assert res.norm_factor.shape == (adata.n_obs,)
    np.testing.assert_array_equal(adata.X, _orig(adata))

    sc.pp.normalize_total(adata, target_sum=1, out=dict(norm_factor=A.obs["nf"]))
    np.testing.assert_allclose(adata.obs["nf"], res.norm_factor)
    np.testing.assert_array_equal(adata.X, _orig(adata))  # `X` not named

    sc.pp.normalize_total(adata, target_sum=1, out=sc.Prefix("p"))
    np.testing.assert_allclose(adata.obs["p_norm_factor"], res.norm_factor)
    np.testing.assert_allclose(adata.layers["p_normalized"], res.normalized)
    np.testing.assert_array_equal(adata.X, _orig(adata))


def test_normalize_total_bare_acc_rejected(adata: AnnData) -> None:
    from anndata.acc import A

    with pytest.raises(
        TypeError, match=r"`out` must be a `bool`, `sc.Prefix`, or a mapping"
    ):
        sc.pp.normalize_total(adata, out=A.layers["x"])


def test_scale_out_true_writes_nothing(adata: AnnData) -> None:
    var_before = adata.var.copy()
    res = sc.pp.scale(adata, zero_center=True, out=True)
    assert isinstance(res, sc.pp.ScaleOut)
    assert res.scaled.shape == adata.shape
    assert res.mean.shape == res.std.shape == (adata.n_vars,)
    np.testing.assert_array_equal(adata.X, _orig(adata))
    assert adata.var.equals(var_before)


def test_scale_out_mapping(adata: AnnData) -> None:
    from anndata.acc import A

    sc.pp.scale(adata, zero_center=True, out=dict(mean=A.var["m"]))
    assert list(adata.var.columns) == ["m"]
    np.testing.assert_array_equal(adata.X, _orig(adata))


def test_pearson_residuals_inplace_deprecated(adata: AnnData) -> None:
    """`inplace=False` still returns a dict, `out=True` the matrix."""
    res = sc.experimental.pp.normalize_pearson_residuals(adata, out=True)
    with pytest.warns(FutureWarning, match=r"argument inplace is deprecated"):
        legacy = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    assert {"X", "theta", "clip", "computed_on"} <= set(legacy)
    np.testing.assert_array_equal(legacy["X"], res)
    assert not adata.uns


def test_pearson_residuals_pca_kwargs_out(adata: AnnData) -> None:
    """`kwargs_pca` takes `out`, and the deprecated `key_added` still works."""
    from anndata.acc import A

    sc.experimental.pp.normalize_pearson_residuals_pca(
        adata, n_comps=2, kwargs_pca=dict(out=dict(embedding=A.obsm["e"]))
    )
    assert set(adata.obsm) == {"e"}
    assert not adata.varm
    assert "e" in adata.uns
    with pytest.warns(FutureWarning, match=r"argument key_added is deprecated"):
        sc.experimental.pp.normalize_pearson_residuals_pca(
            adata, n_comps=2, kwargs_pca=dict(key_added="k")
        )
    assert adata.varm["k"].shape == (adata.n_vars, 2)
    assert "k" in adata.uns


# Every function writing an accompanying `.uns` entry next to its primary output,
# as (call, `out` naming the primary output `k`, `.uns` key that should follow it,
# `.uns` key for `out=sc.Prefix('p')`).
# When adding an explicit `.uns` slot to `out`, extend these: they list all such functions.
UNS_FOLLOWS_OUT = [
    pytest.param(sc.tl.umap, "obsm.k", "k", "p_umap", id="umap"),
    pytest.param(
        lambda a, **kw: sc.tl.tsne(a, perplexity=5, **kw),
        "obsm.k",
        "k",
        "p_tsne",
        id="tsne",
    ),
    pytest.param(
        lambda a, **kw: sc.tl.diffmap(a, n_comps=3, **kw),
        "obsm.k",
        "k",
        "p_diffmap",
        id="diffmap",
    ),
    pytest.param(
        lambda a, **kw: sc.tl.draw_graph(a, layout="fr", **kw),
        "obsm.k",
        "k",
        "p_graph_fr",
        id="draw_graph",
        marks=needs.igraph,
    ),
    pytest.param(
        lambda a, **kw: sc.tl.leiden(a, flavor="igraph", **kw),
        "obs.k",
        "k",
        "p_leiden",
        id="leiden",
        marks=needs.igraph,
    ),
    pytest.param(
        lambda a, **kw: sc.tl.embedding_density(a, basis="umap", **kw),
        "obs.k",
        "k_params",
        "p_umap_density_params",
        id="embedding_density",
    ),
    pytest.param(
        sc.experimental.pp.normalize_pearson_residuals,
        "layers.k",
        "k",
        "p_pearson_residuals",
        id="normalize_pearson_residuals",
    ),
    pytest.param(
        lambda a, **kw: sc.pp.pca(a, n_comps=2, **kw),
        dict(embedding="obsm.k"),
        "k",
        "p_pca",
        id="pca",
    ),
    pytest.param(
        lambda a, **kw: sc.pp.neighbors(a, n_neighbors=3, **kw),
        dict(connectivities="obsp.k"),
        "k",
        "p_neighbors",
        id="neighbors",
    ),
]


@pytest.mark.parametrize(("call", "out", "uns_key", "uns_key_prefix"), UNS_FOLLOWS_OUT)
def test_uns_follows_out(
    adata_neighbors: AnnData,
    call: Callable[..., object],
    out: str | dict[str, str],
    uns_key: str,
    uns_key_prefix: str,
) -> None:
    """The accompanying `.uns` entry is named after the key `out` gives the primary output."""
    sc.tl.umap(adata_neighbors)  # for `embedding_density`
    adata = adata_neighbors.copy()
    before = set(adata.uns)
    call(adata, out=out)
    assert set(adata.uns) - before == {uns_key}

    adata = adata_neighbors.copy()
    before = set(adata.uns)
    call(adata, out=sc.Prefix("p"))
    assert {k for k in set(adata.uns) - before if k.startswith("p_")} == {
        uns_key_prefix
    }


@pytest.mark.parametrize(
    ("kwargs", "marks_x"),
    [
        pytest.param({}, True, id="default"),
        pytest.param(dict(out="X"), True, id="out_X"),
        pytest.param(dict(use="layers.counts"), False, id="layer_in_place"),
        pytest.param(dict(out="layers.log"), False, id="out_layer"),
        pytest.param(dict(out=sc.Prefix("p")), False, id="prefix"),
        pytest.param(dict(use="layers.counts", out="X"), True, id="layer_to_X"),
    ],
)
def test_log1p_marks_only_x(
    adata: AnnData, kwargs: dict[str, object], *, marks_x: bool
) -> None:
    """`.uns['log1p']` says `.X` is log-transformed, so it’s only set when writing `.X`."""
    sc.pp.log1p(adata, base=2, **kwargs)
    assert ("log1p" in adata.uns) is marks_x
