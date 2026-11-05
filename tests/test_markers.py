from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs
from testing.scanpy._pytest.params import ARRAY_TYPES_MEM

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


pytestmark = pytest.mark.filterwarnings(
    r"ignore:The function (rank_genes_groups|rank_genes_groups_df) is deprecated:FutureWarning"
)

COLUMNS = [
    "group",
    "reference",
    "gene",
    "score",
    "log_fc",
    "p_value",
    "adj_p_value",
    "frac_group",
    "frac_reference",
]
GROUPS = ["CD14+ Monocyte", "Dendritic", "CD19+ B"]


@pytest.fixture
def adata() -> AnnData:
    return pbmc68k_reduced().raw.to_adata()


def test_table_layout(adata: AnnData) -> None:
    df = sc.tl.markers.wilcoxon(adata, "bulk_labels")

    assert df.columns.tolist() == COLUMNS
    assert df["group"].cat.categories.tolist() == (
        adata.obs["bulk_labels"].cat.categories.tolist()
    )
    assert (df["reference"] == "rest").all()
    assert len(df) == adata.n_vars * adata.obs["bulk_labels"].nunique()
    for _, group_df in df.groupby("group", observed=True):
        assert group_df["score"].is_monotonic_decreasing


def test_n_top(adata: AnnData) -> None:
    full = sc.tl.markers.ttest(adata, "bulk_labels")
    top = sc.tl.markers.ttest(adata, "bulk_labels", n_top=3)

    expected = full.groupby("group", observed=True).head(3).reset_index(drop=True)
    pd.testing.assert_frame_equal(top, expected)


@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
@pytest.mark.parametrize("reference", ["rest", "Dendritic"])
@pytest.mark.parametrize(
    ("method", "func", "kwargs"),
    [
        pytest.param("wilcoxon", "wilcoxon", {}, id="wilcoxon"),
        pytest.param("wilcoxon", "wilcoxon", {"tie_correct": True}, id="tie_correct"),
        pytest.param("t-test", "ttest", {}, id="t-test"),
        pytest.param(
            "t-test_overestim_var",
            "ttest",
            {"overestimate_var": True},
            id="t-test_overestim_var",
        ),
    ],
)
def test_matches_rank_genes_groups(
    adata: AnnData,
    array_type,
    *,
    reference: str,
    method: str,
    func: Literal["wilcoxon", "ttest"],
    kwargs: dict[str, bool],
) -> None:
    adata.X = array_type(adata.X.toarray())
    new = getattr(sc.tl.markers, func)(
        adata, "bulk_labels", groups=GROUPS, reference=reference, **kwargs
    )

    sc.tl.rank_genes_groups(
        adata,
        "bulk_labels",
        groups=GROUPS,
        reference=reference,
        method=method,
        pts=True,
        mean_in_log_space=True,
        tie_correct=kwargs.get("tie_correct", False),
    )
    old = sc.get.rank_genes_groups_df(adata, None)
    merged = old.merge(new, left_on=["group", "names"], right_on=["group", "gene"])

    assert len(merged) == len(new) == len(old)
    for old_col, new_col in [
        ("scores", "score"),
        ("logfoldchanges", "log_fc"),
        ("pvals", "p_value"),
        ("pvals_adj", "adj_p_value"),
        ("pct_nz_group", "frac_group"),
    ]:
        np.testing.assert_allclose(
            merged[old_col].astype(float), merged[new_col], rtol=1e-5, atol=1e-6
        )


def test_ttest_mean_in_log_space_only_affects_log_fc(adata: AnnData) -> None:
    log_space = sc.tl.markers.ttest(adata, "bulk_labels", mean_in_log_space=True)
    linear = sc.tl.markers.ttest(adata, "bulk_labels", mean_in_log_space=False)

    pd.testing.assert_series_equal(log_space["score"], linear["score"])
    assert not np.allclose(log_space["log_fc"], linear["log_fc"])


def test_logreg_two_groups_have_opposite_scores(adata: AnnData) -> None:
    df = sc.tl.markers.logreg(
        adata, "bulk_labels", groups=["CD14+ Monocyte", "Dendritic"], max_iter=500
    )
    scores = df.pivot_table(index="gene", columns="group", values="score")

    assert "p_value" not in df.columns
    np.testing.assert_allclose(
        scores["CD14+ Monocyte"], -scores["Dendritic"], rtol=0, atol=0
    )


def test_does_not_modify_adata(adata: AnnData) -> None:
    adata.obs["labels"] = adata.obs["bulk_labels"].astype(str)
    adata.X.data[:10] = 0
    nnz = adata.X.nnz

    sc.tl.markers.wilcoxon(adata, "labels")

    assert adata.obs["labels"].dtype != "category"
    assert adata.X.nnz == nnz
    assert "labels_masks" not in adata.uns


def test_mask_and_layer(adata: AnnData) -> None:
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[:20] = True
    adata.layers["shuffled"] = adata.X[:, ::-1]

    masked = sc.tl.markers.wilcoxon(adata, "bulk_labels", mask=mask)
    layered = sc.tl.markers.wilcoxon(adata, "bulk_labels", layer="shuffled")
    plain = sc.tl.markers.wilcoxon(adata, "bulk_labels")

    assert set(masked["gene"]) == set(adata.var_names[mask])
    key = ["group", "gene"]
    layered_scores = layered.set_index(key)["score"]
    plain_scores = plain.set_index(key)["score"].loc[layered_scores.index]
    assert not np.allclose(layered_scores, plain_scores)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param(
            dict(groups=["Dendritic"], reference="Dendritic"),
            r"at least one group other than `reference`",
            id="only_reference",
        ),
        pytest.param(dict(groups=["nope"]), r"not categories", id="unknown_group"),
        pytest.param(dict(reference="nope"), r"needs to be one of", id="bad_ref"),
        pytest.param(dict(groups="Dendritic"), r"sequence of groups", id="str"),
    ],
)
def test_invalid_groups(adata: AnnData, kwargs: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        sc.tl.markers.wilcoxon(adata, "bulk_labels", **kwargs)


@needs.illico
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
@pytest.mark.parametrize("reference", ["rest", "Dendritic"])
def test_illico_backend(adata: AnnData, reference: str) -> None:
    numba_res = sc.tl.markers.wilcoxon(
        adata, "bulk_labels", reference=reference, backend="numba"
    )
    illico_res = sc.tl.markers.wilcoxon(
        adata, "bulk_labels", reference=reference, backend="illico"
    )

    merged = numba_res.merge(illico_res, on=["group", "gene"])
    np.testing.assert_allclose(merged["p_value_x"], merged["p_value_y"], atol=1e-6)


def test_preset_defaults() -> None:
    v1, v2 = sc.Preset.ScanpyV1.markers, sc.Preset.ScanpyV2Preview.markers
    assert (v1.mean_in_log_space, v1.wilcoxon_backend) == (True, "numba")
    assert (v2.mean_in_log_space, v2.wilcoxon_backend) == (False, "illico")
