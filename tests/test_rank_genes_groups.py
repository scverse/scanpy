from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict, cast

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.random import binomial, negative_binomial, seed
from scipy.stats import mannwhitneyu

import scanpy as sc
from scanpy._compat import CSBase
from scanpy._utils import select_groups
from scanpy.get import rank_genes_groups_df
from scanpy.tools import rank_genes_groups
from scanpy.tools._rank_genes_groups import _RankGenes
from testing.scanpy._helpers import random_mask
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.params import ARRAY_TYPES, ARRAY_TYPES_MEM

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    from numpy.lib.npyio import NpzFile
    from numpy.typing import NDArray

HERE = Path(__file__).parent
DATA_PATH = HERE / "_data"


# We test results for a simple generic example
# Tests are conducted for sparse and non-sparse AnnData objects.
# Due to minor changes in multiplication implementation for sparse and non-sparse objects,
# results differ (very) slightly


@pytest.mark.parametrize("array_type", ARRAY_TYPES)
def get_example_data(array_type: Callable[[np.ndarray], Any]) -> AnnData:
    # create test object
    adata = AnnData(
        np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20)))
    )
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(
        binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5))
    )

    adata.X = array_type(adata.X)

    # Create cluster according to groups
    adata.obs["true_groups"] = pd.Categorical(
        np.concatenate((np.zeros((10,), dtype=int), np.ones((90,), dtype=int)))
    )

    return adata


class Expected(TypedDict):
    names: NDArray[np.str_]
    scores: NDArray[np.floating]


def get_true_scores(method: Literal["t-test", "wilcoxon"]) -> Expected:
    path = DATA_PATH / f"objs-{method}.npz"
    with (
        path.open("rb") as f,
        cast("NpzFile", np.load(f, allow_pickle=False)) as z,
    ):
        expected = dict(z)
    return Expected(names=expected["names"].astype("T"), scores=expected["scores"])


# TODO: Make dask compatible
@pytest.mark.parametrize("method", ["t-test", "wilcoxon"])
@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_results(
    subtests: pytest.Subtests, array_type, method: Literal["t-test", "wilcoxon"]
) -> None:
    seed(1234)
    adata = get_example_data(array_type)
    assert adata.raw is None  # Assumption for later checks
    expected = get_true_scores(method)
    # no clue why we did this: https://github.com/scverse/scanpy/commit/7f10fa3138374bbc664776c6aae1c0e05cf2c5cf
    n = 7 if method == "wilcoxon" else None

    rank_genes_groups(adata, "true_groups", n_genes=20, method=method)
    results = adata.uns["rank_genes_groups"]

    for g in range(expected["names"].shape[0]):
        with subtests.test(group=g):
            assert np.allclose(expected["scores"][g, :n], results["scores"][str(g)][:n])
            assert np.array_equal(
                expected["names"][g, :n], results["names"][str(g)][:n]
            )
    assert results["params"]["use_raw"] is False


@pytest.mark.parametrize("method", ["t-test", "wilcoxon"])
@pytest.mark.parametrize("array_type", ARRAY_TYPES_MEM)
def test_results_layers(
    subtests: pytest.Subtests, array_type, method: Literal["t-test", "wilcoxon"]
) -> None:
    seed(1234)
    adata = get_example_data(array_type)
    adata.layers["to_test"] = adata.X.copy()
    x = adata.X.tolil() if isinstance(adata.X, CSBase) else adata.X
    mask = np.random.randint(0, 2, adata.shape, dtype=bool)
    x[mask] = 0
    adata.X = array_type(x)
    scores = get_true_scores(method)["scores"]

    with subtests.test("layer"):
        rank_genes_groups(
            adata,
            "true_groups",
            method=method,
            layer="to_test",
            use_raw=None if method == "wilcoxon" else False,
            n_genes=20,
        )
        assert adata.uns["rank_genes_groups"]["params"]["use_raw"] is False
        for g in range(scores.shape[0]):
            np.testing.assert_allclose(
                scores[g, :7],
                adata.uns["rank_genes_groups"]["scores"][str(g)][:7],
                rtol=1e-5,  # default of np.allclose
            )

    with subtests.test("X"):
        rank_genes_groups(adata, "true_groups", method=method, n_genes=20)
        for g in range(scores.shape[0]):
            assert not np.allclose(
                scores[g, :7], adata.uns["rank_genes_groups"]["scores"][str(g)][:7]
            )


def test_rank_genes_groups_use_raw():
    # https://github.com/scverse/scanpy/issues/1929
    pbmc = pbmc68k_reduced()
    assert pbmc.raw is not None

    sc.tl.rank_genes_groups(pbmc, groupby="bulk_labels", use_raw=True)

    pbmc = pbmc68k_reduced()
    del pbmc.raw
    assert pbmc.raw is None

    with pytest.raises(
        ValueError, match=r"Received `use_raw=True`, but `adata\.raw` is empty"
    ):
        sc.tl.rank_genes_groups(pbmc, groupby="bulk_labels", use_raw=True)


def test_singlets():
    pbmc = pbmc68k_reduced()
    pbmc.obs["louvain"] = pbmc.obs["louvain"].cat.add_categories(["11"])
    pbmc.obs[0, "louvain"] = "11"

    with pytest.raises(ValueError, match=rf"Could not calculate statistics.*{'11'}"):
        rank_genes_groups(pbmc, groupby="louvain")


def test_emptycat():
    pbmc = pbmc68k_reduced()
    pbmc.obs["louvain"] = pbmc.obs["louvain"].cat.add_categories(["11"])

    with pytest.raises(ValueError, match=rf"Could not calculate statistics.*{'11'}"):
        rank_genes_groups(pbmc, groupby="louvain")


def test_log1p_save_restore(tmp_path):
    """Tests the sequence log1p→save→load→rank_genes_groups."""
    from anndata import read_h5ad

    pbmc = pbmc68k_reduced()
    pbmc.X = pbmc.raw.X
    sc.pp.log1p(pbmc)

    path = tmp_path / "test.h5ad"
    pbmc.write(path)

    pbmc = read_h5ad(path)

    sc.tl.rank_genes_groups(pbmc, groupby="bulk_labels", use_raw=True)


def test_wilcoxon_symmetry():
    pbmc = pbmc68k_reduced()

    rank_genes_groups(
        pbmc,
        groupby="bulk_labels",
        groups=["CD14+ Monocyte", "Dendritic"],
        reference="Dendritic",
        method="wilcoxon",
        rankby_abs=True,
    )
    assert pbmc.uns["rank_genes_groups"]["params"]["use_raw"] is True

    stats_mono = (
        rank_genes_groups_df(pbmc, group="CD14+ Monocyte")
        .drop(columns="names")
        .to_numpy()
    )

    rank_genes_groups(
        pbmc,
        groupby="bulk_labels",
        groups=["CD14+ Monocyte", "Dendritic"],
        reference="CD14+ Monocyte",
        method="wilcoxon",
        rankby_abs=True,
    )

    stats_dend = (
        rank_genes_groups_df(pbmc, group="Dendritic").drop(columns="names").to_numpy()
    )

    assert np.allclose(np.abs(stats_mono), np.abs(stats_dend))


@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
@pytest.mark.parametrize("reference", [True, False], ids=["ref", "rest"])
def test_wilcoxon_tie_correction(*, reference: bool) -> None:
    pbmc = pbmc68k_reduced()

    groups = ["CD14+ Monocyte", "Dendritic"]
    groupby = "bulk_labels"

    _, groups_masks = select_groups(pbmc, groups, groupby)

    if reference:
        ref = groups[1]
        mask_rest = groups_masks[1]
    else:
        ref = "rest"
        mask_rest = ~groups_masks[0]
        groups = groups[:1]

    assert isinstance(pbmc.raw.X, CSBase)
    x = pbmc.raw.X[groups_masks[0]].toarray()
    y = pbmc.raw.X[mask_rest].toarray()

    pvals = mannwhitneyu(x, y, use_continuity=False, alternative="two-sided").pvalue
    pvals[np.isnan(pvals)] = 1.0

    test_obj = _RankGenes(pbmc, groups, groupby, reference=ref)
    test_obj.compute_statistics("wilcoxon", tie_correct=True)

    np.testing.assert_allclose(test_obj.stats[groups[0]]["pvals"], pvals, atol=1e-5)


def test_rank_gene_groups_violin_gene_symbols():
    adata = sc.read_h5ad(DATA_PATH / "t-cells.h5ad")

    t_celltypes = ["Naive CD4+ T cells", "Naive CD8+ T cells"]
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        groups=t_celltypes,
        reference="rest",
        method="wilcoxon",
        use_raw=False,
    )

    # We get the KeyError here
    sc.pl.rank_genes_groups_violin(
        adata, groups=t_celltypes, n_genes=15, strip=False, gene_symbols="gene_symbols"
    )


def test_wilcoxon_huge_data(monkeypatch):
    max_size = 300
    adata = pbmc68k_reduced()
    monkeypatch.setattr(sc.tl._rank_genes_groups, "_CONST_MAX_SIZE", max_size)
    rank_genes_groups(adata, groupby="bulk_labels", method="wilcoxon")


@pytest.mark.parametrize(
    ("n_genes_add", "n_genes_out_add"),
    [pytest.param(0, 0, id="equal"), pytest.param(2, 1, id="more")],
)
def test_mask_n_genes(n_genes_add, n_genes_out_add):
    """Check if no. genes in output is correct.

    1. =n_genes when n_genes<sum(mask)
    2. =sum(mask) when n_genes>sum(mask)
    """
    pbmc = pbmc68k_reduced()
    mask_var = np.zeros(pbmc.shape[1]).astype(bool)
    mask_var[:6].fill(True)  # noqa: FBT003
    no_genes = sum(mask_var) - 1

    rank_genes_groups(
        pbmc,
        mask_var=mask_var,
        groupby="bulk_labels",
        groups=["CD14+ Monocyte", "Dendritic"],
        reference="CD14+ Monocyte",
        n_genes=no_genes + n_genes_add,
        method="wilcoxon",
    )

    assert len(pbmc.uns["rank_genes_groups"]["scores"]) == no_genes + n_genes_out_add


def test_mask_not_equal():
    """Check that mask is applied successfully to data set where test statistics are already available (test stats overwritten)."""
    pbmc = pbmc68k_reduced()
    mask_var = random_mask(pbmc.shape[1])
    n_genes = sum(mask_var)

    run = partial(
        rank_genes_groups,
        pbmc,
        groupby="bulk_labels",
        groups=["CD14+ Monocyte", "Dendritic"],
        reference="CD14+ Monocyte",
        method="wilcoxon",
    )

    run(n_genes=n_genes)
    no_mask = pbmc.uns["rank_genes_groups"]["names"]

    run(mask_var=mask_var)
    with_mask = pbmc.uns["rank_genes_groups"]["names"]

    assert not np.array_equal(no_mask, with_mask)
