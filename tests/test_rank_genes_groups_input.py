from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.stats import ttest_ind

import scanpy as sc
from testing.scanpy._pytest.marks import needs


@needs.anndata_acc
@pytest.mark.parametrize(
    ("use", "has_raw", "expected"),
    [
        (None, True, "raw"),
        (None, False, "X"),
        ("X", True, "X"),
        ("layers.alt", True, "alt"),
    ],
)
@pytest.mark.parametrize("mean_in_log_space", [False, True])
def test_rank_genes_groups_input(use, has_raw, expected, mean_in_log_space):
    # Each source has a distinct Welch statistic and fold change for A vs B.
    counts = {
        "X": [1, 2, 3, 8, 10, 12],
        "raw": [8, 10, 12, 1, 2, 3],
        "alt": [1, 5, 15, 2, 3, 4],
    }
    matrices = {k: np.log1p(v)[:, None] for k, v in counts.items()}
    adata = AnnData(matrices["X"])
    adata.obs["condition"] = pd.Categorical(["A"] * 3 + ["B"] * 3)
    if has_raw:
        adata.raw = AnnData(matrices["raw"], obs=adata.obs.copy())
    adata.layers["alt"] = matrices["alt"]
    sc.tl.rank_genes_groups(
        adata,
        "condition",
        groups=["A"],
        reference="B",
        method="t-test",
        use=use,
        mean_in_log_space=mean_in_log_space,
    )
    result = sc.get.rank_genes_groups_df(adata, group="A").iloc[0]
    a, b = np.split(matrices[expected].ravel(), 2)
    test = ttest_ind(a, b, equal_var=False)
    np.testing.assert_allclose(
        [result.scores, result.pvals], [test.statistic, test.pvalue], rtol=1e-6
    )
    if mean_in_log_space:
        a_mean, b_mean = np.expm1([a.mean(), b.mean()])
    else:
        a_mean, b_mean = np.expm1(a).mean(), np.expm1(b).mean()
    np.testing.assert_allclose(
        result.logfoldchanges, np.log2((a_mean + 1e-9) / (b_mean + 1e-9)), rtol=1e-6
    )
