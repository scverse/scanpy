from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scanpy as sc
from testing.scanpy._pytest.marks import needs

pytestmark = [
    needs.pydeseq2,
    pytest.mark.filterwarnings(
        "ignore:The dispersion trend curve fitting did not converge:UserWarning"
    ),
]

N_GENES = 300
DE_GENES = [str(i) for i in range(5)]


@pytest.fixture(scope="module")
def adata() -> AnnData:
    rng = np.random.default_rng(0)
    base = rng.uniform(-1, 4, N_GENES)
    rows, counts = [], []
    for donor in range(6):
        donor_effect = rng.normal(0, 0.3, N_GENES)
        for condition in ["ctrl", "stim"]:
            for cell_type in ["B cells", "T cells"]:
                mu = np.exp(base + donor_effect)
                if condition == "stim" and cell_type == "T cells":
                    mu[:5] *= 4
                counts.append(rng.negative_binomial(5, 5 / (5 + mu), (60, N_GENES)))
                rows += [(f"d{donor}", condition, cell_type, f"batch{donor % 2}")] * 60
    obs = pd.DataFrame(rows, columns=["donor", "label", "cell_type", "batch"])
    obs.index = obs.index.astype(str)
    return AnnData(np.vstack(counts).astype(np.float32), obs=obs)


def test_detects_effect_per_group(adata: AnnData) -> None:
    res = sc.tl.de.deseq2(
        adata,
        "label",
        reference="ctrl",
        sample_key="donor",
        groupby="cell_type",
        covariates=["donor"],
        n_jobs=1,
    )

    assert res.columns.tolist() == [
        "group",
        "condition",
        "reference",
        "gene",
        "base_mean",
        "log_fc",
        "log_fc_se",
        "statistic",
        "p_value",
        "adj_p_value",
    ]
    hits = res[res["adj_p_value"] < 0.01]
    assert set(hits.loc[hits["group"] == "T cells", "gene"]) == set(DE_GENES)
    assert (hits["group"] == "B cells").sum() == 0
    t_de = res[(res["group"] == "T cells") & res["gene"].isin(DE_GENES)]
    np.testing.assert_allclose(t_de["log_fc"], 2, atol=0.4)


def test_matches_pydeseq2(adata: AnnData) -> None:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    t_cells = adata[adata.obs["cell_type"] == "T cells"]
    res = sc.tl.de.deseq2(
        t_cells, "label", reference="ctrl", sample_key="donor", n_jobs=1
    )

    pb = sc.get.aggregate(t_cells, by=["donor", "label"], func="sum")
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=pd.DataFrame(
            pb.layers["sum"].astype(np.int64), index=pb.obs_names, columns=pb.var_names
        ),
        metadata=pd.DataFrame(
            {"condition": pb.obs["label"].astype(str).to_numpy()}, index=pb.obs_names
        ),
        design="~condition",
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds, contrast=["condition", "stim", "ctrl"], inference=inference, quiet=True
    )
    stats.summary()

    np.testing.assert_array_equal(res["gene"], stats.results_df.index)
    np.testing.assert_allclose(res["log_fc"], stats.results_df["log2FoldChange"])
    np.testing.assert_allclose(res["p_value"], stats.results_df["pvalue"])
    assert "group" not in res.columns


def test_min_cells_skips_groups(adata: AnnData) -> None:
    with pytest.raises(ValueError, match=r"No group had at least two samples"):
        sc.tl.de.deseq2(
            adata,
            "label",
            reference="ctrl",
            sample_key="donor",
            groupby="cell_type",
            min_cells=1000,
            n_jobs=1,
        )


@pytest.mark.parametrize(
    ("kwargs", "exc", "match"),
    [
        pytest.param(
            dict(sample_key="nope"), KeyError, r"Could not find", id="missing_key"
        ),
        pytest.param(
            dict(reference="nope"), ValueError, r"not a value", id="bad_reference"
        ),
        pytest.param(
            dict(covariates=["cell_type"]),
            ValueError,
            r"vary within samples",
            id="varying_covariate",
        ),
    ],
)
def test_invalid_input(
    adata: AnnData, kwargs: dict, exc: type[Exception], match: str
) -> None:
    kwargs = dict(reference="ctrl", sample_key="donor") | kwargs
    with pytest.raises(exc, match=match):
        sc.tl.de.deseq2(adata, "label", n_jobs=1, **kwargs)


def test_requires_counts(adata: AnnData) -> None:
    lognorm = adata.copy()
    sc.pp.log1p(lognorm)
    with pytest.raises(ValueError, match=r"needs raw counts"):
        sc.tl.de.deseq2(lognorm, "label", reference="ctrl", sample_key="donor")
