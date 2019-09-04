from itertools import repeat, chain

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc


@pytest.fixture
def adata():
    return AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"]),
        var=pd.DataFrame({"gene_symbols": ["genesymbol1", "genesymbol2"]}, index=["gene1", "gene2"]),
        layers={"double": np.ones((2, 2)) * 2},
    )


def test_obs_df(adata):
    adata.obsm["eye"] = np.eye(2)
    adata.obsm["sparse"] = sparse.csr_matrix(np.eye(2))

    adata.raw = AnnData(
        X=np.zeros((2, 2)),
        var=pd.DataFrame({"gene_symbols": ["raw1", "raw2"]}, index=["gene1", "gene2"]),
    )
    assert np.all(np.equal(
        sc.get.obs_df(adata, keys=["gene2", "obs1"], obsm_keys=[("eye", 0), ("sparse", 1)]),
        pd.DataFrame({"gene2": [1, 1], "obs1": [0, 1], "eye-0": [1, 0], "sparse-1": [0, 1]}, index=adata.obs_names),
    ))
    assert np.all(np.equal(
        sc.get.obs_df(adata, keys=["genesymbol2", "obs1"], obsm_keys=[("eye", 0), ("sparse", 1)], gene_symbols="gene_symbols"),
        pd.DataFrame({"genesymbol2": [1, 1], "obs1": [0, 1], "eye-0": [1, 0], "sparse-1": [0, 1]}, index=adata.obs_names),
    ))
    assert np.all(np.equal(
        sc.get.obs_df(adata, keys=["gene2", "obs1"], layer="double"),
        pd.DataFrame({"gene2": [2, 2], "obs1": [0, 1]}, index=adata.obs_names),
    ))
    assert np.all(np.equal(
        sc.get.obs_df(adata, keys=["raw2", "obs1"], gene_symbols="gene_symbols", use_raw=True),
        pd.DataFrame({"raw2": [0, 0], "obs1": [0, 1]}, index=adata.obs_names),
    ))
    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        sc.get.obs_df(adata, keys=badkeys)
    with pytest.raises(AssertionError):
        sc.get.obs_df(adata, keys=["gene1"], use_raw=True, layer="double")
    assert all(badkey_err.match(k) for k in badkeys)


def test_var_df(adata):
    adata.varm["eye"] = np.eye(2)
    adata.varm["sparse"] = sparse.csr_matrix(np.eye(2))

    assert np.all(np.equal(
        sc.get.var_df(adata, keys=["cell2", "gene_symbols"], varm_keys=[("eye", 0), ("sparse", 1)]),
        pd.DataFrame({"cell2": [1, 1], "gene_symbols": ["genesymbol1", "genesymbol2"], "eye-0": [1, 0], "sparse-1": [0, 1]}, index=adata.obs_names),
    ))
    assert np.all(np.equal(
        sc.get.var_df(adata, keys=["cell1", "gene_symbols"], layer="double"),
        pd.DataFrame({"cell1": [2, 2], "gene_symbols": ["genesymbol1", "genesymbol2"]}, index=adata.obs_names),
    ))
    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        sc.get.var_df(adata, keys=badkeys)
    assert all(badkey_err.match(k) for k in badkeys)


def test_rank_genes_groups_df():
    a = np.zeros((20, 3))
    a[:10, 0] = 5
    adata = AnnData(
        a,
        obs=pd.DataFrame(
            {"celltype": list(chain(repeat("a", 10), repeat("b", 10)))},
            index=[f"cell{i}" for i in range(a.shape[0])],
        ),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(a.shape[1])]),
    )
    sc.tl.rank_genes_groups(adata, groupby="celltype", method="wilcoxon")
    dedf = sc.get.rank_genes_groups_df(adata, "a")
    assert dedf["pvals"].value_counts()[1.] == 2
    assert sc.get.rank_genes_groups_df(adata, "a", log2fc_max=.1).shape[0] == 2
    assert sc.get.rank_genes_groups_df(adata, "a", log2fc_min=.1).shape[0] == 1
    assert sc.get.rank_genes_groups_df(adata, "a", pval_cutoff=.9).shape[0] == 1
    del adata.uns["rank_genes_groups"]
    sc.tl.rank_genes_groups(
        adata, groupby="celltype", method="wilcoxon", key_added="different_key"
    )
    with pytest.raises(KeyError):
        sc.get.rank_genes_groups_df(adata, "a")
    dedf2 = sc.get.rank_genes_groups_df(adata, "a", key="different_key")
    pd.testing.assert_frame_equal(dedf, dedf2)
