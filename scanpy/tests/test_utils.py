from types import ModuleType
from itertools import repeat, chain

import numpy as np
import pandas as pd
import pytest

import scanpy as sc
from scanpy.utils import descend_classes_and_funcs, obs_values_df, rank_genes_groups_df


def test_descend_classes_and_funcs():
    # create module hierarchy
    a = ModuleType('a')
    a.b = ModuleType('a.b')

    # populate with classes
    a.A = type('A', (), {})
    a.A.__module__ = a.__name__
    a.b.B = type('B', (), {})
    a.b.B.__module__ = a.b.__name__

    # create a loop to check if that gets caught
    a.b.a = a

    assert {a.A, a.b.B} == set(descend_classes_and_funcs(a, 'a'))


def test_obs_values():
    obs = pd.DataFrame({"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"])
    var = pd.DataFrame({"gene_symbols": ["genesymbol1", "genesymbol2"]}, index=["gene1", "gene2"])
    obsm = {"eye": np.eye(2)}
    adata = sc.AnnData(X=np.ones((2, 2)), obs=obs, var=var, obsm=obsm)
    test_df = pd.DataFrame({"gene2": [1, 1], "obs1": [0, 1], "eye0": [1, 0]}, index=adata.obs_names)
    assert np.all(obs_values_df(adata, keys=["gene2", "obs1"], obsm_keys=[("eye", 0)]) == test_df)
    test_gene_name_df = pd.DataFrame({"genesymbol2": [1, 1], "obs1": [0, 1], "eye0": [1, 0]}, index=adata.obs_names)
    assert np.all(obs_values_df(adata, keys=["genesymbol2", "obs1"], obsm_keys=[("eye", 0)], gene_symbols="gene_symbols") == test_gene_name_df)
    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        obs_values_df(adata, keys=badkeys)
    assert all(badkey_err.match(k) for k in badkeys)


def test_rank_genes_groups_df():
    a = np.zeros((20, 3))
    a[:10, 0] = 5
    adata = sc.AnnData(
        a,
        obs=pd.DataFrame(
            {"celltype": list(chain(repeat("a", 10), repeat("b", 10)))},
            index=[f"cell{i}" for i in range(a.shape[0])]
        ),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(a.shape[1])]),
    )
    sc.tl.rank_genes_groups(adata, groupby="celltype", method="wilcoxon")
    dedf = rank_genes_groups_df(adata, "a")
    assert dedf["pvals"].value_counts()[1.] == 2
    assert rank_genes_groups_df(adata, "a", log2fc_max=.1).shape[0] == 2
    assert rank_genes_groups_df(adata, "a", log2fc_min=.1).shape[0] == 1
    assert rank_genes_groups_df(adata, "a", pval_cutoff=.9).shape[0] == 1