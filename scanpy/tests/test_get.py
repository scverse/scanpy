from itertools import repeat, chain

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc


@pytest.fixture
def adata():
    """
    adata.X is np.ones((2, 2))
    adata.layers['double'] is sparse np.ones((2,2)) * 2 to also test sparse matrices
    """
    return AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame(
            {"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"]
        ),
        var=pd.DataFrame(
            {"gene_symbols": ["genesymbol1", "genesymbol2"]}, index=["gene1", "gene2"]
        ),
        layers={"double": sparse.csr_matrix(np.ones((2, 2)), dtype=int) * 2},
        dtype=int,
    )


def test_obs_df(adata):
    adata.obsm["eye"] = np.eye(2, dtype=int)
    adata.obsm["sparse"] = sparse.csr_matrix(np.eye(2), dtype='float64')

    # make raw with different genes than adata
    adata.raw = AnnData(
        X=np.array([[1, 2, 3], [2, 4, 6]]),
        var=pd.DataFrame(
            {"gene_symbols": ["raw1", "raw2", 'raw3']},
            index=["gene2", "gene3", "gene4"],
        ),
        dtype='float64',
    )
    pd.testing.assert_frame_equal(
        sc.get.obs_df(
            adata, keys=["gene2", "obs1"], obsm_keys=[("eye", 0), ("sparse", 1)]
        ),
        pd.DataFrame(
            {"gene2": [1, 1], "obs1": [0, 1], "eye-0": [1, 0], "sparse-1": [0.0, 1.0]},
            index=adata.obs_names,
        ),
    )
    pd.testing.assert_frame_equal(
        sc.get.obs_df(
            adata,
            keys=["genesymbol2", "obs1"],
            obsm_keys=[("eye", 0), ("sparse", 1)],
            gene_symbols="gene_symbols",
        ),
        pd.DataFrame(
            {
                "genesymbol2": [1, 1],
                "obs1": [0, 1],
                "eye-0": [1, 0],
                "sparse-1": [0.0, 1.0],
            },
            index=adata.obs_names,
        ),
    )
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=["gene2", "obs1"], layer="double"),
        pd.DataFrame({"gene2": [2, 2], "obs1": [0, 1]}, index=adata.obs_names),
    )

    pd.testing.assert_frame_equal(
        sc.get.obs_df(
            adata,
            keys=["raw2", "raw3", "obs1"],
            gene_symbols="gene_symbols",
            use_raw=True,
        ),
        pd.DataFrame(
            {"raw2": [2.0, 4.0], "raw3": [3.0, 6.0], "obs1": [0, 1]},
            index=adata.obs_names,
        ),
    )
    # test only obs
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=["obs1", "obs2"]),
        pd.DataFrame({"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"]),
    )
    # test only var
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=["gene1", "gene2"]),
        pd.DataFrame({"gene1": [1, 1], "gene2": [1, 1]}, index=adata.obs_names),
    )
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=["gene1", "gene2"]),
        pd.DataFrame({"gene1": [1, 1], "gene2": [1, 1]}, index=adata.obs_names),
    )

    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        sc.get.obs_df(adata, keys=badkeys)
    with pytest.raises(AssertionError):
        sc.get.obs_df(adata, keys=["gene1"], use_raw=True, layer="double")
    assert all(badkey_err.match(k) for k in badkeys)


def test_backed_vs_memory():
    "compares backed vs. memory"
    from pathlib import Path

    # get location test h5ad file in datasets
    HERE = Path(sc.__file__).parent
    adata_file = HERE / "datasets/10x_pbmc68k_reduced.h5ad"
    adata_backed = sc.read(adata_file, backed='r')
    adata = sc.read_h5ad(adata_file,)

    # use non-sequential list of genes
    genes = list(adata.var_names[20::-2])
    obs_names = ['bulk_labels', 'n_genes']
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=genes + obs_names),
        sc.get.obs_df(adata_backed, keys=genes + obs_names),
    )

    # use non-sequential list of cell indices
    cell_indices = list(adata.obs_names[30::-2])
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=cell_indices + ["highly_variable"]),
        sc.get.var_df(adata_backed, keys=cell_indices + ["highly_variable"]),
    )


def test_column_content():
    "uses a larger dataset to test column order and content"
    adata = sc.datasets.pbmc68k_reduced()

    # test that columns content is correct for obs_df
    query = ['CST3', 'NKG7', 'GNLY', 'louvain', 'n_counts', 'n_genes']
    df = sc.get.obs_df(adata, query)
    for col in query:
        assert col in df
        np.testing.assert_array_equal(query, df.columns)
        np.testing.assert_array_equal(df[col].values, adata.obs_vector(col))

    # test that columns content is correct for var_df
    cell_ids = list(adata.obs.sample(5).index)
    query = cell_ids + ['highly_variable', 'dispersions_norm', 'dispersions']
    df = sc.get.var_df(adata, query)
    np.testing.assert_array_equal(query, df.columns)
    for col in query:
        np.testing.assert_array_equal(df[col].values, adata.var_vector(col))


def test_var_df(adata):
    adata.varm["eye"] = np.eye(2, dtype=int)
    adata.varm["sparse"] = sparse.csr_matrix(np.eye(2), dtype='float64')

    pd.testing.assert_frame_equal(
        sc.get.var_df(
            adata,
            keys=["cell2", "gene_symbols"],
            varm_keys=[("eye", 0), ("sparse", 1)],
        ),
        pd.DataFrame(
            {
                "cell2": [1, 1],
                "gene_symbols": ["genesymbol1", "genesymbol2"],
                "eye-0": [1, 0],
                "sparse-1": [0.0, 1.0],
            },
            index=adata.var_names,
        ),
    )
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=["cell1", "gene_symbols"], layer="double"),
        pd.DataFrame(
            {"cell1": [2, 2], "gene_symbols": ["genesymbol1", "genesymbol2"]},
            index=adata.var_names,
        ),
    )
    # test only cells
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=["cell1", "cell2"]),
        pd.DataFrame({"cell1": [1, 1], "cell2": [1, 1]}, index=adata.var_names,),
    )
    # test only var columns
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=["gene_symbols"]),
        pd.DataFrame(
            {"gene_symbols": ["genesymbol1", "genesymbol2"]}, index=adata.var_names,
        ),
    )

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
    sc.tl.rank_genes_groups(adata, groupby="celltype", method="wilcoxon", pts=True)
    dedf = sc.get.rank_genes_groups_df(adata, "a")
    assert dedf["pvals"].value_counts()[1.0] == 2
    assert sc.get.rank_genes_groups_df(adata, "a", log2fc_max=0.1).shape[0] == 2
    assert sc.get.rank_genes_groups_df(adata, "a", log2fc_min=0.1).shape[0] == 1
    assert sc.get.rank_genes_groups_df(adata, "a", pval_cutoff=0.9).shape[0] == 1
    del adata.uns["rank_genes_groups"]
    sc.tl.rank_genes_groups(
        adata, groupby="celltype", method="wilcoxon", key_added="different_key", pts=True,
    )
    with pytest.raises(KeyError):
        sc.get.rank_genes_groups_df(adata, "a")
    dedf2 = sc.get.rank_genes_groups_df(adata, "a", key="different_key")
    pd.testing.assert_frame_equal(dedf, dedf2)
    assert 'pct_nz_group' in dedf2.columns
    assert 'pct_nz_reference' in dedf2.columns
    
    # get all groups
    dedf3 = sc.get.rank_genes_groups_df(adata, group=None, key="different_key")
    assert 'a' in dedf3['group'].unique()
    assert 'b' in dedf3['group'].unique()    
    adata.var_names.name = 'pr1388'
    sc.get.rank_genes_groups_df(adata, group=None, key="different_key")
