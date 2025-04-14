from __future__ import annotations

from functools import partial
from itertools import chain, repeat

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

import scanpy as sc
from scanpy.datasets._utils import filter_oldformatwarning
from testing.scanpy._helpers import anndata_v0_8_constructor_compat
from testing.scanpy._helpers.data import pbmc68k_reduced


# Override so warning gets caught
def transpose_adata(adata: AnnData, *, expect_duplicates: bool = False) -> AnnData:
    if not expect_duplicates:
        return adata.T
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        return adata.T


TRANSPOSE_PARAMS = pytest.mark.parametrize(
    ("dim", "transform", "func"),
    [
        ("obs", lambda x, expect_duplicates=False: x, sc.get.obs_df),
        ("var", transpose_adata, sc.get.var_df),
    ],
    ids=["obs_df", "var_df"],
)


@pytest.fixture
def adata() -> AnnData:
    """Create a tiny AnnData.

    `adata.X` is `np.ones((2, 2))`.
    `adata.layers['double']` is sparse `np.ones((2,2)) * 2` to also test sparse matrices.
    """
    return anndata_v0_8_constructor_compat(
        X=np.ones((2, 2), dtype=int),
        obs=pd.DataFrame(
            {"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"]
        ),
        var=pd.DataFrame(
            {"gene_symbols": ["genesymbol1", "genesymbol2"]}, index=["gene1", "gene2"]
        ),
        layers={"double": sparse.csr_matrix(np.ones((2, 2)), dtype=int) * 2},  # noqa: TID251
    )


########################
# obs_df, var_df tests #
########################


def test_obs_df(adata: AnnData):
    adata.obsm["eye"] = np.eye(2, dtype=int)
    adata.obsm["sparse"] = sparse.csr_matrix(np.eye(2), dtype="float64")  # noqa: TID251

    # make raw with different genes than adata
    adata.raw = anndata_v0_8_constructor_compat(
        X=np.array([[1, 2, 3], [2, 4, 6]], dtype=np.float64),
        var=pd.DataFrame(
            {"gene_symbols": ["raw1", "raw2", "raw3"]},
            index=["gene2", "gene3", "gene4"],
        ),
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
    # test handling of duplicated keys (in this case repeated gene names)
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=["gene1", "gene2", "gene1", "gene1"]),
        pd.DataFrame(
            {"gene1": [1, 1], "gene2": [1, 1]},
            index=adata.obs_names,
        )[["gene1", "gene2", "gene1", "gene1"]],
    )

    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        sc.get.obs_df(adata, keys=badkeys)
    with pytest.raises(AssertionError):
        sc.get.obs_df(adata, keys=["gene1"], use_raw=True, layer="double")
    assert all(badkey_err.match(k) for k in badkeys)

    # test non unique index
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata = sc.AnnData(
            np.arange(16).reshape(4, 4),
            obs=pd.DataFrame(index=["a", "a", "b", "c"]),
            var=pd.DataFrame(index=[f"gene{i}" for i in range(4)]),
        )
    df = sc.get.obs_df(adata, ["gene1"])
    pd.testing.assert_index_equal(df.index, adata.obs_names)


def test_repeated_gene_symbols():
    """Gene symbols column allows repeats, but we can't unambiguously get data for these values."""
    gene_symbols = [f"symbol_{i}" for i in ["a", "b", "b", "c"]]
    var_names = pd.Index([f"id_{i}" for i in ["a", "b.1", "b.2", "c"]])
    adata = sc.AnnData(
        np.arange(3 * 4, dtype=np.float32).reshape((3, 4)),
        var=pd.DataFrame({"gene_symbols": gene_symbols}, index=var_names),
    )

    with pytest.raises(KeyError, match="symbol_b"):
        sc.get.obs_df(adata, ["symbol_b"], gene_symbols="gene_symbols")

    expected = pd.DataFrame(
        np.arange(3 * 4).reshape((3, 4))[:, [0, 3]].astype(np.float32),
        index=adata.obs_names,
        columns=["symbol_a", "symbol_c"],
    )
    result = sc.get.obs_df(adata, ["symbol_a", "symbol_c"], gene_symbols="gene_symbols")

    pd.testing.assert_frame_equal(expected, result)


@filter_oldformatwarning
def test_backed_vs_memory():
    """Compares backed vs. memory."""
    from pathlib import Path

    # get location test h5ad file in datasets
    HERE = Path(sc.__file__).parent
    adata_file = HERE / "datasets/10x_pbmc68k_reduced.h5ad"
    adata_backed = sc.read(adata_file, backed="r")
    adata = sc.read_h5ad(adata_file)

    # use non-sequential list of genes
    genes = list(adata.var_names[20::-2])
    obs_names = ["bulk_labels", "n_genes"]
    pd.testing.assert_frame_equal(
        sc.get.obs_df(adata, keys=genes + obs_names),
        sc.get.obs_df(adata_backed, keys=genes + obs_names),
    )

    # use non-sequential list of cell indices
    cell_indices = list(adata.obs_names[30::-2])
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=[*cell_indices, "highly_variable"]),
        sc.get.var_df(adata_backed, keys=[*cell_indices, "highly_variable"]),
    )


def test_column_content():
    """Uses a larger dataset to test column order and content."""
    adata = pbmc68k_reduced()

    # test that columns content is correct for obs_df
    query = ["CST3", "NKG7", "GNLY", "louvain", "n_counts", "n_genes"]
    df = sc.get.obs_df(adata, query)
    for col in query:
        assert col in df
        np.testing.assert_array_equal(query, df.columns)
        np.testing.assert_array_equal(df[col].values, adata.obs_vector(col))

    # test that columns content is correct for var_df
    cell_ids = list(adata.obs.sample(5).index)
    query = [*cell_ids, "highly_variable", "dispersions_norm", "dispersions"]
    df = sc.get.var_df(adata, query)
    np.testing.assert_array_equal(query, df.columns)
    for col in query:
        np.testing.assert_array_equal(df[col].values, adata.var_vector(col))


def test_var_df(adata: AnnData):
    adata.varm["eye"] = np.eye(2, dtype=int)
    adata.varm["sparse"] = sparse.csr_matrix(np.eye(2), dtype="float64")  # noqa: TID251

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
        pd.DataFrame(
            {"cell1": [1, 1], "cell2": [1, 1]},
            index=adata.var_names,
        ),
    )
    # test only var columns
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=["gene_symbols"]),
        pd.DataFrame(
            {"gene_symbols": ["genesymbol1", "genesymbol2"]},
            index=adata.var_names,
        ),
    )

    # test handling of duplicated keys (in this case repeated cell names)
    pd.testing.assert_frame_equal(
        sc.get.var_df(adata, keys=["cell1", "cell2", "cell2", "cell1"]),
        pd.DataFrame(
            {"cell1": [1, 1], "cell2": [1, 1]},
            index=adata.var_names,
        )[["cell1", "cell2", "cell2", "cell1"]],
    )

    badkeys = ["badkey1", "badkey2"]
    with pytest.raises(KeyError) as badkey_err:
        sc.get.var_df(adata, keys=badkeys)
    assert all(badkey_err.match(k) for k in badkeys)


@TRANSPOSE_PARAMS
def test_just_mapping_keys(dim, transform, func):
    # https://github.com/scverse/scanpy/issues/1634
    # Test for error where just passing obsm_keys, but not keys, would cause error.
    mapping_attr = f"{dim}m"
    kwargs = {f"{mapping_attr}_keys": [("array", 0), ("array", 1)]}

    adata = transform(
        sc.AnnData(
            X=np.zeros((5, 5)),
            obsm={
                "array": np.arange(10).reshape((5, 2)),
            },
        )
    )

    expected = pd.DataFrame(
        np.arange(10).reshape((5, 2)),
        index=getattr(adata, f"{dim}_names"),
        columns=["array-0", "array-1"],
    )
    result = func(adata, **kwargs)

    pd.testing.assert_frame_equal(expected, result)


##################################
# Test errors for obs_df, var_df #
##################################


def test_non_unique_cols_value_error():
    M, N = 5, 3
    adata = sc.AnnData(
        X=np.zeros((M, N)),
        obs=pd.DataFrame(
            np.arange(M * 2).reshape((M, 2)),
            columns=["repeated_col", "repeated_col"],
            index=[f"cell_{i}" for i in range(M)],
        ),
        var=pd.DataFrame(
            index=[f"gene_{i}" for i in range(N)],
        ),
    )
    with pytest.raises(ValueError, match=r"adata\.obs contains duplicated columns"):
        sc.get.obs_df(adata, ["repeated_col"])


def test_non_unique_var_index_value_error():
    adata = sc.AnnData(
        X=np.ones((2, 3)),
        obs=pd.DataFrame(index=["cell-0", "cell-1"]),
        var=pd.DataFrame(index=["gene-0", "gene-0", "gene-1"]),
    )
    with pytest.raises(ValueError, match=r"adata\.var_names contains duplicated items"):
        sc.get.obs_df(adata, ["gene-0"])


def test_keys_in_both_obs_and_var_index_value_error():
    M, N = 5, 3
    adata = sc.AnnData(
        X=np.zeros((M, N)),
        obs=pd.DataFrame(
            np.arange(M),
            columns=["var_id"],
            index=[f"cell_{i}" for i in range(M)],
        ),
        var=pd.DataFrame(
            index=["var_id"] + [f"gene_{i}" for i in range(N - 1)],
        ),
    )
    with pytest.raises(KeyError, match="var_id"):
        sc.get.obs_df(adata, ["var_id"])


@TRANSPOSE_PARAMS
def test_repeated_cols(dim, transform, func):
    adata = transform(
        sc.AnnData(
            np.ones((5, 10)),
            obs=pd.DataFrame(
                np.ones((5, 2)), columns=["a_column_name", "a_column_name"]
            ),
            var=pd.DataFrame(index=[f"gene-{i}" for i in range(10)]),
        )
    )
    # (?s) is inline re.DOTALL
    with pytest.raises(ValueError, match=rf"(?s)^adata\.{dim}.*a_column_name.*$"):
        func(adata, ["gene_5"])


@TRANSPOSE_PARAMS
def test_repeated_index_vals(dim, transform, func):
    # This one could be reverted, see:
    # https://github.com/scverse/scanpy/pull/1583#issuecomment-770641710
    alt_dim = ["obs", "var"][dim == "obs"]

    adata = transform(
        sc.AnnData(
            np.ones((5, 10)),
            var=pd.DataFrame(
                index=["repeated_id"] * 2 + [f"gene-{i}" for i in range(8)]
            ),
        ),
        expect_duplicates=True,
    )

    with pytest.raises(
        ValueError,
        match=rf"(?s)adata\.{alt_dim}_names.*{alt_dim}_names_make_unique",
    ):
        func(adata, "gene_5")


@pytest.fixture(
    params=[
        "obs_df",
        "var_df",
        "obs_df:use_raw",
        "obs_df:gene_symbols",
        "obs_df:gene_symbols,use_raw",
    ]
)
def shared_key_adata(request):
    kind = request.param
    adata = sc.AnnData(
        np.arange(50).reshape((5, 10)),
        obs=pd.DataFrame(np.zeros((5, 1)), columns=["var_id"]),
        var=pd.DataFrame(index=["var_id"] + [f"gene_{i}" for i in range(1, 10)]),
    )
    if kind == "obs_df":
        return (
            adata,
            sc.get.obs_df,
            r"'var_id'.* adata\.obs .* adata.var_names",
        )
    elif kind == "var_df":
        return (
            adata.T,
            sc.get.var_df,
            r"'var_id'.* adata\.var .* adata.obs_names",
        )
    elif kind == "obs_df:use_raw":
        adata.raw = adata
        adata.var_names = [f"gene_{i}" for i in range(10)]
        return (
            adata,
            partial(sc.get.obs_df, use_raw=True),
            r"'var_id'.* adata\.obs .* adata\.raw\.var_names",
        )
    elif kind == "obs_df:gene_symbols":
        adata.var["gene_symbols"] = adata.var_names
        adata.var_names = [f"gene_{i}" for i in range(10)]
        return (
            adata,
            partial(sc.get.obs_df, gene_symbols="gene_symbols"),
            r"'var_id'.* adata\.obs .* adata\.var\['gene_symbols'\]",
        )
    elif kind == "obs_df:gene_symbols,use_raw":
        base = adata.copy()
        adata.var["gene_symbols"] = adata.var_names
        adata.var_names = [f"gene_{i}" for i in range(10)]
        base.raw = adata
        return (
            base,
            partial(
                sc.get.obs_df,
                gene_symbols="gene_symbols",
                use_raw=True,
            ),
            r"'var_id'.* adata\.obs .* adata\.raw\.var\['gene_symbols'\]",
        )
    else:
        pytest.fail("add branch for new kind")


def test_shared_key_errors(shared_key_adata):
    adata, func, regex = shared_key_adata

    # This should error
    with pytest.raises(KeyError, match=regex):
        func(adata, keys=["var_id"])

    # This shouldn't error
    _ = func(adata, keys=["gene_2"])


##############################
# rank_genes_groups_df tests #
##############################


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
        adata,
        groupby="celltype",
        method="wilcoxon",
        key_added="different_key",
        pts=True,
    )
    with pytest.raises(KeyError):
        sc.get.rank_genes_groups_df(adata, "a")
    dedf2 = sc.get.rank_genes_groups_df(adata, "a", key="different_key")
    pd.testing.assert_frame_equal(dedf, dedf2)
    assert "pct_nz_group" in dedf2.columns
    assert "pct_nz_reference" in dedf2.columns

    # get all groups
    dedf3 = sc.get.rank_genes_groups_df(adata, group=None, key="different_key")
    assert "a" in dedf3["group"].unique()
    assert "b" in dedf3["group"].unique()
    adata.var_names.name = "pr1388"
    sc.get.rank_genes_groups_df(adata, group=None, key="different_key")
