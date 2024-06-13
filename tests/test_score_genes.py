from __future__ import annotations

import pickle
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

import scanpy as sc
from testing.scanpy._helpers.data import paul15

if TYPE_CHECKING:
    from typing import Literal

    from numpy.typing import NDArray


HERE = Path(__file__).parent / "_data"


def _create_random_gene_names(n_genes, name_length) -> NDArray[np.str_]:
    """
    creates a bunch of random gene names (just CAPS letters)
    """
    return np.array(
        [
            "".join(map(chr, np.random.randint(65, 90, name_length)))
            for _ in range(n_genes)
        ]
    )


def _create_sparse_nan_matrix(rows, cols, percent_zero, percent_nan):
    """
    creates a sparse matrix, with certain amounts of NaN and Zeros
    """
    A = np.random.randint(0, 1000, rows * cols).reshape((rows, cols)).astype("float32")
    maskzero = np.random.rand(rows, cols) < percent_zero
    masknan = np.random.rand(rows, cols) < percent_nan
    if np.any(maskzero):
        A[maskzero] = 0
    if np.any(masknan):
        A[masknan] = np.nan
    S = csr_matrix(A)
    return S


def _create_adata(n_obs, n_var, p_zero, p_nan):
    """
    creates an AnnData with random data, sparseness and some NaN values
    """
    X = _create_sparse_nan_matrix(n_obs, n_var, p_zero, p_nan)
    adata = AnnData(X)
    gene_names = _create_random_gene_names(n_var, name_length=6)
    adata.var_names = gene_names
    return adata


def test_score_with_reference():
    """
    Checks if score_genes output agrees with pre-computed reference values.
    The reference values had been generated using the same code
    and stored as a pickle object in ./data
    """

    adata = paul15()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    sc.pp.scale(adata)

    sc.tl.score_genes(adata, gene_list=adata.var_names[:100], score_name="Test")
    with (HERE / "score_genes_reference_paul2015.pkl").open("rb") as file:
        reference = pickle.load(file)
    # np.testing.assert_allclose(reference, adata.obs["Test"].to_numpy())
    np.testing.assert_array_equal(reference, adata.obs["Test"].to_numpy())


def test_add_score():
    """
    check the dtype of the scores
    check that non-existing genes get ignored
    """
    # TODO: write a test that costs less resources and is more meaningful
    adata = _create_adata(100, 1000, p_zero=0, p_nan=0)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

    # the actual genes names are all 6letters
    # create some non-estinsting names with 7 letters:
    non_existing_genes = _create_random_gene_names(n_genes=3, name_length=7)
    some_genes = np.r_[
        np.unique(np.random.choice(adata.var_names, 10)), np.unique(non_existing_genes)
    ]
    sc.tl.score_genes(adata, some_genes, score_name="Test")
    assert adata.obs["Test"].dtype == "float64"


def test_sparse_nanmean():
    """
    check that _sparse_nanmean() is equivalent to np.nanmean()
    """
    from scanpy.tools._score_genes import _sparse_nanmean

    R, C = 60, 50

    # sparse matrix, no NaN
    S = _create_sparse_nan_matrix(R, C, percent_zero=0.3, percent_nan=0)
    # col/col sum
    np.testing.assert_allclose(
        S.toarray().mean(0), np.array(_sparse_nanmean(S, 0)).flatten()
    )
    np.testing.assert_allclose(
        S.toarray().mean(1), np.array(_sparse_nanmean(S, 1)).flatten()
    )

    # sparse matrix with nan
    S = _create_sparse_nan_matrix(R, C, percent_zero=0.3, percent_nan=0.3)
    np.testing.assert_allclose(
        np.nanmean(S.toarray(), 1), np.array(_sparse_nanmean(S, 1)).flatten()
    )
    np.testing.assert_allclose(
        np.nanmean(S.toarray(), 0), np.array(_sparse_nanmean(S, 0)).flatten()
    )

    # edge case of only NaNs per row
    A = np.full((10, 1), np.nan)

    meanA = np.array(_sparse_nanmean(csr_matrix(A), 0)).flatten()
    np.testing.assert_allclose(np.nanmean(A, 0), meanA)


def test_sparse_nanmean_on_dense_matrix():
    """
    TypeError must be thrown when calling _sparse_nanmean with a dense matrix
    """
    from scanpy.tools._score_genes import _sparse_nanmean

    with pytest.raises(TypeError):
        _sparse_nanmean(np.random.rand(4, 5), 0)


def test_score_genes_sparse_vs_dense():
    """
    score_genes() should give the same result for dense and sparse matrices
    """
    adata_sparse = _create_adata(100, 1000, p_zero=0.3, p_nan=0.3)

    adata_dense = adata_sparse.copy()
    adata_dense.X = adata_dense.X.toarray()

    gene_set = adata_dense.var_names[:10]

    sc.tl.score_genes(adata_sparse, gene_list=gene_set, score_name="Test")
    sc.tl.score_genes(adata_dense, gene_list=gene_set, score_name="Test")

    np.testing.assert_allclose(
        adata_sparse.obs["Test"].values, adata_dense.obs["Test"].values
    )


def test_score_genes_deplete():
    """
    deplete some cells from a set of genes.
    their score should be <0 since the sum of markers is 0 and
    the sum of random genes is >=0

    check that for both sparse and dense matrices
    """
    adata_sparse = _create_adata(100, 1000, p_zero=0.3, p_nan=0.3)

    adata_dense = adata_sparse.copy()
    adata_dense.X = adata_dense.X.toarray()

    # here's an arbitary gene set
    gene_set = adata_dense.var_names[:10]

    for adata in [adata_sparse, adata_dense]:
        # deplete these genes in 50 cells,
        ix_obs = np.random.choice(adata.shape[0], 50)
        adata[ix_obs][:, gene_set].X = 0

        sc.tl.score_genes(adata, gene_list=gene_set, score_name="Test")
        scores = adata.obs["Test"].values

        np.testing.assert_array_less(scores[ix_obs], 0)


def test_npnanmean_vs_sparsemean(monkeypatch):
    """
    another check that _sparsemean behaves like np.nanmean!

    monkeypatch the _score_genes._sparse_nanmean function to np.nanmean
    and check that the result is the same as the non-patched (i.e. sparse_nanmean)
    function
    """

    adata = _create_adata(100, 1000, p_zero=0.3, p_nan=0.3)
    gene_set = adata.var_names[:10]

    # the unpatched, i.e. _sparse_nanmean version
    sc.tl.score_genes(adata, gene_list=gene_set, score_name="Test")
    sparse_scores = adata.obs["Test"].values.tolist()

    # now patch _sparse_nanmean by np.nanmean inside sc.tools
    def mock_fn(x: csr_matrix, axis: Literal[0, 1]):
        return np.nanmean(x.toarray(), axis, dtype="float64")

    monkeypatch.setattr(sc.tl._score_genes, "_sparse_nanmean", mock_fn)
    sc.tl.score_genes(adata, gene_list=gene_set, score_name="Test")
    dense_scores = adata.obs["Test"].values

    np.testing.assert_allclose(sparse_scores, dense_scores)


def test_missing_genes():
    adata = _create_adata(100, 1000, p_zero=0, p_nan=0)
    # These genes have a different length of name
    non_extant_genes = _create_random_gene_names(n_genes=3, name_length=7)

    with pytest.raises(ValueError):
        sc.tl.score_genes(adata, non_extant_genes)


def test_one_gene():
    # https://github.com/scverse/scanpy/issues/1395
    adata = _create_adata(100, 1000, p_zero=0, p_nan=0)
    sc.tl.score_genes(adata, [adata.var_names[0]])


def test_use_raw_None():
    adata = _create_adata(100, 1000, p_zero=0, p_nan=0)
    adata_raw = adata.copy()
    adata_raw.var_names = [str(i) for i in range(adata_raw.n_vars)]
    adata.raw = adata_raw

    sc.tl.score_genes(adata, adata_raw.var_names[:3], use_raw=None)


@pytest.mark.parametrize("gene_pool", [[], ["foo", "bar"]])
def test_invalid_gene_pool(gene_pool):
    adata = _create_adata(100, 1000, p_zero=0, p_nan=0)

    with pytest.raises(ValueError, match="reference set"):
        sc.tl.score_genes(adata, adata.var_names[:3], gene_pool=gene_pool)
