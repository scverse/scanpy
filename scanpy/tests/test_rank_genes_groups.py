from __future__ import annotations

import pickle
from pathlib import Path
from collections.abc import Callable
from typing import Any

import pytest
import numpy as np
import pandas as pd
import scipy
from numpy.typing import NDArray
from anndata import AnnData
from packaging import version
from scipy import sparse
from scipy.stats import mannwhitneyu
from anndata.tests.helpers import asarray
from numpy.random import negative_binomial, binomial, seed

import scanpy as sc
from scanpy.testing._helpers.data import pbmc68k_reduced
from scanpy.tools import rank_genes_groups
from scanpy.tools._rank_genes_groups import _RankGenes
from scanpy.get import rank_genes_groups_df
from scanpy._utils import select_groups
from scanpy._compat import DaskArray


HERE = Path(__file__).parent
DATA_PATH = HERE / '_data/'


# We test results for a simple generic example
# Tests are conducted for sparse and non-sparse AnnData objects.
# Due to minor changes in multiplication implementation for sparse and non-sparse objects,
# results differ (very) slightly


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
    adata.obs['true_groups'] = pd.Categorical(
        np.concatenate((np.zeros((10,), dtype=int), np.ones((90,), dtype=int)))
    )

    return adata


def get_true_scores() -> (
    tuple[
        NDArray[np.object_],
        NDArray[np.object_],
        NDArray[np.floating],
        NDArray[np.floating],
    ]
):
    with (DATA_PATH / 'objs_t_test.pkl').open('rb') as f:
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with (DATA_PATH / 'objs_wilcoxon.pkl').open('rb') as f:
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)

    return (
        true_names_t_test,
        true_names_wilcoxon,
        true_scores_t_test,
        true_scores_wilcoxon,
    )


@pytest.fixture(
    params=[
        pytest.param(asarray, id='numpy-ndarray'),
        pytest.param(sparse.csr_matrix, id='scipy-csr'),
        pytest.param(sparse.csc_matrix, id='scipy-csc'),
    ]
)
def array_type(request: pytest.FixtureRequest):
    """Overrided array_type without dask, for the time being"""
    return request.param


def test_results(array_type):
    seed(1234)

    adata = get_example_data(array_type)
    assert adata.raw is None  # Assumption for later checks

    (
        true_names_t_test,
        true_names_wilcoxon,
        true_scores_t_test,
        true_scores_wilcoxon,
    ) = get_true_scores()

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='t-test')

    adata.uns['rank_genes_groups']['names'] = adata.uns['rank_genes_groups'][
        'names'
    ].astype(true_names_t_test.dtype)

    for name in true_scores_t_test.dtype.names:
        assert np.allclose(
            true_scores_t_test[name], adata.uns['rank_genes_groups']['scores'][name]
        )
    assert np.array_equal(true_names_t_test, adata.uns['rank_genes_groups']['names'])
    assert adata.uns["rank_genes_groups"]["params"]["use_raw"] is False

    rank_genes_groups(adata, 'true_groups', n_genes=20, method='wilcoxon')

    adata.uns['rank_genes_groups']['names'] = adata.uns['rank_genes_groups'][
        'names'
    ].astype(true_names_wilcoxon.dtype)

    for name in true_scores_t_test.dtype.names:
        assert np.allclose(
            true_scores_wilcoxon[name][:7],
            adata.uns['rank_genes_groups']['scores'][name][:7],
        )
    assert np.array_equal(
        true_names_wilcoxon[:7], adata.uns['rank_genes_groups']['names'][:7]
    )
    assert adata.uns["rank_genes_groups"]["params"]["use_raw"] is False


def elem_mul(x: np.ndarray | sparse.spmatrix | DaskArray, y: np.ndarray):
    if isinstance(x, sparse.spmatrix):
        # returns coo_matrix, so cast back to input type
        return type(x)(x.multiply(y))
    else:
        return x * y


def test_results_layers(array_type):
    seed(1234)

    adata = get_example_data(array_type)
    adata.layers["to_test"] = adata.X.copy()
    adata.X = elem_mul(adata.X, np.random.randint(0, 2, adata.shape, dtype=bool))

    _, _, true_scores_t_test, true_scores_wilcoxon = get_true_scores()

    # Wilcoxon
    rank_genes_groups(
        adata,
        'true_groups',
        method='wilcoxon',
        layer="to_test",
        n_genes=20,
    )
    assert adata.uns["rank_genes_groups"]["params"]["use_raw"] is False
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(
            true_scores_wilcoxon[name][:7],
            adata.uns['rank_genes_groups']['scores'][name][:7],
        )

    rank_genes_groups(adata, 'true_groups', method='wilcoxon', n_genes=20)
    for name in true_scores_t_test.dtype.names:
        assert not np.allclose(
            true_scores_wilcoxon[name][:7],
            adata.uns['rank_genes_groups']['scores'][name][:7],
        )

    # t-test
    rank_genes_groups(
        adata,
        'true_groups',
        method='t-test',
        layer="to_test",
        use_raw=False,
        n_genes=20,
    )
    for name in true_scores_t_test.dtype.names:
        assert np.allclose(
            true_scores_t_test[name][:7],
            adata.uns['rank_genes_groups']['scores'][name][:7],
        )

    rank_genes_groups(adata, 'true_groups', method='t-test', n_genes=20)
    for name in true_scores_t_test.dtype.names:
        assert not np.allclose(
            true_scores_t_test[name][:7],
            adata.uns['rank_genes_groups']['scores'][name][:7],
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
        ValueError, match="Received `use_raw=True`, but `adata.raw` is empty"
    ):
        sc.tl.rank_genes_groups(pbmc, groupby="bulk_labels", use_raw=True)


def test_singlets():
    pbmc = pbmc68k_reduced()
    pbmc.obs['louvain'] = pbmc.obs['louvain'].cat.add_categories(['11'])
    pbmc.obs['louvain'][0] = '11'

    with pytest.raises(ValueError, match=rf"Could not calculate statistics.*{'11'}"):
        rank_genes_groups(pbmc, groupby='louvain')


def test_emptycat():
    pbmc = pbmc68k_reduced()
    pbmc.obs['louvain'] = pbmc.obs['louvain'].cat.add_categories(['11'])

    with pytest.raises(ValueError, match=rf"Could not calculate statistics.*{'11'}"):
        rank_genes_groups(pbmc, groupby='louvain')


def test_log1p_save_restore(tmp_path):
    """tests the sequence log1p→save→load→rank_genes_groups"""
    from anndata import read

    pbmc = pbmc68k_reduced()
    sc.pp.log1p(pbmc)

    path = tmp_path / 'test.h5ad'
    pbmc.write(path)

    pbmc = read(path)

    sc.tl.rank_genes_groups(pbmc, groupby='bulk_labels', use_raw=True)


def test_wilcoxon_symmetry():
    pbmc = pbmc68k_reduced()

    rank_genes_groups(
        pbmc,
        groupby="bulk_labels",
        groups=["CD14+ Monocyte", "Dendritic"],
        reference="Dendritic",
        method='wilcoxon',
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
        method='wilcoxon',
        rankby_abs=True,
    )

    stats_dend = (
        rank_genes_groups_df(pbmc, group="Dendritic").drop(columns="names").to_numpy()
    )

    assert np.allclose(np.abs(stats_mono), np.abs(stats_dend))


@pytest.mark.parametrize('reference', [True, False])
def test_wilcoxon_tie_correction(reference):
    pbmc = pbmc68k_reduced()

    groups = ['CD14+ Monocyte', 'Dendritic']
    groupby = 'bulk_labels'

    _, groups_masks = select_groups(pbmc, groups, groupby)

    X = pbmc.raw.X[groups_masks[0]].toarray()

    mask_rest = groups_masks[1] if reference else ~groups_masks[0]
    Y = pbmc.raw.X[mask_rest].toarray()

    # Handle scipy versions
    if version.parse(scipy.__version__) >= version.parse("1.7.0"):
        pvals = mannwhitneyu(X, Y, use_continuity=False, alternative='two-sided').pvalue
        pvals[np.isnan(pvals)] = 1.0
    else:
        # Backwards compat, to drop once we drop scipy < 1.7
        n_genes = X.shape[1]
        pvals = np.zeros(n_genes)

        for i in range(n_genes):
            try:
                _, pvals[i] = mannwhitneyu(
                    X[:, i], Y[:, i], use_continuity=False, alternative='two-sided'
                )
            except ValueError:
                pvals[i] = 1

    if reference:
        ref = groups[1]
    else:
        ref = 'rest'
        groups = groups[:1]

    test_obj = _RankGenes(pbmc, groups, groupby, reference=ref)
    test_obj.compute_statistics('wilcoxon', tie_correct=True)

    np.testing.assert_allclose(test_obj.stats[groups[0]]['pvals'], pvals)
