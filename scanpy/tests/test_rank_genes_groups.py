import pytest

from packaging import version
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
import scipy
from scipy import sparse as sp
from scipy.stats import mannwhitneyu
from numpy.random import negative_binomial, binomial, seed

import scanpy as sc
from anndata import AnnData
from scanpy.tools import rank_genes_groups
from scanpy.tools._rank_genes_groups import _RankGenes
from scanpy.get import rank_genes_groups_df
from scanpy.tests._data._cached_datasets import pbmc68k_reduced
from scanpy._utils import select_groups


HERE = Path(__file__).parent / Path('_data/')


# We test results for a simple generic example
# Tests are conducted for sparse and non-sparse AnnData objects.
# Due to minor changes in multiplication implementation for sparse and non-sparse objects,
# results differ (very) slightly


def get_example_data(*, sparse=False):
    # create test object
    adata = AnnData(
        np.multiply(binomial(1, 0.15, (100, 20)), negative_binomial(2, 0.25, (100, 20)))
    )
    # adapt marker_genes for cluster (so as to have some form of reasonable input
    adata.X[0:10, 0:5] = np.multiply(
        binomial(1, 0.9, (10, 5)), negative_binomial(1, 0.5, (10, 5))
    )

    # The following construction is inefficient, but makes sure that the same data is used in the sparse case
    if sparse:
        adata.X = sp.csr_matrix(adata.X)

    # Create cluster according to groups
    adata.obs['true_groups'] = pd.Categorical(
        np.concatenate(
            (
                np.zeros((10,), dtype=int),
                np.ones((90,), dtype=int),
            )
        )
    )

    return adata


def get_true_scores():
    with Path(HERE, 'objs_t_test.pkl').open('rb') as f:
        true_scores_t_test, true_names_t_test = pickle.load(f)
    with Path(HERE, 'objs_wilcoxon.pkl').open('rb') as f:
        true_scores_wilcoxon, true_names_wilcoxon = pickle.load(f)

    return (
        true_names_t_test,
        true_names_wilcoxon,
        true_scores_t_test,
        true_scores_wilcoxon,
    )


def test_results_dense():
    seed(1234)

    adata = get_example_data()
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


def test_results_sparse():
    seed(1234)

    adata = get_example_data(sparse=True)

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


def test_results_layers():
    seed(1234)

    adata = get_example_data(sparse=False)
    adata.layers["to_test"] = adata.X.copy()
    adata.X = adata.X * np.random.randint(0, 2, adata.shape, dtype=bool)

    (
        true_names_t_test,
        true_names_wilcoxon,
        true_scores_t_test,
        true_scores_wilcoxon,
    ) = get_true_scores()

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
