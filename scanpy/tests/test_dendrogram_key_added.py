import scanpy as sc
import numpy as np
import pytest

n_neighbors = 5
key = 'test'


@pytest.fixture
def adata():
    return sc.AnnData(sc.datasets.pbmc68k_reduced())


@pytest.mark.parametrize('groupby', ['bulk_labels', ['bulk_labels', 'phase']])
@pytest.mark.parametrize('key_added', [None, "custom_key"])
def test_dendrogram_key_added(adata, groupby, key_added):
    sc.tl.dendrogram(adata, groupby=groupby, key_added=key_added, use_rep="X_pca")
    if isinstance(groupby, list):
        dendrogram_key = f'dendrogram_{"_".join(groupby)}'
    else:
        dendrogram_key = f'dendrogram_{groupby}'

    if key_added is None:
        key_added = dendrogram_key
    assert key_added in adata.uns
