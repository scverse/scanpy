from pathlib import Path

import anndata
import pytest

import scanpy as sc


HERE = Path(__file__).parent


@pytest.fixture
def adata_kbet_sim():
    return anndata.read_h5ad(HERE / '_data' / 'kbet-sim.h5ad')


def test_kbet_needs_neighbors(adata_kbet_sim):
    with pytest.raises(ValueError):
        sc.tl.kbet(adata_kbet_sim)


def test_kbet_basic(adata_kbet_sim):
    sc.pp.pca(adata_kbet_sim)
    # Heuristic gives: k=75
    sc.pp.neighbors(adata_kbet_sim, n_neighbors=75)
    alpha = .05
    acceptance = sc.tl.kbet(adata_kbet_sim, alpha=alpha)
    assert acceptance > 1 - alpha
