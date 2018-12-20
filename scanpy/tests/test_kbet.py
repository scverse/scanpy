from pathlib import Path

import anndata
import pytest

import scanpy as sc


HERE = Path(__file__).parent


@pytest.fixture
def adata_kbet_sim() -> anndata.AnnData:
    return anndata.read_h5ad(HERE / '_data' / 'kbet-sim.h5ad')


def test_kbet_needs_neighbors(adata_kbet_sim):
    with pytest.raises(ValueError):
        sc.tl.kbet(adata_kbet_sim)


def test_kbet_basic(adata_kbet_sim):
    sc.pp.pca(adata_kbet_sim)
    # Heuristic gives: k=75
    sc.pp.neighbors(adata_kbet_sim, n_neighbors=75)
    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha)
    assert acceptance > 1 - alpha


def test_kbet_with_adj(adata_kbet_sim):
    import numpy as np
    from scipy import sparse
    adj = sparse.csr_matrix(np.loadtxt(HERE / '_data' / 'adj.tsv', delimiter='\t'))
    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha, adjacency=adj)
    assert acceptance > 1 - alpha


def test_kbet_heuristic(adata_kbet_sim):
    sc.pp.pca(adata_kbet_sim)
    sc.pp.kbet_neighbors(adata_kbet_sim)
    assert 72 < adata_kbet_sim.uns['neighbors']['params']['n_neighbors'] < 78
    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha)
    assert acceptance > 1 - alpha
