from pathlib import Path

import anndata
import pytest
import numpy as np
from scipy import sparse

import scanpy as sc


HERE = Path(__file__).parent


@pytest.fixture
def adata_kbet_sim() -> anndata.AnnData:
    return anndata.read_h5ad(HERE / '_data' / 'kbet-sim.h5ad')


@pytest.fixture
def adj_r():
    return sparse.csr_matrix(np.loadtxt(HERE / '_data' / 'adj.tsv', delimiter='\t'))


def test_kbet_neighbors_match_r(adata_kbet_sim, adj_r):
    def nearest_neighbors(adj: sparse.spmatrix):
        def which_min(row: np.ndarray):
            nz = row != 0
            idx = np.arange(len(row))
            return idx[row == row[nz].min()][0]
        return np.fromiter(map(which_min, adj.A), int, adj.shape[0])
    sc.pp.pca(adata_kbet_sim)
    sc.pp.neighbors(adata_kbet_sim, n_neighbors=50)  # some higher number than necessary
    adj_py = adata_kbet_sim.uns['neighbors']['distances']
    # adj_py2 = adata_kbet_sim.uns['neighbors']['connectivities']

    nn_r = nearest_neighbors(adj_r)
    nn_py = nearest_neighbors(adj_py)

    assert (nn_r == nn_py).all()


def test_kbet_needs_neighbors(adata_kbet_sim):
    with pytest.raises(ValueError):
        sc.tl.kbet(adata_kbet_sim)


def test_kbet_basic(adata_kbet_sim):
    sc.pp.pca(adata_kbet_sim)
    # Heuristic gives: k=75
    sc.pp.neighbors(adata_kbet_sim, n_neighbors=75)
    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha)
    assert 1 - alpha < acceptance <= 1, adata_kbet_sim.obs['kbet']


def test_kbet_with_adj(adata_kbet_sim, adj_r):
    assert np.all((adj_r != 0).sum(axis=1) == 75)

    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha, adjacency=adj_r)
    assert 1 - alpha < acceptance <= 1, adata_kbet_sim.obs['kbet']


def test_kbet_heuristic(adata_kbet_sim):
    sc.pp.pca(adata_kbet_sim)
    sc.pp.kbet_neighbors(adata_kbet_sim)
    assert 65 < adata_kbet_sim.uns['neighbors']['params']['n_neighbors'] < 80
    alpha = .05
    acceptance, _ = sc.tl.kbet(adata_kbet_sim, alpha=alpha)
    assert 1 - alpha < acceptance <= 1, adata_kbet_sim.obs['kbet']
