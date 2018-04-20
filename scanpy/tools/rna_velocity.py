import numba
import numpy as np
from .. import logging as logg

def compute_velocity_graph(adata, adata_u, X_du):
    if (adata.shape[0] != adata_u.shape[0]
        or adata_u.shape[0] != X_du.shape[0]
        or X_du.shape[0] != adata.shape[0]):
        raise ValueError('Number of cells do not match.')

    from scanpy.neighbors import Neighbors, get_indices_distances_from_sparse_matrix
    neigh = Neighbors(adata)
    knn_indices, knn_distances = get_indices_distances_from_sparse_matrix(
        neigh.distances, neigh.n_neighbors)
    n_obs = adata.n_obs
    n_neighbors = neigh.n_neighbors

    from numpy.linalg import norm
    X_u = adata_u.X.toarray()
    X_du = X_du.astype('float32').toarray()

    def fill_graph():
        rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
        cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
        vals = np.zeros((n_obs * n_neighbors), dtype=np.float32)
        for i in range(n_obs):
            if i % 1000 == 0:
                logg.msg('{}/{},'.format(i, n_obs), end=' ', v=4)
            for nj in range(n_neighbors):
                j = knn_indices[i, nj]
                if j != i:
                    du_i = X_du[i]
                    du_ji = X_u[j] - X_u[i]
                    subset = np.logical_or(du_ji != 0, du_i != 0)
                    du_i = du_i[subset]
                    du_ji = du_ji[subset]
                    val = np.dot(du_ji, du_i) / norm(du_ji) / norm(du_i)
                    # if val > 0, this means transitioning from i to j,
                    # convention of standard stochastic matrices even though
                    # this isn't one
                    rows[i * n_neighbors + nj] = j
                    cols[i * n_neighbors + nj] = i
                    vals[i * n_neighbors + nj] = val
        return rows, cols, vals

    rows, cols, vals = fill_graph()
    from scipy.sparse import coo_matrix
    graph = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    graph.eliminate_zeros()
    return graph.tocsr()


def rna_velocity(adata, loomfile, copy=False):
    """Estimate RNA velocity [LaManno17]_

    This requires generating a loom file with Velocyto, which will store the
    counts of spliced, unspliced and ambiguous RNA for every cell and every
    gene.

    In contrast to Velocyto, here, we neither use RNA velocities for
    extrapolation nor for constructing a Markov process. Instead, we directly
    orient and weight edges in the nearest neighbor graph by computing
    ``cosine_similarity((x_i - x_j), v_i)``, where `i` labels a cell, `j` a
    neighbor of the cell, `x` a gene expression vector and `v` a velocity
    vector.
    """
    import loompy
    adata = adata.copy() if copy else adata

    # this is n_genes x n_cells
    ds = loompy.connect(loomfile)
    row_attrs = dict(ds.row_attrs.items())
    col_attrs = dict(ds.col_attrs.items())
    gene_names = [gene for gene in row_attrs['Gene'] if gene in adata.var_names]
    cell_names = [cell for cell in col_attrs['CellID'] if cell in adata.obs_names]

    # subset the spliced and unspliced matrices to the genes in adata
    from anndata.base import _normalize_index
    gene_index = _normalize_index(gene_names, adata.var_names)
    cell_index = _normalize_index(cell_names, adata.obs_names)
    if len(cell_index) == 0:
        raise ValueError(
            'Cell names in loom file do not match cell names in AnnData.')

    # subset to cells and genes present in adata
    X_s = ds.layer['spliced'][:, :].T[cell_index][:, gene_index].copy()
    X_u = ds.layer['unspliced'][:, :].T[cell_index][:, gene_index].copy()

    # for now, take non-normalized values
    from ..preprocessing.simple import normalize_per_cell
    X_s = normalize_per_cell(X_s, copy=True)
    X_u = normalize_per_cell(X_u, copy=True)

    # loop over genes
    import scipy.optimize as opt
    offset = np.zeros(X_s.shape[1], dtype='float32')
    gamma = np.zeros(X_s.shape[1], dtype='float32')
    for i in range(X_s.shape[1]):
        gamma[i], offset[i] = opt.minimize(
            lambda m: np.sum((-X_u[:, i] + X_s[:, i] * m[0] + m[1])**2),
            x0=(0.1, 1e-16),
            method='L-BFGS-B',
            bounds=[(1e-8, 30), (0, 1.5)]).x

    velocity = X_u.T - (gamma[:, None] * X_s.T + offset[:, None])
    velocity = velocity.T

    from ..neighbors import Neighbors, get_indices_distances_from_sparse_matrix
    neigh = Neighbors(adata)
    knn_indices, knn_distances = get_indices_distances_from_sparse_matrix(
        neigh.distances, neigh.n_neighbors)

    n_obs = adata.n_obs
    n_neighbors = neigh.n_neighbors

    from scipy.sparse import dok_matrix
    graph = dok_matrix((n_obs, n_obs), dtype='float32')
    from scipy.spatial.distance import cosine
    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] != i:
                val = 1 - cosine((X_s[i] - X_s[j]), velocity[i])
                if val > 0:
                    # transition from i to j
                    graph[j, i] = val
                else:
                    # transition from j to i
                    graph[i, j] = val

    graph = graph.tocoo().tocsr()
    adata.uns['rna_velocity'] = {}
    adata.uns['rna_velocity']['graph'] = graph
    adata.var['rna_velocity_gamma'] = gamma
    adata.var['rna_velocity_offset'] = offset
    return adata if copy else None
