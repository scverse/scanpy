import numpy as np
from numpy.linalg import norm

from .. import logging as logg
from ..neighbors import Neighbors, get_indices_distances_from_sparse_matrix


def compute_velocity_graph(adata, adata_u, X_du):
    if (adata.shape[0] != adata_u.shape[0]
        or adata_u.shape[0] != X_du.shape[0]
        or X_du.shape[0] != adata.shape[0]):
        raise ValueError('Number of cells do not match.')

    neigh = Neighbors(adata)
    knn_indices, knn_distances = get_indices_distances_from_sparse_matrix(
        neigh.distances, neigh.n_neighbors)
    n_obs = adata.n_obs
    n_neighbors = neigh.n_neighbors

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
                    # dividing this by norm(du_i) doesn't make much of a difference
                    val = np.dot(du_ji, du_i) / norm(du_ji) / norm(du_i)
                    # if val > 0, this means transitioning from i to j,
                    # convention of standard stochastic matrices even though
                    # this isn't one
                    # the problem with velocities at the boundaries of the knn
                    # graph is that, no matter in which direction they point,
                    # they anticorrelate with all neighbors: hence, one will
                    # always observe "out-going" velocity even if there is no
                    # indication for that
                    if True:  # val > 0:
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
    # also see notebook "scanpy_rna_velocity_all"

    # the following only works in theory...
    
    # this is n_genes x n_cells
    ds = loompy.connect('./all_sgete_4GU75.loom')
    row_attrs = dict(ds.row_attrs.items())
    col_attrs = dict(ds.col_attrs.items())
    gene_names = [gene for gene in row_attrs['Gene'] if gene in adata.var_names]
    cell_names = [cell for cell in col_attrs['CellID'] if cell in adata.obs_names]

    # subset the s and u matrices to the genes in adata
    from anndata.base import _normalize_index
    gene_index = _normalize_index(gene_names, adata.var_names)
    cell_index = _normalize_index(cell_names, adata.obs_names)
    if len(cell_index) == 0:
        raise ValueError(
            'Cell names in loom file do not match cell names in AnnData.')
    # subset to cells and genes present in adata
    ad_s = AnnData(ds.layer['spliced'].sparse(gene_index, cell_index).tocsr().T)
    ad_u = AnnData(ds.layer['unspliced'].sparse(gene_index, cell_index).tocsr().T)
    ds.close()

    subset, _ = sc.pp.filter_genes(ad_u.X, min_cells=50)
    print(np.sum(subset))
    ad_s = ad_s[:, subset]
    ad_u = ad_u[:, subset]
    ad_s.var_names = np.array(gene_names)[subset]

    # loop over genes
    from scipy.sparse import dok_matrix
    offset = np.zeros(ad_s.shape[1], dtype='float32')
    gamma = np.zeros(ad_s.shape[1], dtype='float32')
    X_du = dok_matrix(ad_s.shape, dtype='float32')
    for i in range(ad_s.shape[1]):
        x = ad_s.X[:, i].toarray()
        y = ad_u.X[:, i].toarray()
        subset = np.logical_and(x > 0, y > 0)
        x = x[subset]
        y = y[subset]
        X = np.c_[np.ones(len(x)), x]
        offset[i], gamma[i] = np.linalg.pinv(X).dot(y)
        subset_indices = np.flatnonzero(subset)
        index = subset_indices, np.array([i for dummy in subset_indices])
        X_du[index] = y - gamma[i]*x - offset[i]
        # pl.scatter(x, y)
        # pl.scatter(x, gamma[i]*x + offset[i])
        # pl.scatter(x, X_du[index].toarray()[0])
        # pl.show()
    X_du = X_du.tocoo().tocsr()

    # sc.pp.neighbors(adata, n_neighbors=100)

    graph = compute_velocity_graph(adata, ad_u, X_du)

    return adata if copy else None
