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


def comppute_arrows_embedding():
    if 'rna_velocity' not in adata.uns:
        raise ValueError('`arrows=True` requires `tl.rna_velocity` to be run before.')
    adjacency = adata.uns['rna_velocity']['graph']
    # loop over columns of adjacency, this is where transitions start from
    from numpy.linalg import norm
    V = np.zeros((adjacency.shape[0], 2), dtype='float32')
    for i, n in enumerate(adjacency.T):  # loop over columns (note the transposition)
        for j in n.nonzero()[1]:  # these are row indices
            diff = adata.obsm['X_' + basis][j] - adata.obsm['X_' + basis][i]
            # need the normalized difference vector: the distance in the embedding
            # might be completely meaningless
            diff /= norm(diff)
            V[i] += adjacency[j, i] * diff
    logg.hint('added \'V_{}\' to `.obsm`'.format(basis))
    adata.obsm['V_' + basis] = V
    X = adata.obsm['X_' + basis]
    for ax in axs:
        quiver_kwds = arrows_kwds if arrows_kwds is not None else {}
        ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], **quiver_kwds)    


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
    # see notebook "scanpy_rna_velocity_all"
    return adata if copy else None
