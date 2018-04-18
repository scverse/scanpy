

def rna_velocity(adata, loomfile, copy=False):
    """Estimate RNA velocity. [LaManno17]_

    This requires generating a loom file with Velocyto, which stores the counts
    of spliced, unspliced and ambiguous RNA for every cell and every gene.

    In contrast to Velocyto, here, we neither use RNA velocities for
    extrapolation nor for constructing a Markov process. Instead, we directly
    orient and weight edges in the nearest neighbor graph by computing
    ``cosine_similarity((x_i - x_j), v_i)``, where `i` labels a cell, `j` a
    neighbor of the cell, `x` a gene expression vector and `v` a velocity
    vector.
    """
    adata = adata.copy() if copy else adata    

    # this is n_genes x n_cells
    ds = loompy.connect(self.loom_filepath)
    X_spliced = ds.layer['spliced'][:, :]
    X_unspliced = ds.layer['unspliced'][:, :]
    # X_ambiguous = ds.layer['ambiguous'][:, :]
    row_attrs = dict(ds.row_attrs.items())
    gene_names = row_attrs['Gene']

    # subset the spliced and unspliced matrices to the genes in adata
    from anndata.base import _normalize_index
    gene_index = _normalize_index(gene_names, adata.var_names)
    X_spliced = X_spliced[gene_index]
    X_unspliced = X_unspliced[gene_index]

    # for now, take non-normalized values
    from ..preprocessing.simple import normalize_per_cell
    normalize_per_cell(X_spliced.T)
    normalize_per_cell(X_unspliced.T)

    # loop over genes
    offset = np.zeros(s.shape[0], dtype='float32')
    gamma = np.zeros(s.shape[0], dtype='float32')
    for i in range(s.shape[0]):
        gamma[i], offset[i] = opt.minimize(
            lambda m: np.sum((-X_unspliced[i] + X_spliced[i] * m[0] + m[1])**2),
            x0=(0.1, 1e-16),
            method='L-BFGS-B',
            bounds=[(1e-8, 30), (0, 1.5)]).x

    velocity = X_unspliced - (gamma[:, None] * X_spliced + offset[:, None])

    from ..neighbors import Neighbors, get_indices_distances_from_sparse_matrix
    neigh = Neighbors(adata)
    knn_indices, knn_distances = get_indices_distances_from_sparse_matrix(
        neigh.distances, self.n_neighbors)
    
    n_obs = adata.n_obs
    n_neighbors = neigh.n_neighbors

    from scipy.sparse import dok_matrix
    graph = dok_matrix((n_obs, n_obs), dtype='float32')
    from scipy.spatial.distance import cosine
    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] != i:
                val = 1 - cosine((X_spliced[:, i] - X_spliced[:, j]), velocity[:, i])
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
    
