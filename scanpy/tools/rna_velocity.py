

def rna_velocity(adata, loomfile):

    # this is n_genes x n_cells
    ds = loompy.connect(self.loom_filepath)
    X_spliced = ds.layer['spliced'][:, :]
    X_unspliced = ds.layer['unspliced'][:, :]
    # X_ambiguous = ds.layer['ambiguous'][:, :]

    # for now, take non-normalized values
    s = X_spliced
    u = X_unspliced

    # loop over genes
    q = np.zeros(s.shape[0], dtype='float32')
    gammas = np.zeros(s.shape[0], dtype='float32')
    R2 = np.zeros(s.shape[0], dtype='float32')
    for i in range(s.shape[0]):
        m = opt.minimize(
            lambda m: np.sum((-u[i] + s[i] * m[0] + m[1])**2), x0=(0.1, 1e-16),
            method='L-BFGS-B', bounds=[(1e-8, 30), (0, 1.5)]).x
        gammas[i] = m[0]
        q[i] = m[1]
        R2[i] = 1 - (np.sum((gammas[i] * s[i] + q[i] - u[i])**2)
                     / np.sum((u[i].mean() - u[i])**2))

    velocity = u - (gammas[:, None] * s + q[:, None])
