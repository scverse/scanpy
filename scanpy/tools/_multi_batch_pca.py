from ..preprocessing import pca
import numpy as np

def multi_batch_pca(adata, weights=None, batch_key='batch', **kwargs):
    batches = list(adata.obs[batch_key].unique())

    if weights is None:
        weights = np.ones(len(batches))

    means = np.zeros(adata.n_vars)
    scales = np.ones(adata.n_obs)

    for i, batch in enumerate(batches):
        batch_inds = adata.obs[batch_key] == batch
        adata_batch = adata[batch_inds]

        means += weights[i] * np.ravel(adata_batch.X.mean(0))
        scales[batch_inds] = np.sqrt(adata_batch.n_obs / weights[i])

    means /= sum(weights)

    X_proc = (adata.X - means)/scales[:, None]

    comps = pca(X_proc, return_info=True, zero_center=False, **kwargs)[1]

    adata.obsm['X_multi_batch_pca'] = adata.X.dot(comps.T)
