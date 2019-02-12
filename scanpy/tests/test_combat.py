import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

import scanpy as sc
from scanpy.preprocessing._combat import stand_data


def test_norm():
    # this test trivially checks wether mean normalisation worked
    
    # load in data
    adata = sc.datasets.blobs()
    key = 'blobs'
    data = pd.DataFrame(data=adata.X.T, index=adata.var_names,
                        columns=adata.obs_names)
    
    # construct a pandas series of the batch annotation
    batch = pd.Series(adata.obs[key])
    model = pd.DataFrame({'batch': batch})
    
    # standardize the data
    s_data, design, var_pooled,  stand_mean = stand_data(model, data) 

    assert np.allclose(s_data.mean(axis = 1), np.zeros(s_data.shape[0]))


def test_shilhouette():
    # this test checks wether combat can align data from several gaussians
    # it checks this by computing the silhouette coefficient in a pca embedding
    
    # load in data
    adata = sc.datasets.blobs()
    
    # apply combat
    sc.pp.combat(adata, 'blobs')
    
    # compute pca
    sc.tl.pca(adata)
    X_pca = adata.obsm['X_pca']
    
    # compute silhouette coefficient in pca
    sh = silhouette_score(X_pca[:, :2], adata.obs['blobs'].values )
    
    assert sh < 0.1
