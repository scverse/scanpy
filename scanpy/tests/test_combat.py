import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

import scanpy as sc
from scanpy.preprocessing._combat import _standardize_data, _design_matrix


def test_norm():
    # this test trivially checks whether mean normalisation worked

    # load in data
    adata = sc.datasets.blobs()
    key = 'blobs'
    data = pd.DataFrame(data=adata.X.T, index=adata.var_names,
                        columns=adata.obs_names)

    # construct a pandas series of the batch annotation
    batch = pd.Series(adata.obs[key])
    model = pd.DataFrame({'batch': batch})

    # standardize the data
    s_data, design, var_pooled,  stand_mean = _standardize_data(model, data, 'batch')

    assert np.allclose(s_data.mean(axis = 1), np.zeros(s_data.shape[0]))


def test_covariates():
    adata = sc.datasets.blobs()
    key = 'blobs'

    X1 = sc.pp.combat(adata, key=key, inplace=False)

    np.random.seed(0)
    adata.obs['cat1'] = np.random.binomial(3, 0.5, size=(adata.n_obs))
    adata.obs['cat2'] = np.random.binomial(2, 0.1, size=(adata.n_obs))
    adata.obs['num1'] = np.random.normal(size=(adata.n_obs))

    X2 = sc.pp.combat(adata, key=key, covariates=['cat1', 'cat2', 'num1'], inplace=False)
    sc.pp.combat(adata, key=key, covariates=['cat1', 'cat2', 'num1'], inplace=True)

    assert X1.shape == X2.shape

    df = adata.obs[['cat1', 'cat2', 'num1', key]]
    batch_cats = adata.obs[key].cat.categories
    design = _design_matrix(df, key, batch_cats)

    assert len(design.columns) == 4 + len(batch_cats) - 1


def test_silhouette():
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
