import pytest

import numpy as np
import pandas as pd
from anndata import AnnData

import scanpy as sc
import scanpy.external as sce

pytest.importorskip("phenograph")


def test_phenograph():
    df = np.random.rand(1000, 40)
    dframe = pd.DataFrame(df)
    dframe.index, dframe.columns = (map(str, dframe.index), map(str, dframe.columns))
    adata = AnnData(dframe)
    sc.tl.pca(adata, n_comps=20)
    sce.tl.phenograph(adata, clustering_algo="leiden", k=50)
    assert adata.obs['pheno_leiden'].shape[0], "phenograph_Community Detection Error!"
