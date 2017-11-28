import numpy as np
import pandas as pd
from anndata import AnnData

from scanpy.readwrite import read, write


def test_read_write_hdf5():
    adata = AnnData(
        data=np.array([[1, 0], [3, 0], [5, 6]]),
        smp={'row_names': ['name1', 'name2', 'name3'],
             'sanno1': ['cat1', 'cat2', 'cat2'],  # categorical anno
             'sanno2': ['s1', 's2', 's3'],  # string annotation
             'sanno3': [2.1, 2.2, 2.3]},  # float annotation
        var={'vanno1': [3.1, 3.2]},
        uns={'sanno1_colors': ['#000000', '#FFFFFF'],
             'uns2': ['some annotation']})
    assert pd.api.types.is_string_dtype(adata.smp['sanno1'])
    write('./test.h5', adata)
    adata = read('./test.h5')
    assert pd.api.types.is_categorical(adata.smp['sanno1'])
    assert pd.api.types.is_string_dtype(adata.smp['sanno2'])
    assert adata.smp.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.smp['sanno1'].cat.categories.tolist() == ['cat1', 'cat2']


def test_read_write_hdf5_sparse():
    from scipy.sparse import csr_matrix
    adata = AnnData(
        data=csr_matrix([[1, 0], [3, 0], [5, 6]]),
        smp={'row_names': ['name1', 'name2', 'name3'],
             'sanno1': ['cat1', 'cat2', 'cat2'],
             'sanno2': [2.1, 2.2, 2.3]},
        var={'vanno1': [3.1, 3.2]},
        uns={'sanno1_colors': ['#000000', '#FFFFFF'],
             'uns2_sparse': csr_matrix([[1, 0], [3, 0]])})
    assert pd.api.types.is_string_dtype(adata.smp['sanno1'])
    write('./test.h5', adata)
    adata = read('./test.h5')
    assert pd.api.types.is_categorical(adata.smp['sanno1'])
    assert adata.smp.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.smp['sanno1'].cat.categories.tolist() == ['cat1', 'cat2']


def test_compression():
    adata = AnnData(np.ones((10000, 1000)))
    write('./test_compr.h5', adata, compression='gzip')
    write('./test_no_compr.h5', adata)
