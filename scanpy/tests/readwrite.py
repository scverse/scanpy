import numpy as np
import pandas as pd
from anndata import AnnData

from scanpy.readwrite import read, write

def test_read_write_hdf5():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(row_names=['abc', 'bcd'], scol1=['Aaa', 'Bbb'], scol2=[1.0, 2]),
        dict(vcol1=[3, 4, 4.1], vcol2=[1, 2, 3]))
    assert pd.api.types.is_string_dtype(adata.smp['scol1'])
    write('./test.h5', adata)
    adata = read('./test.h5')
    assert pd.api.types.is_categorical(adata.smp['scol1'])
    assert adata.smp.index.tolist() == ['abc', 'bcd']
    assert adata.smp['scol1'].cat.categories.tolist() == ['Aaa', 'Bbb']
