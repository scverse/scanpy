import numpy as np
from numpy import ma
from scipy import sparse as sp
# we don’t need this in requirements.txt, as it’s only needed for testing
from pytest import mark

from scanpy.data_structs.ann_data import AnnData, BoundStructArray, SMP_INDEX


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(ma.array([[1, 2], [3, 4]]), add={'mask': [0, 1, 1, 0]})
    AnnData(sp.eye(2))
    AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(Smp=['A', 'B']),
        dict(Feat=['a', 'b', 'c']))

    assert AnnData(np.array([1, 2])).X.shape == (2,)

    from pytest import raises
    raises(ValueError, AnnData,
           np.array([[1, 2], [3, 4]]),
           dict(TooLong=[1, 2, 3, 4]))


def test_ddata():
    ddata = dict(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        row_names=['A', 'B'],
        col_names=['a', 'b', 'c'])
    AnnData(ddata)


def test_names():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata.smp_names.tolist() == 'A B'.split()
    assert adata.var_names.tolist() == 'a b c'.split()

    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]),
                    var={'var_names': ['a', 'b']})
    assert adata.var_names.tolist() == ['a', 'b']


def test_indices_dtypes():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    # this assignment is nice
    adata.smp_names = ['d', 'b']
    from pytest import raises
    with raises(ValueError):
        # this is not possible currently as we store
        # datatypes of fixed length
        adata.smp_names = ['hello', 'b']
        # unicode not allowed for annotation
        adata.smp_names = ['ö', 'a']

def test_creation_from_vector():
    adata = AnnData(np.array([1, 2, 3]))
    adata = AnnData(np.array([[1], [2], [3]]))


def test_slicing():
    adata = AnnData(np.array([[1, 2, 3],
                              [4, 5, 6]]))

    assert np.all(adata[:, 0].X == adata.X[:, 0])

    assert adata[0, 0].X.tolist() == 1
    assert adata[0, :].X.tolist() == [1, 2, 3]
    assert adata[:, 0].X.tolist() == [1, 4]

    assert adata[:, [0, 1]].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array([0, 2])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, np.array([False, True, True])].X.tolist() == [[2, 3], [5, 6]]
    assert adata[:, 1:3].X.tolist() == [[2, 3], [5, 6]]


def test_slicing_strings():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata['A', 'a'].X.tolist() == 1
    assert adata['A', :].X.tolist() == [1, 2, 3]
    assert adata[:, 'a'].X.tolist() == [1, 4]
    assert adata[:, ['a', 'b']].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array(['a', 'c'])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, 'b':'c'].X.tolist() == [[2, 3], [5, 6]]

    from pytest import raises
    with raises(IndexError): _ = adata[:, 'X']
    with raises(IndexError): _ = adata['X', :]
    with raises(IndexError): _ = adata['A':'X', :]
    with raises(IndexError): _ = adata[:, 'a':'X']


def test_get_subset_add():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                    dict(Smp=['A', 'B']),
                    dict(Feat=['a', 'b', 'c']))

    assert adata[0, 0].smp['Smp'].tolist() == ['A']
    assert adata[0, 0].var['Feat'].tolist() == ['a']


def test_transpose():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(smp_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    adata1 = adata.T

    # make sure to not modify the original!
    assert adata.smp_names.tolist() == ['A', 'B']
    assert adata.var_names.tolist() == ['a', 'b', 'c']

    assert SMP_INDEX in adata1.smp.dtype.names
    assert adata1.smp_names.tolist() == ['a', 'b', 'c']
    assert adata1.var_names.tolist() == ['A', 'B']
    assert adata1.X.shape == adata.X.T.shape

    adata2 = adata.transpose()
    assert np.array_equal(adata1.X, adata2.X)
    assert np.array_equal(adata1.smp, adata2.smp)
    assert np.array_equal(adata1.var, adata2.var)

    assert adata1.smp._is_attr_of[1] == 'smp' == adata2.smp._is_attr_of[1]
    assert adata1.var._is_attr_of[1] == 'var' == adata2.var._is_attr_of[1]


def test_append_add_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.smp['new'] = [1, 2]
    adata.smp[['new2', 'new3']] = [['A', 'B'], ['c', 'd']]

    from pytest import raises
    with raises(ValueError):
        adata.smp['new4'] = 'far too long'.split()


def test_set_add():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.smp = {SMP_INDEX: ['1', '2']}
    assert isinstance(adata.smp, BoundStructArray)
    assert len(adata.smp.dtype) == 1

    adata.smp = dict(a=[3, 4])
    assert isinstance(adata.smp, BoundStructArray)
    assert len(adata.smp.dtype) == 2
    assert adata.smp_names.tolist() == ['1', '2']  # still the same smp_names

    from pytest import raises
    with raises(ValueError):
        adata.smp = dict(a=[1, 2, 3])


def test_print():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                    dict(foo=['A', 'B']),
                    dict(bar=['a', 'b', 'c']))
    print(adata)
    print('>>> print(adata.smp)')
    print(adata.smp)


def test_multicol_getitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    # 'a' and 'b' label columns
    adata.smp[['a', 'b']] = [[0, 1], [2, 3]]
    assert 0 in adata.smp['a']
    assert adata.smp['b'].tolist() == [1, 3]
    assert adata.smp[['a', 'b']].tolist() == [[0, 1], [2, 3]]


def test_multicol_single_key_setitem():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    # 'c' keeps the columns as should be
    adata.smp['c'] = np.array([[0, 1], [2, 3]])
    assert adata.smp.dtype.names == (SMP_INDEX, 'c001of002', 'c002of002')
    assert adata.smp.keys() == ['c']
    assert adata.smp['c'].tolist() == [[0, 1], [2, 3]]


def test_structdict_keys():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.smp['foo'] = np.array([[0, 1], [2, 3]])
    assert adata.smp_keys() == ['foo']
    assert adata.smp.keys() == ['foo']
    assert adata.smp.dtype.names == (SMP_INDEX, 'foo001of002', 'foo002of002')

    adata.smp['d'] = np.array([[0, 1], [2, 3]])
    assert adata.smp.keys() == ['foo', 'd']
    assert 'd' in adata.smp
    from pytest import raises
    with raises(KeyError):
        adata.smp['e']


def test_n_smps():
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    assert adata.n_smps == 3
    adata1 = adata[:2, ]
    assert adata1.n_smps == 2


def test_structdict_index():
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))

    assert adata.var.index.tolist() == ['0', '1']
    adata.var.index = ['1', '2']
    assert adata.var.index.tolist() == ['1', '2']
    adata.var_names = ['3', '4']
    assert adata.var.index.tolist() == ['3', '4']


def test_struct_dict_copy():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    scp = adata.smp.copy()

    assert adata.smp.__dict__.keys() == scp.__dict__.keys()


#TODO: remove logging and actually test values
from scanpy import logging as logg

def test_profile_memory():
    import gc
    dim = 10  # increase this when profiling
    print()
    logg.print_memory_usage('start profiling')
    X = np.random.rand(dim, dim).astype('float32')
    logg.print_memory_usage('allocated X')
    var_filter = np.array([0, 1])
    X = X[:, var_filter]
    logg.print_memory_usage('sliced X')
    X = np.random.rand(dim, dim).astype('float32')
    logg.print_memory_usage('allocated X')
    adata = AnnData(X)
    logg.print_memory_usage('init adata with reference to X')
    adata.var['multi'] = np.random.rand(dim, 3)
    logg.print_memory_usage('added some annotation')
    # ------------------------------------------------
    # compare adata.__getitem__ with adata.filter_var
    # ------------------------------------------------
    # here, it doesn't make a difference in other scenarios
    # (e.g. sc.preprocess.weinreb16), filter_var seems to invoke earlier garbage
    # collection than slicing
    # adata.filter_var(var_filter)  # inplace
    adata = adata[:, var_filter]  # with copy
    logg.print_memory_usage('sliced adata')
    gc.collect()
    logg.print_memory_usage('after calling gc.collect()')
    return adata


def test_profile_memory_2():
    adata = test_profile_memory()
    logg.print_memory_usage('after leaving function')
