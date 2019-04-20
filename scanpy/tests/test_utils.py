from types import ModuleType

import numpy as np
import pandas as pd

import scanpy as sc
from scanpy.utils import descend_classes_and_funcs, obs_values_df


def test_descend_classes_and_funcs():
    # create module hierarchy
    a = ModuleType('a')
    a.b = ModuleType('a.b')

    # populate with classes
    a.A = type('A', (), {})
    a.A.__module__ = a.__name__
    a.b.B = type('B', (), {})
    a.b.B.__module__ = a.b.__name__

    # create a loop to check if that gets caught
    a.b.a = a

    assert {a.A, a.b.B} == set(descend_classes_and_funcs(a, 'a'))


def test_obs_values():
    obs = pd.DataFrame({"obs1": [0, 1], "obs2": ["a", "b"]}, index=["cell1", "cell2"])
    var = pd.DataFrame(index=["gene1", "gene2"])
    adata = sc.AnnData(X=np.ones((2, 2)), obs=obs, var=var)
    test_df = pd.DataFrame({"gene2": [1, 1], "obs1": [0, 1]}, index=adata.obs_names)
    assert np.all(obs_values_df(adata, ["gene2", "obs1"]) == test_df)
