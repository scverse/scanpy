from string import ascii_lowercase
from time import sleep

from memory_profiler import memory_usage
import numpy as np
import pandas as pd
from scipy import sparse
import sys
import gc

import anndata
from anndata import AnnData


def get_actualsize(input_obj):
    """Using Python Garbage Collector to calculate the size of all elements attached to an object"""

    memory_size = 0
    ids = set()
    objects = [input_obj]
    while objects:
        new = []
        for obj in objects:
            if id(obj) not in ids:
                ids.add(id(obj))
                memory_size += sys.getsizeof(obj)
                new.append(obj)
        objects = gc.get_referents(*new)
    return memory_size


def get_anndata_memsize(adata):
    recording = memory_usage(
        (sedate(adata.copy, naplength=0.005), (adata,)), interval=0.001
    )
    diff = recording[-1] - recording[0]
    return diff


def get_peak_mem(op, interval=0.001):
    recording = memory_usage(op, interval=interval)
    return np.max(recording) - np.min(recording)


def sedate(func, naplength=0.05):
    """Make a function sleepy, so we can sample the start and end state."""

    def wrapped_function(*args, **kwargs):
        sleep(naplength)
        val = func(*args, **kwargs)
        sleep(naplength)
        return val

    return wrapped_function


# TODO: Factor out the time it takes to generate these


def gen_indexer(adata, dim, index_kind, ratio):
    dimnames = ("obs", "var")
    index_kinds = {"slice", "intarray", "boolarray", "strarray"}

    if index_kind not in index_kinds:
        raise ValueError(
            f"Argument 'index_kind' must be one of {index_kinds}. Was {index_kind}."
        )

    axis = dimnames.index(dim)
    subset = [slice(None), slice(None)]
    axis_size = adata.shape[axis]

    if index_kind == "slice":
        subset[axis] = slice(0, int(np.round(axis_size * ratio)))
    elif index_kind == "intarray":
        subset[axis] = np.random.choice(
            np.arange(axis_size), int(np.round(axis_size * ratio)), replace=False
        )
        subset[axis].sort()
    elif index_kind == "boolarray":
        pos = np.random.choice(
            np.arange(axis_size), int(np.round(axis_size * ratio)), replace=False
        )
        a = np.zeros(axis_size, dtype=bool)
        a[pos] = True
        subset[axis] = a
    elif index_kind == "strarray":
        subset[axis] = np.random.choice(
            getattr(adata, dim).index, int(np.round(axis_size * ratio)), replace=False
        )
    else:
        raise ValueError()
    return tuple(subset)


def take_view(adata, *, dim, index_kind, ratio=0.5, nviews=100):
    subset = gen_indexer(adata, dim, index_kind, ratio)
    views = []
    for i in range(nviews):
        views.append(adata[subset])


def take_repeated_view(adata, *, dim, index_kind, ratio=0.9, nviews=10):
    v = adata
    views = []
    for i in range(nviews):
        subset = gen_indexer(v, dim, index_kind, ratio)
        v = v[subset]
        views.append(v)


def gen_adata(n_obs, n_var, attr_set):
    if "X-csr" in attr_set:
        X = sparse.random(n_obs, n_var, density=0.1, format="csr")
    elif "X-dense" in attr_set:
        X = sparse.random(n_obs, n_var, density=0.1, format="csr")
        X = X.toarray()
    else:
        # TODO: Theres probably a better way to do this
        X = sparse.random(n_obs, n_var, density=0, format="csr")
    adata = AnnData(X)
    if "obs,var" in attr_set:
        adata.obs = pd.DataFrame(
            {k: np.random.randint(0, 100, n_obs) for k in ascii_lowercase},
            index=["cell{}".format(i) for i in range(n_obs)],
        )
        adata.var = pd.DataFrame(
            {k: np.random.randint(0, 100, n_var) for k in ascii_lowercase},
            index=["gene{}".format(i) for i in range(n_var)],
        )
    return adata
