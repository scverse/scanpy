from __future__ import annotations

import gc
import sys
from functools import wraps
from string import ascii_lowercase
from time import sleep
from typing import TYPE_CHECKING, Literal, TypeVar, get_args

import numpy as np
import pandas as pd
from anndata import AnnData
from memory_profiler import memory_usage
from scipy import sparse

if TYPE_CHECKING:
    from collections.abc import Callable, Set

    from numpy.typing import NDArray

    OneDIndex = NDArray[np.bool_] | NDArray[np.intp] | slice
    C = TypeVar("C", bound=Callable)


def get_actualsize(input_obj: object) -> int:
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


def get_anndata_memsize(adata: AnnData) -> float:
    recording = memory_usage(
        (sedate(adata.copy, seconds=0.005), (adata,)), interval=0.001
    )
    diff = recording[-1] - recording[0]
    return diff


def get_peak_mem(op, interval: float = 0.001) -> float:
    recording = memory_usage(op, interval=interval)
    return np.max(recording) - np.min(recording)


def sedate(func: C, *, seconds: float = 0.05) -> C:
    """Make a function sleepy, so we can sample the start and end state."""

    @wraps(func)
    def wrapped_function(*args, **kwargs):
        sleep(seconds)
        val = func(*args, **kwargs)
        sleep(seconds)
        return val

    return wrapped_function  # type: ignore


# TODO: Factor out the time it takes to generate these

IndexKind = Literal["slice", "intarray", "boolarray", "strarray"]
INDEX_KINDS: frozenset[IndexKind] = frozenset(get_args(IndexKind))


def gen_indexer(
    adata: AnnData,
    dim_name: Literal["obs", "var"],
    *,
    index_kind: IndexKind,
    ratio: float,
) -> tuple[OneDIndex, OneDIndex]:
    if index_kind not in INDEX_KINDS:
        raise ValueError(
            f"Argument 'index_kind' must be one of {INDEX_KINDS}. Was {index_kind}."
        )

    subset: dict[Literal["obs", "var"], OneDIndex] = dict(
        obs=slice(None), var=slice(None)
    )
    axis_size = adata.shape[("obs", "var").index(dim_name)]

    if index_kind == "slice":
        subset[dim_name] = slice(0, int(np.round(axis_size * ratio)))
    elif index_kind == "intarray":
        subset[dim_name] = np.random.choice(
            np.arange(axis_size), int(np.round(axis_size * ratio)), replace=False
        )
        subset[dim_name].sort()
    elif index_kind == "boolarray":
        pos = np.random.choice(
            np.arange(axis_size), int(np.round(axis_size * ratio)), replace=False
        )
        a = np.zeros(axis_size, dtype=bool)
        a[pos] = True
        subset[dim_name] = a
    elif index_kind == "strarray":
        dim: pd.DataFrame = getattr(adata, dim_name)
        subset[dim_name] = np.random.choice(
            dim.index, int(np.round(axis_size * ratio)), replace=False
        )
    else:
        raise ValueError()
    return (subset["obs"], subset["var"])


def take_view(
    adata: AnnData,
    dim_name: Literal["obs", "var"],
    *,
    index_kind: IndexKind,
    ratio: float = 0.5,
    nviews: int = 100,
):
    subset = gen_indexer(adata, dim_name, index_kind=index_kind, ratio=ratio)
    views = []
    for _ in range(nviews):
        views.append(adata[subset])


def take_repeated_view(
    adata: AnnData,
    dim_name: Literal["obs", "var"],
    *,
    index_kind: IndexKind,
    ratio: float = 0.9,
    nviews: int = 10,
):
    v = adata
    views = []
    for _ in range(nviews):
        subset = gen_indexer(v, dim_name, index_kind=index_kind, ratio=ratio)
        v = v[subset]
        views.append(v)


def gen_adata(
    n_obs: int, n_var: int, attr_set: Set[Literal["X-csr", "X-dense", "obs,var"]]
) -> AnnData:
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
            index=[f"cell{i}" for i in range(n_obs)],
        )
        adata.var = pd.DataFrame(
            {k: np.random.randint(0, 100, n_var) for k in ascii_lowercase},
            index=[f"gene{i}" for i in range(n_var)],
        )
    return adata
