"""Utility functions and classes

This file largely consists of the old _utils.py file. Over time, these functions
should be moved of this file.
"""

from __future__ import annotations

import importlib.util
import inspect
import random
import re
import sys
import warnings
from collections import namedtuple
from contextlib import contextmanager
from enum import Enum
from functools import partial, singledispatch, wraps
from operator import mul, truediv
from textwrap import dedent
from types import MethodType, ModuleType
from typing import TYPE_CHECKING, overload
from weakref import WeakSet

import h5py
import numpy as np
from anndata import __version__ as anndata_version
from packaging.version import Version
from scipy import sparse
from sklearn.utils import check_random_state

from .. import logging as logg
from .._compat import DaskArray
from .._settings import settings
from .compute.is_constant import is_constant  # noqa: F401

if Version(anndata_version) >= Version("0.10.0"):
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )
else:
    from anndata._core.sparse_dataset import SparseDataset

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path
    from typing import Any, Callable, Literal, TypeVar, Union

    from anndata import AnnData
    from numpy.typing import DTypeLike, NDArray

    # e.g. https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    # maybe in the future random.Generator
    AnyRandom = Union[int, np.random.RandomState, None]


class Empty(Enum):
    token = 0

    def __repr__(self) -> str:
        return "_empty"


_empty = Empty.token


class RNGIgraph:
    """
    Random number generator for ipgraph so global seed is not changed.
    See :func:`igraph.set_random_number_generator` for the requirements.
    """

    def __init__(self, random_state: int = 0) -> None:
        self._rng = check_random_state(random_state)

    def __getattr__(self, attr: str):
        return getattr(self._rng, "normal" if attr == "gauss" else attr)


def ensure_igraph() -> None:
    if importlib.util.find_spec("igraph"):
        return
    raise ImportError(
        "Please install the igraph package: "
        "`conda install -c conda-forge python-igraph` or "
        "`pip3 install igraph`."
    )


@contextmanager
def set_igraph_random_state(random_state: int):
    ensure_igraph()
    import igraph

    rng = RNGIgraph(random_state)
    try:
        igraph.set_random_number_generator(rng)
        yield None
    finally:
        igraph.set_random_number_generator(random)


EPS = 1e-15


def check_versions():
    if Version(anndata_version) < Version("0.6.10"):
        from .. import __version__

        raise ImportError(
            f"Scanpy {__version__} needs anndata version >=0.6.10, "
            f"not {anndata_version}.\nRun `pip install anndata -U --no-deps`."
        )


def getdoc(c_or_f: Callable | type) -> str | None:
    if getattr(c_or_f, "__doc__", None) is None:
        return None
    doc = inspect.getdoc(c_or_f)
    if isinstance(c_or_f, type) and hasattr(c_or_f, "__init__"):
        sig = inspect.signature(c_or_f.__init__)
    else:
        sig = inspect.signature(c_or_f)

    def type_doc(name: str):
        param: inspect.Parameter = sig.parameters[name]
        cls = getattr(param.annotation, "__qualname__", repr(param.annotation))
        if param.default is not param.empty:
            return f"{cls}, optional (default: {param.default!r})"
        else:
            return cls

    return "\n".join(
        f"{line} : {type_doc(line)}" if line.strip() in sig.parameters else line
        for line in doc.split("\n")
    )


def renamed_arg(old_name, new_name, *, pos_0: bool = False):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if old_name in kwargs:
                f_name = func.__name__
                pos_str = (
                    (
                        f" at first position. Call it as `{f_name}(val, ...)` "
                        f"instead of `{f_name}({old_name}=val, ...)`"
                    )
                    if pos_0
                    else ""
                )
                msg = (
                    f"In function `{f_name}`, argument `{old_name}` "
                    f"was renamed to `{new_name}`{pos_str}."
                )
                warnings.warn(msg, FutureWarning, stacklevel=3)
                if pos_0:
                    args = (kwargs.pop(old_name), *args)
                else:
                    kwargs[new_name] = kwargs.pop(old_name)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def _import_name(name: str) -> Any:
    from importlib import import_module

    parts = name.split(".")
    obj = import_module(parts[0])
    for i, name in enumerate(parts[1:]):
        try:
            obj = import_module(f"{obj.__name__}.{name}")
        except ModuleNotFoundError:
            break
    else:
        i = len(parts)
    for name in parts[i + 1 :]:
        try:
            obj = getattr(obj, name)
        except AttributeError:
            raise RuntimeError(f"{parts[:i]}, {parts[i + 1:]}, {obj} {name}")
    return obj


def _one_of_ours(obj, root: str):
    return (
        hasattr(obj, "__name__")
        and not obj.__name__.split(".")[-1].startswith("_")
        and getattr(
            obj, "__module__", getattr(obj, "__qualname__", obj.__name__)
        ).startswith(root)
    )


def descend_classes_and_funcs(mod: ModuleType, root: str, encountered=None):
    if encountered is None:
        encountered = WeakSet()
    for obj in vars(mod).values():
        if not _one_of_ours(obj, root) or obj in encountered:
            continue
        encountered.add(obj)
        if callable(obj) and not isinstance(obj, MethodType):
            yield obj
            if isinstance(obj, type):
                for m in vars(obj).values():
                    if callable(m) and _one_of_ours(m, root):
                        yield m
        elif isinstance(obj, ModuleType):
            if obj.__name__.startswith("scanpy.tests"):
                # Python’s import mechanism seems to add this to `scanpy`’s attributes
                continue
            yield from descend_classes_and_funcs(obj, root, encountered)


def annotate_doc_types(mod: ModuleType, root: str):
    for c_or_f in descend_classes_and_funcs(mod, root):
        c_or_f.getdoc = partial(getdoc, c_or_f)


def _doc_params(**kwds):
    """\
    Docstrings should start with ``\\`` in the first line for proper formatting.
    """

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj

    return dec


def _check_array_function_arguments(**kwargs):
    """Checks for invalid arguments when an array is passed.

    Helper for functions that work on either AnnData objects or array-likes.
    """
    # TODO: Figure out a better solution for documenting dispatched functions
    invalid_args = [k for k, v in kwargs.items() if v is not None]
    if len(invalid_args) > 0:
        raise TypeError(
            f"Arguments {invalid_args} are only valid if an AnnData object is passed."
        )


def _check_use_raw(adata: AnnData, use_raw: None | bool) -> bool:
    """
    Normalize checking `use_raw`.

    My intentention here is to also provide a single place to throw a deprecation warning from in future.
    """
    if use_raw is not None:
        return use_raw
    return adata.raw is not None


# --------------------------------------------------------------------------------
# Graph stuff
# --------------------------------------------------------------------------------


def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es["weight"] = weights
    except KeyError:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f"The constructed graph has only {g.vcount()} nodes. "
            "Your adjacency matrix contained redundant nodes."
        )
    return g


# --------------------------------------------------------------------------------
# Group stuff
# --------------------------------------------------------------------------------


def compute_association_matrix_of_groups(
    adata: AnnData,
    prediction: str,
    reference: str,
    *,
    normalization: Literal["prediction", "reference"] = "prediction",
    threshold: float = 0.01,
    max_n_names: int | None = 2,
):
    """Compute overlaps between groups.

    See ``identify_groups`` for identifying the groups.

    Parameters
    ----------
    adata
    prediction
        Field name of adata.obs.
    reference
        Field name of adata.obs.
    normalization
        Whether to normalize with respect to the predicted groups or the
        reference groups.
    threshold
        Do not consider associations whose overlap is below this fraction.
    max_n_names
        Control how many reference names you want to be associated with per
        predicted name. Set to `None`, if you want all.

    Returns
    -------
    asso_names
        List of associated reference names
        (`max_n_names` for each predicted name).
    asso_matrix
        Matrix where rows correspond to the predicted labels and columns to the
        reference labels, entries are proportional to degree of association.
    """
    if normalization not in {"prediction", "reference"}:
        raise ValueError(
            '`normalization` needs to be either "prediction" or "reference".'
        )
    sanitize_anndata(adata)
    cats = adata.obs[reference].cat.categories
    for cat in cats:
        if cat in settings.categories_to_ignore:
            logg.info(
                f"Ignoring category {cat!r} "
                "as it’s in `settings.categories_to_ignore`."
            )
    asso_names = []
    asso_matrix = []
    for ipred_group, pred_group in enumerate(adata.obs[prediction].cat.categories):
        if "?" in pred_group:
            pred_group = str(ipred_group)
        # starting from numpy version 1.13, subtractions of boolean arrays are deprecated
        mask_pred = adata.obs[prediction].values == pred_group
        mask_pred_int = mask_pred.astype(np.int8)
        asso_matrix += [[]]
        for ref_group in adata.obs[reference].cat.categories:
            mask_ref = (adata.obs[reference].values == ref_group).astype(np.int8)
            mask_ref_or_pred = mask_ref.copy()
            mask_ref_or_pred[mask_pred] = 1
            # e.g. if the pred group is contained in mask_ref, mask_ref and
            # mask_ref_or_pred are the same
            if normalization == "prediction":
                # compute which fraction of the predicted group is contained in
                # the ref group
                ratio_contained = (
                    np.sum(mask_pred_int) - np.sum(mask_ref_or_pred - mask_ref)
                ) / np.sum(mask_pred_int)
            else:
                # compute which fraction of the reference group is contained in
                # the predicted group
                ratio_contained = (
                    np.sum(mask_ref) - np.sum(mask_ref_or_pred - mask_pred_int)
                ) / np.sum(mask_ref)
            asso_matrix[-1] += [ratio_contained]
        name_list_pred = [
            cats[i] if cats[i] not in settings.categories_to_ignore else ""
            for i in np.argsort(asso_matrix[-1])[::-1]
            if asso_matrix[-1][i] > threshold
        ]
        asso_names += ["\n".join(name_list_pred[:max_n_names])]
    Result = namedtuple(
        "compute_association_matrix_of_groups", ["asso_names", "asso_matrix"]
    )
    return Result(asso_names=asso_names, asso_matrix=np.array(asso_matrix))


def get_associated_colors_of_groups(reference_colors, asso_matrix):
    return [
        {
            reference_colors[i_ref]: asso_matrix[i_pred, i_ref]
            for i_ref in range(asso_matrix.shape[1])
        }
        for i_pred in range(asso_matrix.shape[0])
    ]


def identify_groups(ref_labels, pred_labels, return_overlaps=False):
    """Which predicted label explains which reference label?

    A predicted label explains the reference label which maximizes the minimum
    of ``relative_overlaps_pred`` and ``relative_overlaps_ref``.

    Compare this with ``compute_association_matrix_of_groups``.

    Returns
    -------
    A dictionary of length ``len(np.unique(ref_labels))`` that stores for each
    reference label the predicted label that best explains it.

    If ``return_overlaps`` is ``True``, this will in addition return the overlap
    of the reference group with the predicted group; normalized with respect to
    the reference group size and the predicted group size, respectively.
    """
    ref_unique, ref_counts = np.unique(ref_labels, return_counts=True)
    ref_dict = dict(zip(ref_unique, ref_counts))
    pred_unique, pred_counts = np.unique(pred_labels, return_counts=True)
    pred_dict = dict(zip(pred_unique, pred_counts))
    associated_predictions = {}
    associated_overlaps = {}
    for ref_label in ref_unique:
        sub_pred_unique, sub_pred_counts = np.unique(
            pred_labels[ref_label == ref_labels], return_counts=True
        )
        relative_overlaps_pred = [
            sub_pred_counts[i] / pred_dict[n] for i, n in enumerate(sub_pred_unique)
        ]
        relative_overlaps_ref = [
            sub_pred_counts[i] / ref_dict[ref_label]
            for i, n in enumerate(sub_pred_unique)
        ]
        relative_overlaps = np.c_[relative_overlaps_pred, relative_overlaps_ref]
        relative_overlaps_min = np.min(relative_overlaps, axis=1)
        pred_best_index = np.argsort(relative_overlaps_min)[::-1]
        associated_predictions[ref_label] = sub_pred_unique[pred_best_index]
        associated_overlaps[ref_label] = relative_overlaps[pred_best_index]
    if return_overlaps:
        return associated_predictions, associated_overlaps
    else:
        return associated_predictions


# --------------------------------------------------------------------------------
# Other stuff
# --------------------------------------------------------------------------------


# backwards compat... remove this in the future
def sanitize_anndata(adata: AnnData) -> None:
    """Transform string annotations to categoricals."""
    adata._sanitize()


def view_to_actual(adata: AnnData) -> None:
    if adata.is_view:
        warnings.warn(
            "Received a view of an AnnData. Making a copy.",
            stacklevel=2,
        )
        adata._init_as_actual(adata.copy())


def moving_average(a: np.ndarray, n: int):
    """Moving average over one-dimensional array.

    Parameters
    ----------
    a
        One-dimensional array.
    n
        Number of entries to average over. n=2 means averaging over the currrent
        the previous entry.

    Returns
    -------
    An array view storing the moving average.
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1 :] / n


def get_random_state(seed: AnyRandom) -> np.random.RandomState:
    if isinstance(seed, np.random.RandomState):
        return seed
    return np.random.RandomState(seed)


# --------------------------------------------------------------------------------
# Deal with tool parameters
# --------------------------------------------------------------------------------


def update_params(
    old_params: Mapping[str, Any],
    new_params: Mapping[str, Any],
    check=False,
) -> dict[str, Any]:
    """\
    Update old_params with new_params.

    If check==False, this merely adds and overwrites the content of old_params.

    If check==True, this only allows updating of parameters that are already
    present in old_params.

    Parameters
    ----------
    old_params
    new_params
    check

    Returns
    -------
    updated_params
    """
    updated_params = dict(old_params)
    if new_params:  # allow for new_params to be None
        for key, val in new_params.items():
            if key not in old_params and check:
                raise ValueError(
                    "'"
                    + key
                    + "' is not a valid parameter key, "
                    + "consider one of \n"
                    + str(list(old_params.keys()))
                )
            if val is not None:
                updated_params[key] = val
    return updated_params


# --------------------------------------------------------------------------------
# Others
# --------------------------------------------------------------------------------


if TYPE_CHECKING:
    _SparseMatrix = Union[sparse.csr_matrix, sparse.csc_matrix]
    _MemoryArray = Union[NDArray, _SparseMatrix]
    _SupportedArray = Union[_MemoryArray, DaskArray]


@singledispatch
def elem_mul(x: _SupportedArray, y: _SupportedArray) -> _SupportedArray:
    raise NotImplementedError


@elem_mul.register(np.ndarray)
@elem_mul.register(sparse.spmatrix)
def _elem_mul_in_mem(x: _MemoryArray, y: _MemoryArray) -> _MemoryArray:
    if isinstance(x, sparse.spmatrix):
        # returns coo_matrix, so cast back to input type
        return type(x)(x.multiply(y))
    return x * y


@elem_mul.register(DaskArray)
def _elem_mul_dask(x: DaskArray, y: DaskArray) -> DaskArray:
    import dask.array as da

    return da.map_blocks(elem_mul, x, y)


if TYPE_CHECKING:
    Scaling_T = TypeVar("Scaling_T", DaskArray, np.ndarray)


def broadcast_axis(divisor: Scaling_T, axis: Literal[0, 1]) -> Scaling_T:
    divisor = np.ravel(divisor)
    if axis:
        return divisor[None, :]
    return divisor[:, None]


def check_op(op):
    if op not in {truediv, mul}:
        raise ValueError(f"{op} not one of truediv or mul")


@singledispatch
def axis_mul_or_truediv(
    X: np.ndarray,
    scaling_array: np.ndarray,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: np.ndarray | None = None,
) -> np.ndarray:
    check_op(op)
    scaling_array = broadcast_axis(scaling_array, axis)
    if op is mul:
        return np.multiply(X, scaling_array, out=out)
    if not allow_divide_by_zero:
        scaling_array = scaling_array.copy() + (scaling_array == 0)
    return np.true_divide(X, scaling_array, out=out)


@axis_mul_or_truediv.register(sparse.csr_matrix)
@axis_mul_or_truediv.register(sparse.csc_matrix)
def _(
    X: sparse.csr_matrix | sparse.csc_matrix,
    scaling_array,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: sparse.csr_matrix | sparse.csc_matrix | None = None,
) -> sparse.csr_matrix | sparse.csc_matrix:
    check_op(op)
    if out is not None:
        if X.data is not out.data:
            raise ValueError(
                "`out` argument provided but not equal to X.  This behavior is not supported for sparse matrix scaling."
            )
    if not allow_divide_by_zero and op is truediv:
        scaling_array = scaling_array.copy() + (scaling_array == 0)

    row_scale = axis == 0
    column_scale = axis == 1
    if row_scale:

        def new_data_op(x):
            return op(x.data, np.repeat(scaling_array, np.diff(x.indptr)))

    elif column_scale:

        def new_data_op(x):
            return op(x.data, scaling_array.take(x.indices, mode="clip"))

    if X.format == "csr":
        indices = X.indices
        indptr = X.indptr
        if out is not None:
            X.data = new_data_op(X)
            return X
        return sparse.csr_matrix(
            (new_data_op(X), indices.copy(), indptr.copy()), shape=X.shape
        )
    transposed = X.T
    return axis_mul_or_truediv(
        transposed,
        scaling_array,
        op=op,
        axis=1 - axis,
        out=transposed,
        allow_divide_by_zero=allow_divide_by_zero,
    ).T


def make_axis_chunks(
    X: DaskArray, axis: Literal[0, 1], pad=True
) -> tuple[tuple[int], tuple[int]]:
    if axis == 0:
        return (X.chunks[axis], (1,))
    return ((1,), X.chunks[axis])


@axis_mul_or_truediv.register(DaskArray)
def _(
    X: DaskArray,
    scaling_array: Scaling_T,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: None = None,
) -> DaskArray:
    check_op(op)
    if out is not None:
        raise TypeError(
            "`out` is not `None`. Do not do in-place modifications on dask arrays."
        )

    import dask.array as da

    scaling_array = broadcast_axis(scaling_array, axis)
    row_scale = axis == 0
    column_scale = axis == 1

    if isinstance(scaling_array, DaskArray):
        if (row_scale and not X.chunksize[0] == scaling_array.chunksize[0]) or (
            column_scale
            and (
                (
                    len(scaling_array.chunksize) == 1
                    and X.chunksize[1] != scaling_array.chunksize[0]
                )
                or (
                    len(scaling_array.chunksize) == 2
                    and X.chunksize[1] != scaling_array.chunksize[1]
                )
            )
        ):
            warnings.warn("Rechunking scaling_array in user operation", UserWarning)
            scaling_array = scaling_array.rechunk(make_axis_chunks(X, axis))
    else:
        scaling_array = da.from_array(
            scaling_array,
            chunks=make_axis_chunks(X, axis),
        )
    return da.map_blocks(
        axis_mul_or_truediv,
        X,
        scaling_array,
        axis,
        op,
        meta=X._meta,
        out=out,
        allow_divide_by_zero=allow_divide_by_zero,
    )


@overload
def axis_sum(
    X: sparse.spmatrix,
    *,
    axis: tuple[Literal[0, 1], ...] | Literal[0, 1] | None = None,
    dtype: DTypeLike | None = None,
) -> np.matrix: ...


@singledispatch
def axis_sum(
    X: np.ndarray,
    *,
    axis: tuple[Literal[0, 1], ...] | Literal[0, 1] | None = None,
    dtype: DTypeLike | None = None,
) -> np.ndarray:
    return np.sum(X, axis=axis, dtype=dtype)


@axis_sum.register(DaskArray)
def _(
    X: DaskArray,
    *,
    axis: tuple[Literal[0, 1], ...] | Literal[0, 1] | None = None,
    dtype: DTypeLike | None = None,
) -> DaskArray:
    import dask.array as da

    if dtype is None:
        dtype = getattr(np.zeros(1, dtype=X.dtype).sum(), "dtype", object)

    if isinstance(X._meta, np.ndarray) and not isinstance(X._meta, np.matrix):
        return X.sum(axis=axis, dtype=dtype)

    def sum_drop_keepdims(*args, **kwargs):
        kwargs.pop("computing_meta", None)
        # masked operations on sparse produce which numpy matrices gives the same API issues handled here
        if isinstance(X._meta, (sparse.spmatrix, np.matrix)) or isinstance(
            args[0], (sparse.spmatrix, np.matrix)
        ):
            kwargs.pop("keepdims", None)
            axis = kwargs["axis"]
            if isinstance(axis, tuple):
                if len(axis) != 1:
                    raise ValueError(
                        f"`axis_sum` can only sum over one axis when `axis` arg is provided but got {axis} instead"
                    )
                kwargs["axis"] = axis[0]
        # returns a np.matrix normally, which is undesireable
        return np.array(np.sum(*args, dtype=dtype, **kwargs))

    def aggregate_sum(*args, **kwargs):
        return np.sum(args[0], dtype=dtype, **kwargs)

    return da.reduction(
        X,
        sum_drop_keepdims,
        aggregate_sum,
        axis=axis,
        dtype=dtype,
        meta=np.array([], dtype=dtype),
    )


@singledispatch
def check_nonnegative_integers(X: _SupportedArray) -> bool | DaskArray:
    """Checks values of X to ensure it is count data"""
    raise NotImplementedError


@check_nonnegative_integers.register(np.ndarray)
@check_nonnegative_integers.register(sparse.spmatrix)
def _check_nonnegative_integers_in_mem(X: _MemoryArray) -> bool:
    from numbers import Integral

    data = X if isinstance(X, np.ndarray) else X.data
    # Check no negatives
    if np.signbit(data).any():
        return False
    # Check all are integers
    elif issubclass(data.dtype.type, Integral):
        return True
    return not np.any((data % 1) != 0)


@check_nonnegative_integers.register(DaskArray)
def _check_nonnegative_integers_dask(X: DaskArray) -> DaskArray:
    return X.map_blocks(check_nonnegative_integers, dtype=bool, drop_axis=(0, 1))


def select_groups(
    adata: AnnData,
    groups_order_subset: list[str] | Literal["all"] = "all",
    key: str = "groups",
) -> tuple[list[str], NDArray[np.bool_]]:
    """Get subset of groups in adata.obs[key]."""
    groups_order = adata.obs[key].cat.categories
    if key + "_masks" in adata.uns:
        groups_masks_obs = adata.uns[key + "_masks"]
    else:
        groups_masks_obs = np.zeros(
            (len(adata.obs[key].cat.categories), adata.obs[key].values.size), dtype=bool
        )
        for iname, name in enumerate(adata.obs[key].cat.categories):
            # if the name is not found, fallback to index retrieval
            if adata.obs[key].cat.categories[iname] in adata.obs[key].values:
                mask_obs = adata.obs[key].cat.categories[iname] == adata.obs[key].values
            else:
                mask_obs = str(iname) == adata.obs[key].values
            groups_masks_obs[iname] = mask_obs
    groups_ids = list(range(len(groups_order)))
    if groups_order_subset != "all":
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(
                np.where(adata.obs[key].cat.categories.values == name)[0][0]
            )
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(
                np.in1d(
                    np.arange(len(adata.obs[key].cat.categories)).astype(str),
                    np.array(groups_order_subset),
                )
            )[0]
        if len(groups_ids) == 0:
            logg.debug(
                f"{np.array(groups_order_subset)} invalid! specify valid "
                f"groups_order (or indices) from {adata.obs[key].cat.categories}",
            )
            from sys import exit

            exit(0)
        groups_masks_obs = groups_masks_obs[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks_obs


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):  # noqa: PLR0917
    """Get full tracebacks when warning is raised by setting

    warnings.showwarning = warn_with_traceback

    See also
    --------
    https://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    import traceback

    traceback.print_stack()
    log = (  # noqa: F841  # TODO Does this need fixing?
        file if hasattr(file, "write") else sys.stderr
    )
    settings.write(warnings.formatwarning(message, category, filename, lineno, line))


def warn_once(msg: str, category: type[Warning], stacklevel: int = 1):
    warnings.warn(msg, category, stacklevel=stacklevel)
    # You'd think `'once'` works, but it doesn't at the repl and in notebooks
    warnings.filterwarnings("ignore", category=category, message=re.escape(msg))


def subsample(
    X: np.ndarray,
    subsample: int = 1,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """\
    Subsample a fraction of 1/subsample samples from the rows of X.

    Parameters
    ----------
    X
        Data array.
    subsample
        1/subsample is the fraction of data sampled, n = X.shape[0]/subsample.
    seed
        Seed for sampling.

    Returns
    -------
    Xsampled
        Subsampled X.
    rows
        Indices of rows that are stored in Xsampled.
    """
    if subsample == 1 and seed == 0:
        return X, np.arange(X.shape[0], dtype=int)
    if seed == 0:
        # this sequence is defined simply by skipping rows
        # is faster than sampling
        rows = np.arange(0, X.shape[0], subsample, dtype=int)
        n = rows.size
        Xsampled = np.array(X[rows])
    else:
        if seed < 0:
            raise ValueError(f"Invalid seed value < 0: {seed}")
        n = int(X.shape[0] / subsample)
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X, n=n)
    logg.debug(f"... subsampled to {n} of {X.shape[0]} data points")
    return Xsampled, rows


def subsample_n(
    X: np.ndarray, n: int = 0, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Subsample n samples from rows of array.

    Parameters
    ----------
    X
        Data array.
    n
        Sample size.
    seed
        Seed for sampling.

    Returns
    -------
    Xsampled
        Subsampled X.
    rows
        Indices of rows that are stored in Xsampled.
    """
    if n < 0:
        raise ValueError("n must be greater 0")
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0], size=n, replace=False)
    Xsampled = X[rows]
    return Xsampled, rows


def check_presence_download(filename: Path, backup_url):
    """Check if file is present otherwise download."""
    if not filename.is_file():
        from ..readwrite import _download

        _download(backup_url, filename)


def lazy_import(full_name):
    """Imports a module in a way that it’s only executed on member access"""
    try:
        return sys.modules[full_name]
    except KeyError:
        spec = importlib.util.find_spec(full_name)
        module = importlib.util.module_from_spec(spec)
        loader = importlib.util.LazyLoader(spec.loader)
        # Make module with proper locking and get it inserted into sys.modules.
        loader.exec_module(module)
        return module


# --------------------------------------------------------------------------------
# Neighbors
# --------------------------------------------------------------------------------


def _fallback_to_uns(dct, conns, dists, conns_key, dists_key):
    if conns is None and conns_key in dct:
        conns = dct[conns_key]
    if dists is None and dists_key in dct:
        dists = dct[dists_key]

    return conns, dists


class NeighborsView:
    """Convenience class for accessing neighbors graph representations.

    Allows to access neighbors distances, connectivities and settings
    dictionary in a uniform manner.

    Parameters
    ----------

    adata
        AnnData object.
    key
        This defines where to look for neighbors dictionary,
        connectivities, distances.

        neigh = NeighborsView(adata, key)
        neigh['distances']
        neigh['connectivities']
        neigh['params']
        'connectivities' in neigh
        'params' in neigh

        is the same as

        adata.obsp[adata.uns[key]['distances_key']]
        adata.obsp[adata.uns[key]['connectivities_key']]
        adata.uns[key]['params']
        adata.uns[key]['connectivities_key'] in adata.obsp
        'params' in adata.uns[key]
    """

    def __init__(self, adata, key=None):
        self._connectivities = None
        self._distances = None

        if key is None or key == "neighbors":
            if "neighbors" not in adata.uns:
                raise KeyError('No "neighbors" in .uns')
            self._neighbors_dict = adata.uns["neighbors"]
            self._conns_key = "connectivities"
            self._dists_key = "distances"
        else:
            if key not in adata.uns:
                raise KeyError(f'No "{key}" in .uns')
            self._neighbors_dict = adata.uns[key]
            self._conns_key = self._neighbors_dict["connectivities_key"]
            self._dists_key = self._neighbors_dict["distances_key"]

        if self._conns_key in adata.obsp:
            self._connectivities = adata.obsp[self._conns_key]
        if self._dists_key in adata.obsp:
            self._distances = adata.obsp[self._dists_key]

        # fallback to uns
        self._connectivities, self._distances = _fallback_to_uns(
            self._neighbors_dict,
            self._connectivities,
            self._distances,
            self._conns_key,
            self._dists_key,
        )

    def __getitem__(self, key):
        if key == "distances":
            if "distances" not in self:
                raise KeyError(f'No "{self._dists_key}" in .obsp')
            return self._distances
        elif key == "connectivities":
            if "connectivities" not in self:
                raise KeyError(f'No "{self._conns_key}" in .obsp')
            return self._connectivities
        else:
            return self._neighbors_dict[key]

    def __contains__(self, key):
        if key == "distances":
            return self._distances is not None
        elif key == "connectivities":
            return self._connectivities is not None
        else:
            return key in self._neighbors_dict


def _choose_graph(adata, obsp, neighbors_key):
    """Choose connectivities from neighbbors or another obsp column"""
    if obsp is not None and neighbors_key is not None:
        raise ValueError(
            "You can't specify both obsp, neighbors_key. " "Please select only one."
        )

    if obsp is not None:
        return adata.obsp[obsp]
    else:
        neighbors = NeighborsView(adata, neighbors_key)
        if "connectivities" not in neighbors:
            raise ValueError(
                "You need to run `pp.neighbors` first "
                "to compute a neighborhood graph."
            )
        return neighbors["connectivities"]


def _resolve_axis(
    axis: Literal["obs", 0, "var", 1],
) -> tuple[Literal[0], Literal["obs"]] | tuple[Literal[1], Literal["var"]]:
    if axis in {0, "obs"}:
        return (0, "obs")
    if axis in {1, "var"}:
        return (1, "var")
    raise ValueError(f"`axis` must be either 0, 1, 'obs', or 'var', was {axis!r}")


def is_backed_type(X: object) -> bool:
    return isinstance(X, (SparseDataset, h5py.File, h5py.Dataset))


def raise_not_implemented_error_if_backed_type(X: object, method_name: str) -> None:
    if is_backed_type(X):
        raise NotImplementedError(
            f"{method_name} is not implemented for matrices of type {type(X)}"
        )
