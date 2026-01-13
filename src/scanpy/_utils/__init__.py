"""Utility functions and classes.

This file largely consists of the old _utils.py file. Over time, these functions
should be moved of this file.
"""

from __future__ import annotations

import importlib.util
import inspect
import re
from contextlib import suppress
from enum import Enum
from functools import partial, reduce, singledispatch, wraps
from operator import mul, or_, truediv
from textwrap import indent
from types import MethodType, ModuleType, UnionType
from typing import (
    TYPE_CHECKING,
    Literal,
    NamedTuple,
    TypeAliasType,
    Union,
    get_args,
    get_origin,
    overload,
)
from weakref import WeakSet

import h5py
import numpy as np
import pandas as pd
from anndata._core.sparse_dataset import BaseCompressedSparseDataset
from packaging.version import Version

from .. import logging as logg
from .._compat import CSBase, DaskArray, _CSArray, pkg_version, warn
from .._settings import settings

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, KeysView, Mapping
    from pathlib import Path
    from typing import Any

    from anndata import AnnData
    from igraph import Graph
    from numpy.typing import ArrayLike, NDArray
    from pandas._typing import Dtype as PdDtype

    from .._compat import CSRBase
    from ..neighbors import NeighborsParams, RPForestDict

    type _MemoryArray = NDArray | CSBase
    type _SupportedArray = _MemoryArray | DaskArray


__all__ = [
    "AssoResult",
    "Empty",
    "NeighborsView",
    "_choose_graph",
    "_doc_params",
    "_empty",
    "_resolve_axis",
    "annotate_doc_types",
    "axis_mul_or_truediv",
    "axis_nnz",
    "check_array_function_arguments",
    "check_nonnegative_integers",
    "check_presence_download",
    "check_use_raw",
    "compute_association_matrix_of_groups",
    "descend_classes_and_funcs",
    "ensure_igraph",
    "get_literal_vals",
    "indent",
    "is_backed_type",
    "is_backed_type",
    "raise_not_implemented_error_if_backed_type",
    "renamed_arg",
    "sanitize_anndata",
    "select_groups",
    "update_params",
    "with_cat_dtype",
]


LegacyUnionType: type = type(Union[int, str])  # noqa: UP007


class Empty(Enum):
    token = 0

    def __repr__(self) -> str:
        return "_empty"


_empty = Empty.token


def ensure_igraph() -> None:
    if importlib.util.find_spec("igraph"):
        return
    msg = (
        "Please install the igraph package: "
        "`conda install -c conda-forge python-igraph` or "
        "`pip install igraph`."
    )
    raise ImportError(msg)


def _getdoc(c_or_f: Callable | type) -> str | None:
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
            __tracebackhide__ = True
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
                warn(msg, FutureWarning)
                if pos_0:
                    args = (kwargs.pop(old_name), *args)
                else:
                    kwargs[new_name] = kwargs.pop(old_name)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def import_name(full_name: str) -> Any:
    from importlib import import_module

    parts = full_name.split(".")
    obj = import_module(parts[0])
    for _i, name in enumerate(parts[1:]):
        i = _i
        try:
            obj = import_module(f"{obj.__name__}.{name}")
        except ModuleNotFoundError:
            break
    else:
        i = len(parts)
    for name in parts[i + 1 :]:
        try:
            obj = getattr(obj, name)
        except AttributeError as e:
            msg = f"{parts[:i]}, {parts[i + 1 :]}, {obj} {name}"
            raise RuntimeError(msg) from e
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
        try:
            encountered.add(obj)
        except TypeError:
            continue  # TypeAliasTypes etc. are not weakref-able
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
        with suppress(AttributeError):
            c_or_f.getdoc = partial(_getdoc, c_or_f)


_leading_whitespace_re = re.compile("(^[ ]*)(?:[^ \n])", re.MULTILINE)


def _doc_params[T: Callable | type](**replacements: str) -> Callable[[T], T]:
    def dec(obj: T) -> T:
        assert obj.__doc__
        assert "\t" not in obj.__doc__

        # The first line of the docstring is unindented,
        # so find indent size starting after it.
        start_line_2 = obj.__doc__.find("\n") + 1
        assert start_line_2 > 0, f"{obj.__name__} has single-line docstring."
        n_spaces = min(
            len(m.group(1))
            for m in _leading_whitespace_re.finditer(obj.__doc__[start_line_2:])
        )

        # The placeholder is already indented, so only indent subsequent lines
        indented_replacements = {
            k: indent(v, " " * n_spaces)[n_spaces:] for k, v in replacements.items()
        }
        obj.__doc__ = obj.__doc__.format_map(indented_replacements)
        return obj

    return dec


def check_array_function_arguments(**kwargs):
    """Check for invalid arguments when an array is passed.

    Helper for functions that work on either AnnData objects or array-likes.
    """
    # TODO: Figure out a better solution for documenting dispatched functions
    invalid_args = [k for k, v in kwargs.items() if v is not None]
    if len(invalid_args) > 0:
        msg = f"Arguments {invalid_args} are only valid if an AnnData object is passed."
        raise TypeError(msg)


def check_use_raw(
    adata: AnnData,
    use_raw: None | bool,  # noqa: FBT001
    *,
    layer: str | None = None,
) -> bool:
    """Normalize checking `use_raw`.

    My intentention here is to also provide a single place to throw a deprecation warning from in future.
    """
    if use_raw is not None:
        return use_raw
    if layer is not None:
        return False
    return adata.raw is not None


# --------------------------------------------------------------------------------
# Graph stuff
# --------------------------------------------------------------------------------


def get_igraph_from_adjacency(adjacency: CSBase, *, directed: bool = False) -> Graph:
    """Get igraph graph from adjacency matrix."""
    import igraph as ig

    sources, targets = adjacency.nonzero()
    weights = dematrix(adjacency[sources, targets]).ravel() if len(sources) else []
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets, strict=True)))
    with suppress(KeyError):
        g.es["weight"] = weights
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f"The constructed graph has only {g.vcount()} nodes. "
            "Your adjacency matrix contained redundant nodes."
        )
    return g


# --------------------------------------------------------------------------------
# Group stuff
# --------------------------------------------------------------------------------


class AssoResult(NamedTuple):
    asso_names: list[str]
    asso_matrix: NDArray[np.floating]


def compute_association_matrix_of_groups(
    adata: AnnData,
    prediction: str,
    reference: str,
    *,
    normalization: Literal["prediction", "reference"] = "prediction",
    threshold: float = 0.01,
    max_n_names: int | None = 2,
) -> AssoResult:
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
        msg = '`normalization` needs to be either "prediction" or "reference".'
        raise ValueError(msg)
    sanitize_anndata(adata)
    cats = adata.obs[reference].cat.categories
    for cat in cats:
        if cat in settings.categories_to_ignore:
            logg.info(
                f"Ignoring category {cat!r} as it’s in `settings.categories_to_ignore`."
            )
    asso_names: list[str] = []
    asso_matrix: list[list[float]] = []
    for ipred_group, pred_group in enumerate(adata.obs[prediction].cat.categories):
        if "?" in pred_group:
            pred_group = str(ipred_group)  # noqa: PLW2901
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
    return AssoResult(asso_names=asso_names, asso_matrix=np.array(asso_matrix))


def get_associated_colors_of_groups(
    reference_colors: Mapping[int, str], asso_matrix: NDArray[np.floating]
) -> list[dict[str, float]]:
    return [
        {
            reference_colors[i_ref]: asso_matrix[i_pred, i_ref]
            for i_ref in range(asso_matrix.shape[1])
        }
        for i_pred in range(asso_matrix.shape[0])
    ]


def identify_groups(ref_labels, pred_labels, *, return_overlaps: bool = False):
    """Identify which predicted label explains which reference label.

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
    ref_dict = dict(zip(ref_unique, ref_counts, strict=True))
    pred_unique, pred_counts = np.unique(pred_labels, return_counts=True)
    pred_dict = dict(zip(pred_unique, pred_counts, strict=True))
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
        msg = "Received a view of an AnnData. Making a copy."
        warn(msg, UserWarning)
        adata._init_as_actual(adata.copy())


def moving_average(a: np.ndarray, n: int):
    """Moving average over one-dimensional array.

    Parameters
    ----------
    a
        One-dimensional array.
    n
        Number of entries to average over. n=2 means averaging over the current
        the previous entry.

    Returns
    -------
    An array view storing the moving average.

    """  # noqa: D401
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1 :] / n


@singledispatch
def with_cat_dtype[X: pd.Series | pd.CategoricalIndex | pd.Categorical](
    x: X, dtype: PdDtype
) -> X:
    raise NotImplementedError


@with_cat_dtype.register(pd.Series)
def _(x: pd.Series, dtype: PdDtype) -> pd.Series:
    return x.cat.set_categories(x.cat.categories.astype(dtype))


@with_cat_dtype.register(pd.Categorical | pd.CategoricalIndex)
def _[X: pd.Categorical | pd.CategoricalIndex](x: X, dtype: PdDtype) -> X:
    return x.set_categories(x.categories.astype(dtype))


# --------------------------------------------------------------------------------
# Deal with tool parameters
# --------------------------------------------------------------------------------


def update_params(
    old_params: Mapping[str, Any],
    new_params: Mapping[str, Any],
    *,
    check: bool = False,
) -> dict[str, Any]:
    """Update `old_params` with `new_params`.

    If check==False, this merely adds and overwrites the content of `old_params`.

    If check==True, this only allows updating of parameters that are already
    present in `old_params`.

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
                msg = (
                    f"{key!r} is not a valid parameter key, "
                    f"consider one of \n{list(old_params.keys())}"
                )
                raise ValueError(msg)
            if val is not None:
                updated_params[key] = val
    return updated_params


# `get_args` returns `tuple[Any]` so I don’t think it’s possible to get the correct type here
def get_literal_vals(typ: UnionType | TypeAliasType | Any) -> KeysView[Any]:
    """Get all literal values from a Literal or Union of … of Literal type."""
    if isinstance(typ, UnionType | LegacyUnionType):
        return reduce(
            or_, (dict.fromkeys(get_literal_vals(t)) for t in get_args(typ))
        ).keys()
    if isinstance(typ, TypeAliasType):
        return get_literal_vals(typ.__value__)
    if get_origin(typ) is Literal:
        return dict.fromkeys(get_args(typ)).keys()
    msg = f"{typ!r} ({type(typ).__name__}) is not a valid Literal"
    raise TypeError(msg)


# --------------------------------------------------------------------------------
# Others
# --------------------------------------------------------------------------------


def _broadcast_axis[T: (DaskArray, np.ndarray)](divisor: T, axis: Literal[0, 1]) -> T:
    divisor = np.ravel(divisor)
    if axis:
        return divisor[None, :]
    return divisor[:, None]


def _check_op(op) -> None:
    if op not in {truediv, mul}:
        msg = f"{op} not one of truediv or mul"
        raise ValueError(msg)


@singledispatch
def axis_mul_or_truediv(
    x: ArrayLike,
    /,
    scaling_array: np.ndarray,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: ArrayLike | None = None,
) -> np.ndarray:
    _check_op(op)
    scaling_array = _broadcast_axis(scaling_array, axis)
    if op is mul:
        return np.multiply(x, scaling_array, out=out)
    if not allow_divide_by_zero:
        scaling_array = scaling_array.copy() + (scaling_array == 0)
    return np.true_divide(x, scaling_array, out=out)


@axis_mul_or_truediv.register(CSBase)
def _(
    x: CSBase,
    /,
    scaling_array: np.ndarray,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: CSBase | None = None,
) -> CSBase:
    _check_op(op)
    if out is not None and x.data is not out.data:
        msg = "`out` argument provided but not equal to X.  This behavior is not supported for sparse matrix scaling."
        raise ValueError(msg)
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

    if x.format == "csr":
        indices = x.indices
        indptr = x.indptr
        if out is not None:
            x.data = new_data_op(x)
            return x
        return type(x)((new_data_op(x), indices.copy(), indptr.copy()), shape=x.shape)
    transposed = x.T
    return axis_mul_or_truediv(
        transposed,
        scaling_array,
        op=op,
        axis=1 - axis,
        out=transposed,
        allow_divide_by_zero=allow_divide_by_zero,
    ).T


def _make_axis_chunks(
    x: DaskArray, axis: Literal[0, 1]
) -> tuple[tuple[int], tuple[int]]:
    if axis == 0:
        return (x.chunks[axis], (1,))
    return ((1,), x.chunks[axis])


@axis_mul_or_truediv.register(DaskArray)
def _[T: (DaskArray, np.ndarray)](
    x: DaskArray,
    /,
    scaling_array: T,
    axis: Literal[0, 1],
    op: Callable[[Any, Any], Any],
    *,
    allow_divide_by_zero: bool = True,
    out: None = None,
) -> DaskArray:
    _check_op(op)
    if out is not None:
        msg = "`out` is not `None`. Do not do in-place modifications on dask arrays."
        raise TypeError(msg)

    import dask.array as da

    scaling_array = _broadcast_axis(scaling_array, axis)
    row_scale = axis == 0
    column_scale = axis == 1

    if isinstance(scaling_array, DaskArray):
        if (row_scale and x.chunksize[0] != scaling_array.chunksize[0]) or (
            column_scale
            and (
                (
                    len(scaling_array.chunksize) == 1
                    and x.chunksize[1] != scaling_array.chunksize[0]
                )
                or (
                    len(scaling_array.chunksize) == 2
                    and x.chunksize[1] != scaling_array.chunksize[1]
                )
            )
        ):
            msg = "Rechunking scaling_array in user operation"
            warn(msg, UserWarning)
            scaling_array = scaling_array.rechunk(_make_axis_chunks(x, axis))
    else:
        scaling_array = da.from_array(
            scaling_array,
            chunks=_make_axis_chunks(x, axis),
        )
    return da.map_blocks(
        axis_mul_or_truediv,
        x,
        scaling_array,
        axis,
        op,
        meta=x._meta,
        out=out,
        allow_divide_by_zero=allow_divide_by_zero,
    )


@singledispatch
def axis_nnz(x: ArrayLike, /, axis: Literal[0, 1]) -> np.ndarray:
    return np.count_nonzero(x, axis=axis)


if pkg_version("scipy") >= Version("1.15"):
    # newer scipy versions support the `axis` argument for count_nonzero
    @axis_nnz.register(CSBase)
    def _(x: CSBase, /, axis: Literal[0, 1]) -> np.ndarray:
        return x.count_nonzero(axis=axis)

else:
    # older scipy versions don’t have any way to get the nnz of a sparse array
    @axis_nnz.register(CSBase)
    def _(x: CSBase, /, axis: Literal[0, 1]) -> np.ndarray:
        if isinstance(x, _CSArray):
            from scipy.sparse import csc_array, csr_array  # noqa: TID251

            x = (csr_array if x.format == "csr" else csc_array)(x)
        return x.getnnz(axis=axis)


@axis_nnz.register(DaskArray)
def _(x: DaskArray, /, axis: Literal[0, 1]) -> DaskArray:
    return x.map_blocks(
        partial(axis_nnz, axis=axis),
        dtype=np.int64,
        meta=np.array([], dtype=np.int64),
        drop_axis=axis,
    )


@singledispatch
def check_nonnegative_integers(x: _SupportedArray, /) -> bool | DaskArray:
    """Check values of X to ensure it is count data."""
    raise NotImplementedError


@check_nonnegative_integers.register(np.ndarray)
@check_nonnegative_integers.register(CSBase)
def _check_nonnegative_integers_in_mem(x: _MemoryArray, /) -> bool:
    from numbers import Integral

    data = x if isinstance(x, np.ndarray) else x.data
    # Check no negatives
    if np.signbit(data).any():
        return False
    # Check all are integers
    elif issubclass(data.dtype.type, Integral):
        return True
    return not np.any((data % 1) != 0)


@check_nonnegative_integers.register(DaskArray)
def _check_nonnegative_integers_dask(x: DaskArray, /) -> DaskArray:
    return x.map_blocks(check_nonnegative_integers, dtype=bool, drop_axis=(0, 1))


def dematrix[SA: _SupportedArray](x: SA | np.matrix) -> SA:
    if isinstance(x, np.matrix):
        return x.A
    if isinstance(x, DaskArray) and isinstance(x._meta, np.matrix):
        return x.map_blocks(np.asarray, meta=np.array([], dtype=x.dtype))
    return x


def raise_if_dask_feature_axis_chunked(x: Any):
    if isinstance(x, DaskArray) and x.chunksize[1] != x.shape[1]:
        msg = (
            "Only dask arrays with chunking along the first axis are supported. "
            f"Got chunksize {x.chunksize} with shape {x.shape}. "
        )
        raise ValueError(msg)


def select_groups(
    adata: AnnData,
    groups_order_subset: Iterable[str] | Literal["all"] = "all",
    key: str = "groups",
) -> tuple[list[str], NDArray[np.bool_]]:
    """Get subset of groups in adata.obs[key]."""
    groups_order = adata.obs[key].cat.categories
    if f"{key}_masks" in adata.uns:
        groups_masks_obs = adata.uns[f"{key}_masks"]
    else:
        groups_masks_obs = np.zeros(
            (len(adata.obs[key].cat.categories), adata.obs[key].values.size), dtype=bool
        )
        for iname, name in enumerate(adata.obs[key].cat.categories):
            # if the name is not found, fallback to index retrieval
            if name in adata.obs[key].values:
                mask_obs = name == adata.obs[key].values
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
                np.isin(
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


def check_presence_download(filename: Path, backup_url):
    """Check if file is present otherwise download."""
    if not filename.is_file():
        from ..readwrite import _download

        _download(backup_url, filename)


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

    Examples
    --------
    >>> neigh = NeighborsView(adata, key)
    >>> neigh["distances"]
    >>> neigh["connectivities"]
    >>> neigh["params"]
    >>> "connectivities" in neigh
    >>> "params" in neigh

    is the same as

    >>> adata.obsp[adata.uns[key]["distances_key"]]
    >>> adata.obsp[adata.uns[key]["connectivities_key"]]
    >>> adata.uns[key]["params"]
    >>> adata.uns[key]["connectivities_key"] in adata.obsp
    >>> "params" in adata.uns[key]

    """

    def __init__(self, adata: AnnData, key=None):
        self._connectivities = None
        self._distances = None

        if key is None or key == "neighbors":
            if "neighbors" not in adata.uns:
                msg = 'No "neighbors" in .uns'
                raise KeyError(msg)
            self._neighbors_dict = adata.uns["neighbors"]
            self._conns_key = "connectivities"
            self._dists_key = "distances"
        else:
            if key not in adata.uns:
                msg = f"No {key!r} in .uns"
                raise KeyError(msg)
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

    @overload
    def __getitem__(self, key: Literal["distances", "connectivities"]) -> CSRBase: ...
    @overload
    def __getitem__(self, key: Literal["params"]) -> NeighborsParams: ...
    @overload
    def __getitem__(self, key: Literal["rp_forest"]) -> RPForestDict: ...
    @overload
    def __getitem__(self, key: Literal["connectivities_key"]) -> str: ...

    def __getitem__(self, key: str):
        if key == "distances":
            if "distances" not in self:
                msg = f"No {self._dists_key!r} in .obsp"
                raise KeyError(msg)
            return self._distances
        elif key == "connectivities":
            if "connectivities" not in self:
                msg = f"No {self._conns_key!r} in .obsp"
                raise KeyError(msg)
            return self._connectivities
        elif key == "connectivities_key":
            return self._conns_key
        else:
            return self._neighbors_dict[key]

    def __contains__(self, key: str) -> bool:
        if key == "distances":
            return self._distances is not None
        elif key == "connectivities":
            return self._connectivities is not None
        else:
            return key in self._neighbors_dict


def _choose_graph(
    adata: AnnData, obsp: str | None, neighbors_key: str | None
) -> CSBase:
    """Choose connectivities from neighbors or another obsp entry."""
    if obsp is not None and neighbors_key is not None:
        msg = "You can't specify both obsp, neighbors_key. Please select only one."
        raise ValueError(msg)

    if obsp is not None:
        return adata.obsp[obsp]
    else:
        neighbors = NeighborsView(adata, neighbors_key)
        if "connectivities" not in neighbors:
            msg = (
                "You need to run `pp.neighbors` first to compute a neighborhood graph."
            )
            raise ValueError(msg)
        return neighbors["connectivities"]


def _resolve_axis(
    axis: Literal["obs", 0, "var", 1],
) -> tuple[Literal[0], Literal["obs"]] | tuple[Literal[1], Literal["var"]]:
    if axis in {0, "obs"}:
        return (0, "obs")
    if axis in {1, "var"}:
        return (1, "var")
    msg = f"`axis` must be either 0, 1, 'obs', or 'var', was {axis!r}"
    raise ValueError(msg)


def is_backed_type(x: object, /) -> bool:
    return isinstance(x, BaseCompressedSparseDataset | h5py.File | h5py.Dataset)


def raise_not_implemented_error_if_backed_type(x: object, method_name: str, /) -> None:
    if is_backed_type(x):
        msg = f"{method_name} is not implemented for matrices of type {type(x)}"
        raise NotImplementedError(msg)
