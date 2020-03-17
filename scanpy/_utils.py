"""Utility functions and classes
"""
import sys
import inspect
import warnings
import importlib.util
from enum import Enum
from pathlib import Path
from weakref import WeakSet
from collections import namedtuple
from functools import partial, wraps
from types import ModuleType, MethodType
from typing import Union, Callable, Optional, Mapping, Any, Dict, Tuple

import numpy as np
from numpy import random
from anndata import AnnData
from textwrap import dedent
from packaging import version

from ._settings import settings
from ._compat import Literal
from . import logging as logg


class Empty(Enum):
    token = 0


_empty = Empty.token

# e.g. https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
AnyRandom = Union[None, int, random.RandomState]  # maybe in the future random.Generator

EPS = 1e-15


def pkg_version(package):
    try:
        from importlib.metadata import version as v
    except ImportError:  # < Python 3.8: Use backport module
        from importlib_metadata import version as v
    return version.parse(v(package))


def check_versions():
    anndata_version = pkg_version("anndata")
    umap_version = pkg_version("umap-learn")

    if anndata_version < version.parse('0.6.10'):
        from . import __version__
        raise ImportError(
            f'Scanpy {__version__} needs anndata version >=0.6.10, '
            f'not {anndata_version}.\nRun `pip install anndata -U --no-deps`.'
        )

    if umap_version < version.parse('0.3.0'):
        from . import __version__
        # make this a warning, not an error
        # it might be useful for people to still be able to run it
        logg.warning(
            f'Scanpy {__version__} needs umap '
            f'version >=0.3.0, not {umap_version}.'
        )


def getdoc(c_or_f: Union[Callable, type]) -> Optional[str]:
    if getattr(c_or_f, '__doc__', None) is None:
        return None
    doc = inspect.getdoc(c_or_f)
    if isinstance(c_or_f, type) and hasattr(c_or_f, '__init__'):
        sig = inspect.signature(c_or_f.__init__)
    else:
        sig = inspect.signature(c_or_f)

    def type_doc(name: str):
        param: inspect.Parameter = sig.parameters[name]
        cls = getattr(param.annotation, '__qualname__', repr(param.annotation))
        if param.default is not param.empty:
            return f'{cls}, optional (default: {param.default!r})'
        else:
            return cls

    return '\n'.join(
        f'{line} : {type_doc(line)}' if line.strip() in sig.parameters else line
        for line in doc.split('\n')
    )


def deprecated_arg_names(arg_mapping: Mapping[str, str]):
    """
    Decorator which marks a functions keyword arguments as deprecated. It will
    result in a warning being emitted when the deprecated keyword argument is
    used, and the function being called with the new argument.

    Parameters
    ----------
    arg_mapping
        Mapping from deprecated argument name to current argument name.
    """
    def decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            warnings.simplefilter(
                'always', DeprecationWarning)  # turn off filter
            for old, new in arg_mapping.items():
                if old in kwargs:
                    warnings.warn(
                        f"Keyword argument '{old}' has been "
                        f"deprecated in favour of '{new}'. "
                        f"'{old}' will be removed in a future version.",
                        category=DeprecationWarning,
                        stacklevel=2,
                    )
                    val = kwargs.pop(old)
                    kwargs[new] = val
            # reset filter
            warnings.simplefilter('default', DeprecationWarning)
            return func(*args, **kwargs)
        return func_wrapper
    return decorator


def _one_of_ours(obj, root: str):
    return (
        hasattr(obj, "__name__")
        and not obj.__name__.split(".")[-1].startswith("_")
        and getattr(obj, '__module__', getattr(obj, '__qualname__', obj.__name__)).startswith(root)
    )


def descend_classes_and_funcs(mod: ModuleType, root: str, encountered=None):
    if encountered is None:
        encountered = WeakSet()
    for obj in vars(mod).values():
        if not _one_of_ours(obj, root):
            continue
        if callable(obj) and not isinstance(obj, MethodType):
            yield obj
            if isinstance(obj, type):
                for m in vars(obj).values():
                    if callable(m) and _one_of_ours(m, root):
                        yield m
        elif isinstance(obj, ModuleType) and obj not in encountered:
            encountered.add(obj)
            yield from descend_classes_and_funcs(obj, root, encountered)


def annotate_doc_types(mod: ModuleType, root: str):
    for c_or_f in descend_classes_and_funcs(mod, root):
        c_or_f.getdoc = partial(getdoc, c_or_f)


def _doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """
    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj
    return dec


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
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g


def get_sparse_from_igraph(graph, weight_attr=None):
    from scipy.sparse import csr_matrix
    edges = graph.get_edgelist()
    if weight_attr is None:
        weights = [1] * len(edges)
    else:
        weights = graph.es[weight_attr]
    if not graph.is_directed():
        edges.extend([(v, u) for u, v in edges])
        weights.extend(weights)
    shape = graph.vcount()
    shape = (shape, shape)
    if len(edges) > 0:
        return csr_matrix((weights, zip(*edges)), shape=shape)
    else:
        return csr_matrix(shape)


# --------------------------------------------------------------------------------
# Group stuff
# --------------------------------------------------------------------------------


def compute_association_matrix_of_groups(
    adata: AnnData,
    prediction: str,
    reference: str,
    normalization: Literal['prediction', 'reference'] = 'prediction',
    threshold: float = 0.01,
    max_n_names: Optional[int] = 2,
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
    if normalization not in {'prediction', 'reference'}:
        raise ValueError('`normalization` needs to be either "prediction" or "reference".')
    sanitize_anndata(adata)
    cats = adata.obs[reference].cat.categories
    for cat in cats:
        if cat in settings.categories_to_ignore:
            logg.info(
                f'Ignoring category {cat!r} '
                'as it’s in `settings.categories_to_ignore`.'
            )
    asso_names = []
    asso_matrix = []
    for ipred_group, pred_group in enumerate(
            adata.obs[prediction].cat.categories):
        if '?' in pred_group: pred_group = str(ipred_group)
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
            if normalization == 'prediction':
                # compute which fraction of the predicted group is contained in
                # the ref group
                ratio_contained = (np.sum(mask_pred_int) -
                    np.sum(mask_ref_or_pred - mask_ref)) / np.sum(mask_pred_int)
            else:
                # compute which fraction of the reference group is contained in
                # the predicted group
                ratio_contained = (np.sum(mask_ref) -
                    np.sum(mask_ref_or_pred - mask_pred_int)) / np.sum(mask_ref)
            asso_matrix[-1] += [ratio_contained]
        name_list_pred = [cats[i] if cats[i] not in settings.categories_to_ignore else ''
                          for i in np.argsort(asso_matrix[-1])[::-1]
                          if asso_matrix[-1][i] > threshold]
        asso_names += ['\n'.join(name_list_pred[:max_n_names])]
    Result = namedtuple('compute_association_matrix_of_groups',
                        ['asso_names', 'asso_matrix'])
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
        sub_pred_unique, sub_pred_counts = np.unique(pred_labels[ref_label == ref_labels], return_counts=True)
        relative_overlaps_pred = [sub_pred_counts[i] / pred_dict[n] for i, n in enumerate(sub_pred_unique)]
        relative_overlaps_ref = [sub_pred_counts[i] / ref_dict[ref_label] for i, n in enumerate(sub_pred_unique)]
        relative_overlaps = np.c_[relative_overlaps_pred, relative_overlaps_ref]
        relative_overlaps_min = np.min(relative_overlaps, axis=1)
        pred_best_index = np.argsort(relative_overlaps_min)[::-1]
        associated_predictions[ref_label] = sub_pred_unique[pred_best_index]
        associated_overlaps[ref_label] = relative_overlaps[pred_best_index]
    if return_overlaps: return associated_predictions, associated_overlaps
    else: return associated_predictions

# --------------------------------------------------------------------------------
# Other stuff
# --------------------------------------------------------------------------------


# backwards compat... remove this in the future
def sanitize_anndata(adata):
    adata._sanitize()


def view_to_actual(adata):
    if adata.is_view:
        warnings.warn(
            "Revieved a view of an AnnData. Making a copy.",
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
    return ret[n - 1:] / n


# --------------------------------------------------------------------------------
# Deal with tool parameters
# --------------------------------------------------------------------------------


def update_params(
    old_params: Mapping[str, Any],
    new_params: Mapping[str, Any],
    check=False,
) -> Dict[str, Any]:
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
                raise ValueError('\'' + key
                                 + '\' is not a valid parameter key, '
                                 + 'consider one of \n'
                                 + str(list(old_params.keys())))
            if val is not None:
                updated_params[key] = val
    return updated_params


# --------------------------------------------------------------------------------
# Others
# --------------------------------------------------------------------------------


def select_groups(adata, groups_order_subset='all', key='groups'):
    """Get subset of groups in adata.obs[key].
    """
    groups_order = adata.obs[key].cat.categories
    if key + '_masks' in adata.uns:
        groups_masks = adata.uns[key + '_masks']
    else:
        groups_masks = np.zeros((len(adata.obs[key].cat.categories),
                                 adata.obs[key].values.size), dtype=bool)
        for iname, name in enumerate(adata.obs[key].cat.categories):
            # if the name is not found, fallback to index retrieval
            if adata.obs[key].cat.categories[iname] in adata.obs[key].values:
                mask = adata.obs[key].cat.categories[iname] == adata.obs[key].values
            else:
                mask = str(iname) == adata.obs[key].values
            groups_masks[iname] = mask
    groups_ids = list(range(len(groups_order)))
    if groups_order_subset != 'all':
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(
                np.where(adata.obs[key].cat.categories.values == name)[0][0])
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(
                np.in1d(np.arange(len(adata.obs[key].cat.categories)).astype(str),
                                          np.array(groups_order_subset)))[0]
        if len(groups_ids) == 0:
            logg.debug(
                f'{np.array(groups_order_subset)} invalid! specify valid '
                f'groups_order (or indices) from {adata.obs[key].cat.categories}',
            )
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks


def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    """Get full tracebacks when warning is raised by setting

    warnings.showwarning = warn_with_traceback

    See also
    --------
    http://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    import traceback
    traceback.print_stack()
    log = file if hasattr(file, 'write') else sys.stderr
    settings.write(warnings.formatwarning(message, category, filename, lineno, line))


def subsample(
    X: np.ndarray,
    subsample: int = 1,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
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
            raise ValueError(f'Invalid seed value < 0: {seed}')
        n = int(X.shape[0]/subsample)
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X, n=n)
    logg.debug(f'... subsampled to {n} of {X.shape[0]} data points')
    return Xsampled, rows


def subsample_n(
    X: np.ndarray, n: int = 0, seed: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
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
        raise ValueError('n must be greater 0')
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0], size=n, replace=False)
    Xsampled = X[rows]
    return Xsampled, rows


def check_presence_download(filename: Path, backup_url):
    """Check if file is present otherwise download."""
    if not filename.is_file():
        from .readwrite import _download
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
