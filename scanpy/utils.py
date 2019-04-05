"""Utility functions and classes
"""

import sys
import inspect
from weakref import WeakSet
from collections import namedtuple
from functools import partial, wraps
from types import ModuleType, MethodType
from typing import Union, Callable, Optional

import numpy as np
import scipy.sparse
from natsort import natsorted
from textwrap import dedent
from pandas.api.types import CategoricalDtype

from ._settings import settings
from . import logging as logg
import warnings

EPS = 1e-15


def check_versions():
    from distutils.version import LooseVersion

    if sys.version_info < (3, 6):
        warnings.warn('Scanpy prefers Python 3.6 or higher. '
                      'Currently, Python 3.5 leads to a bug in `tl.marker_gene_overlap` '
                      'and we might stop supporting it in the future.')

    import anndata
    # NOTE: pytest does not correctly retrieve anndata's version? why?
    #       use the following hack...
    if anndata.__version__ != '0+unknown':
        if anndata.__version__ < LooseVersion('0.6.10'):
            from . import __version__
            raise ImportError('Scanpy {} needs anndata version >=0.6.10, not {}.\n'
                              'Run `pip install anndata -U --no-deps`.'
                              .format(__version__, anndata.__version__))


def getdoc(c_or_f: Union[Callable, type]) -> Optional[str]:
    if getattr(c_or_f, '__doc__', None) is None:
        return None
    doc = inspect.getdoc(c_or_f)
    if isinstance(c_or_f, type) and hasattr(c_or_f, '__init__'):
        sig = inspect.signature(c_or_f.__init__)
    else:
        sig = inspect.signature(c_or_f)

    def type_doc(name: str):
        param = sig.parameters[name]  # type: inspect.Parameter
        cls = getattr(param.annotation, '__qualname__', repr(param.annotation))
        if param.default is not param.empty:
            return '{}, optional (default: {!r})'.format(cls, param.default)
        else:
            return cls

    return '\n'.join(
        '{} : {}'.format(line, type_doc(line)) if line.strip() in sig.parameters else line
        for line in doc.split('\n')
    )


def deprecated_arg_names(arg_mapping):
    """
    Decorator which marks a functions keyword arguments as deprecated. It will
    result in a warning being emitted when the deprecated keyword argument is
    used, and the function being called with the new argument.

    Parameters
    ----------
    arg_mapping : dict[str, str]
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
                        "Keyword argument '{0}' has been deprecated in favour "
                        "of '{1}'. '{0}' will be removed in a future version."
                        .format(old, new),
                        category=DeprecationWarning,
                        stacklevel=2,
                    )
                    val = kwargs.pop(old)
                    kwargs[new] = val
            warnings.simplefilter(
                'default', DeprecationWarning)  # reset filter
            return func(*args, **kwargs)
        return func_wrapper
    return decorator


def descend_classes_and_funcs(mod: ModuleType, root: str, encountered=None):
    if encountered is None:
        encountered = WeakSet()
    for obj in vars(mod).values():
        if not getattr(obj, '__module__', getattr(obj, '__qualname__', getattr(obj, '__name__', ''))).startswith(root):
            continue
        if isinstance(obj, Callable) and not isinstance(obj, MethodType):
            yield obj
            if isinstance(obj, type):
                yield from (m for m in vars(obj).values() if isinstance(m, Callable))
        elif isinstance(obj, ModuleType) and obj not in encountered:
            encountered.add(obj)
            yield from descend_classes_and_funcs(obj, root, encountered)


def annotate_doc_types(mod: ModuleType, root: str):
    for c_or_f in descend_classes_and_funcs(mod, root):
        c_or_f.getdoc = partial(getdoc, c_or_f)


def doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """
    def dec(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwds)
        return obj
    return dec


def merge_groups(adata, key, map_groups, key_added=None, map_colors=None):
    """
    Parameters
    ----------
    map_colors : `dict`
        Dict with color specification for new groups that have no corresponding
        old group.
    """
    if key_added is None:
        key_added = key + '_merged'
    adata.obs[key_added] = adata.obs[key].map(
        map_groups).astype(CategoricalDtype())
    old_categories = adata.obs[key].cat.categories
    new_categories = adata.obs[key_added].cat.categories
    # map_colors is passed
    if map_colors is not None:
        old_colors = None
        if key + '_colors' in adata.uns:
            old_colors = adata.uns[key + '_colors']
        new_colors = []
        for group in adata.obs[key_added].cat.categories:
            if group in map_colors:
                new_colors.append(map_colors[group])
            elif group in old_categories and old_colors is not None:
                new_colors.append(old_colors[old_categories.get_loc(group)])
            else:
                raise ValueError('You didn\'t specify a color for {}.'
                                 .format(group))
        adata.uns[key_added + '_colors'] = new_colors
    # map_colors is not passed
    elif key + '_colors' in adata.uns:
        old_colors = adata.uns[key + '_colors']
        inverse_map_groups = {g: [] for g in new_categories}
        for old_group in old_categories:
            inverse_map_groups[map_groups[old_group]].append(old_group)
        new_colors = []
        for group in new_categories:
            # take the largest of the old groups
            old_group = adata.obs[key][adata.obs[key].isin(
                inverse_map_groups[group])].value_counts().index[0]
            new_colors.append(old_colors[old_categories.get_loc(old_group)])
        adata.uns[key_added + '_colors'] = new_colors


# --------------------------------------------------------------------------------
# Graph stuff
# --------------------------------------------------------------------------------


def cross_entropy_neighbors_in_rep(adata, use_rep, n_points=3):
    """Compare neighborhood graph of representation based on cross entropy.

    `n_points` denotes the number of points to add as highlight annotation.

    Returns
    -------
    The cross entropy and the geodesic-distance-weighted cross entropy as
    ``entropy, geo_entropy_d, geo_entropy_o``.

    Adds the most overlapping or disconnected points as annotation to `adata`.
    """
    # see below why we need this
    if 'X_diffmap' not in adata.obsm.keys():
        raise ValueError('Run `tl.diffmap` on `adata`, first.')

    adata_ref = adata  # simple renaming, don't need copy here
    adata_cmp = adata.copy()
    n_neighbors = adata_ref.uns['neighbors']['params']['n_neighbors']
    from .neighbors import neighbors
    neighbors(adata_cmp, n_neighbors=n_neighbors, use_rep=use_rep)
    from .tools.diffmap import diffmap
    diffmap(adata_cmp)

    graph_ref = adata_ref.uns['neighbors']['connectivities']
    graph_cmp = adata_cmp.uns['neighbors']['connectivities']

    graph_ref = graph_ref.tocoo()  # makes a copy
    graph_cmp = graph_cmp.tocoo()

    edgeset_ref = {e for e in zip(graph_ref.row, graph_ref.col)}
    edgeset_cmp = {e for e in zip(graph_cmp.row, graph_cmp.col)}
    edgeset_union = list(edgeset_ref.union(edgeset_cmp))

    edgeset_union_indices = tuple(zip(*edgeset_union))
    edgeset_union_indices = (np.array(edgeset_union_indices[0]), np.array(edgeset_union_indices[1]))

    n_edges_ref = len(graph_ref.nonzero()[0])
    n_edges_cmp = len(graph_cmp.nonzero()[0])
    n_edges_union = len(edgeset_union)
    logg.msg(
        '... n_edges_ref', n_edges_ref,
        'n_edges_cmp', n_edges_cmp,
        'n_edges_union', n_edges_union)

    graph_ref = graph_ref.tocsr()  # need a copy of the csr graph anyways
    graph_cmp = graph_cmp.tocsr()

    p_ref = graph_ref[edgeset_union_indices].A1
    p_cmp = graph_cmp[edgeset_union_indices].A1

    # the following is how one compares it to log_loss form sklearn
    # p_ref[p_ref.nonzero()] = 1
    # from sklearn.metrics import log_loss
    # print(log_loss(p_ref, p_cmp))
    p_cmp = np.clip(p_cmp, EPS, 1-EPS)
    ratio = np.clip(p_ref / p_cmp, EPS, None)
    ratio_1m = np.clip((1 - p_ref) / (1 - p_cmp), EPS, None)

    entropy = np.sum(p_ref * np.log(ratio) + (1-p_ref) * np.log(ratio_1m))

    n_edges_fully_connected = (graph_ref.shape[0]**2 - graph_ref.shape[0])
    entropy /= n_edges_fully_connected

    fraction_edges = n_edges_ref / n_edges_fully_connected
    naive_entropy = (fraction_edges * np.log(1./fraction_edges)
                     + (1-fraction_edges) * np.log(1./(1-fraction_edges)))
    logg.msg('cross entropy of naive sparse prediction {:.3e}'.format(naive_entropy))
    logg.msg('cross entropy of random prediction {:.3e}'.format(-np.log(0.5)))
    logg.info('cross entropy {:.3e}'.format(entropy))

    # for manifold analysis, restrict to largest connected component in
    # reference
    # now that we clip at a quite high value below, this might not even be
    # necessary
    n_components, labels = scipy.sparse.csgraph.connected_components(graph_ref)
    largest_component = np.arange(graph_ref.shape[0], dtype=int)
    if n_components > 1:
        component_sizes = np.bincount(labels)
        logg.msg('largest component has size', component_sizes.max())
        largest_component = np.where(
            component_sizes == component_sizes.max())[0][0]
        graph_ref_red = graph_ref.tocsr()[labels == largest_component, :]
        graph_ref_red = graph_ref_red.tocsc()[:, labels == largest_component]
        graph_ref_red = graph_ref_red.tocoo()
        graph_cmp_red = graph_cmp.tocsr()[labels == largest_component, :]
        graph_cmp_red = graph_cmp_red.tocsc()[:, labels == largest_component]
        graph_cmp_red = graph_cmp_red.tocoo()
        edgeset_ref_red = {e for e in zip(graph_ref_red.row, graph_ref_red.col)}
        edgeset_cmp_red = {e for e in zip(graph_cmp_red.row, graph_cmp_red.col)}
        edgeset_union_red = edgeset_ref_red.union(edgeset_cmp_red)
        map_indices = np.where(labels == largest_component)[0]
        edgeset_union_red = {
            (map_indices[i], map_indices[j]) for (i, j) in edgeset_union_red}

    from .neighbors import Neighbors
    neigh_ref = Neighbors(adata_ref)
    dist_ref = neigh_ref.distances_dpt  # we expect 'X_diffmap' to be already present

    neigh_cmp = Neighbors(adata_cmp)
    dist_cmp = neigh_cmp.distances_dpt

    d_cmp = np.zeros_like(p_ref)
    d_ref = np.zeros_like(p_ref)
    for i, e in enumerate(edgeset_union):
        # skip contributions that are not in the largest component
        if n_components > 1 and e not in edgeset_union_red:
            continue
        d_cmp[i] = dist_cmp[e]
        d_ref[i] = dist_ref[e]

    MAX_DIST = 1000
    d_cmp = np.clip(d_cmp, 0.1, MAX_DIST)  # we don't want to measure collapsing clusters
    d_ref = np.clip(d_ref, 0.1, MAX_DIST)

    weights = np.array(d_cmp / d_ref)            # disconnected regions
    weights_overlap = np.array(d_ref / d_cmp)    # overlapping regions

    # the following is just for annotation of figures
    if 'highlights' not in adata_ref.uns:
        adata_ref.uns['highlights'] = {}
    else:
        # remove old disconnected and overlapping points
        new_highlights = {}
        for k, v in adata_ref.uns['highlights'].items():
            if v != 'O' and v not in {'D0', 'D1', 'D2', 'D3', 'D4'}:
                new_highlights[k] = v
        adata_ref.uns['highlights'] = new_highlights

    # points that are maximally disconnected
    max_weights = np.argpartition(weights, kth=-n_points)[-n_points:]
    points = list(edgeset_union_indices[0][max_weights])
    points2 = list(edgeset_union_indices[1][max_weights])
    found_disconnected_points = False
    for ip, p in enumerate(points):
        if d_cmp[max_weights][ip] == MAX_DIST:
            adata_ref.uns['highlights'][p] = 'D' + str(ip)
            adata_ref.uns['highlights'][points2[ip]] = 'D' + str(ip)
            found_disconnected_points = True
    if found_disconnected_points:
        logg.msg('most disconnected points', points)
        logg.msg('    with weights', weights[max_weights].round(1))

    max_weights = np.argpartition(
        weights_overlap, kth=-n_points)[-n_points:]
    points = list(edgeset_union_indices[0][max_weights])
    for p in points:
        adata_ref.uns['highlights'][p] = 'O'
    logg.msg('most overlapping points', points)
    logg.msg('    with weights', weights_overlap[max_weights].round(1))
    logg.msg('    with d_rep', d_cmp[max_weights].round(1))
    logg.msg('    with d_ref', d_ref[max_weights].round(1))

    geo_entropy_d = np.sum(weights * p_ref * np.log(ratio))
    geo_entropy_o = np.sum(weights_overlap * (1-p_ref) * np.log(ratio_1m))

    geo_entropy_d /= n_edges_fully_connected
    geo_entropy_o /= n_edges_fully_connected

    logg.info('geodesic cross entropy {:.3e}'.format(geo_entropy_d + geo_entropy_o))
    return entropy, geo_entropy_d, geo_entropy_o


def get_graph_tool_from_adjacency(adjacency, directed=None):
    """Get graph_tool graph from adjacency matrix."""
    import graph_tool as gt
    adjacency_edge_list = adjacency
    if not directed:
        from scipy.sparse import tril
        adjacency_edge_list = tril(adjacency)
    g = gt.Graph(directed=directed)
    g.add_vertex(adjacency.shape[0])  # this adds adjacency.shap[0] vertices
    g.add_edge_list(np.transpose(adjacency_edge_list.nonzero()))
    weights = g.new_edge_property('double')
    for e in g.edges():
        # graph_tool uses the following convention,
        # which is opposite to the rest of scanpy
        weights[e] = adjacency[int(e.source()), int(e.target())]
    g.edge_properties['weight'] = weights
    return g


def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shap[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warn('The constructed graph has only {} nodes. '
                  'Your adjacency matrix contained redundant nodes.'
                  .format(g.vcount()))
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


def compute_association_matrix_of_groups(adata, prediction, reference,
                                         normalization='prediction',
                                         threshold=0.01, max_n_names=2):
    """Compute overlaps between groups.

    See ``identify_groups`` for identifying the groups.

    Parameters
    ----------
    adata : AnnData
    prediction : str
        Field name of adata.obs.
    reference : str
        Field name of adata.obs.
    normalization : {'prediction', 'reference'}
        Whether to normalize with respect to the predicted groups or the
        reference groups.
    threshold : float, optional (default: 0.01)
        Do not consider associations whose overlap is below this fraction.
    max_n_names : int or None, optional (default: 2)
        Control how many reference names you want to be associated with per
        predicted name. Set to `None`, if you want all.

    Returns
    -------
    Tuple of
    asso_names : list of associated reference names (`max_n_names` for each
        predicted name)
    asso_matrix : matrix where rows correspond to the predicted labels and
        columns to the reference labels, entries are proportional to degree of
        association
    """
    if normalization not in {'prediction', 'reference'}:
        raise ValueError('`normalization` needs to be either "prediction" or "reference".')
    sanitize_anndata(adata)
    cats = adata.obs[reference].cat.categories
    for cat in cats:
        if cat in settings.categories_to_ignore:
            logg.info('Ignoring category \'{}\' '
                      'as it\'s in `settings.categories_to_ignore`.'
                      .format(cat))
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
    asso_colors = [{reference_colors[i_ref]: asso_matrix[i_pred, i_ref]
                    for i_ref in range(asso_matrix.shape[1])}
                   for i_pred in range(asso_matrix.shape[0])]
    return asso_colors


def compute_group_overlap_score(ref_labels, pred_labels,
                                threshold_overlap_pred=0.5,
                                threshold_overlap_ref=0.5):
    """How well do the pred_labels explain the ref_labels?

    A predicted cluster explains a reference cluster if it is contained within the reference
    cluster with at least 50% (threshold_overlap_pred) of its points and these correspond
    to at least 50% (threshold_overlap_ref) of the reference cluster.
    """
    ref_unique, ref_counts = np.unique(ref_labels, return_counts=True)
    ref_dict = dict(zip(ref_unique, ref_counts))
    pred_unique, pred_counts = np.unique(pred_labels, return_counts=True)
    pred_dict = dict(zip(pred_unique, pred_counts))
    summary = []
    for true in ref_unique:
        sub_pred_unique, sub_pred_counts = np.unique(pred_labels[true == ref_labels], return_counts=True)
        relative_overlaps_pred = [sub_pred_counts[i] / pred_dict[n] for i, n in enumerate(sub_pred_unique)]
        relative_overlaps_ref = [sub_pred_counts[i] / ref_dict[true] for i, n in enumerate(sub_pred_unique)]
        pred_best_index = np.argmax(relative_overlaps_pred)
        summary.append(1 if (relative_overlaps_pred[pred_best_index] >= threshold_overlap_pred and
                             relative_overlaps_ref[pred_best_index] >= threshold_overlap_ref)
                       else 0)
        # print(true, sub_pred_unique[pred_best_index], relative_overlaps_pred[pred_best_index],
        #       relative_overlaps_ref[pred_best_index], summary[-1])
    return sum(summary)/len(summary)


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


def remove_repetitions_from_list(l):
    return [l[0]] + [e for (i, e) in enumerate(l[1:]) if l[i] != e]


def plot_category_association(adata, prediction, reference, asso_matrix):
    pl.figure(figsize=(5, 5))
    pl.imshow(np.array(asso_matrix)[:], shape=(12, 4))
    pl.xticks(range(len(adata.uns[reference + '_order'])), adata.uns[reference + '_order'], rotation='vertical')
    pl.yticks(range(len(adata.uns[prediction + '_order'])), adata.uns[prediction + '_order'])
    pl.colorbar()


def unique_categories(categories):
    """Pass array-like categories, return sorted cleaned unique categories."""
    categories = np.unique(categories)
    categories = np.setdiff1d(categories, np.array(settings.categories_to_ignore))
    categories = np.array(natsorted(categories, key=lambda v: v.upper()))
    return categories


def fill_in_datakeys(example_parameters, dexdata):
    """Update the 'examples dictionary' _examples.example_parameters.

    If a datakey (key in 'datafile dictionary') is not present in the 'examples
    dictionary' it is used to initialize an entry with that key.

    If not specified otherwise, any 'exkey' (key in 'examples dictionary') is
    used as 'datakey'.
    """
    # default initialization of 'datakey' key with entries from data dictionary
    for exkey in example_parameters:
        if 'datakey' not in example_parameters[exkey]:
            if exkey in dexdata:
                example_parameters[exkey]['datakey'] = exkey
            else:
                example_parameters[exkey]['datakey'] = 'unspecified in dexdata'
    return example_parameters


# backwards compat... remove this in the future
def sanitize_anndata(adata):
    adata._sanitize()


def moving_average(a, n):
    """Moving average over one-dimensional array.

    Parameters
    ----------
    a : np.ndarray
        One-dimensional array.
    n : int
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


def update_params(old_params, new_params, check=False):
    """Update old_params with new_params.

    If check==False, this merely adds and overwrites the content of old_params.

    If check==True, this only allows updating of parameters that are already
    present in old_params.

    Parameters
    ----------
    old_params : dict
    new_params : dict
    check : bool, optional (default: False)

    Returns
    -------
    updated_params : dict
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
# Command-line argument reading and processing
# --------------------------------------------------------------------------------


def read_args_tool(toolkey, example_parameters, tool_add_args=None):
    """Read args for single tool.
    """
    import scanpy as sc
    p = default_tool_argparser(help(toolkey), example_parameters)
    if tool_add_args is None:
        p = add_args(p)
    else:
        p = tool_add_args(p)
    args = vars(p.parse_args())
    args = settings.process_args(args)
    return args


def default_tool_argparser(description, example_parameters):
    """Create default parser for single tools.
    """
    import argparse
    epilog = '\n'
    for k, v in sorted(example_parameters.items()):
        epilog += '  ' + k + '\n'
    p = argparse.ArgumentParser(
        description=description,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=('available values for examples (exkey):'+epilog))
    return p


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
            logg.m(np.array(groups_order_subset),
                   'invalid! specify valid groups_order (or indices) one of',
                   adata.obs[key].cat.categories)
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks


def pretty_dict_string(d, indent=0):
    """Pretty output of nested dictionaries.
    """
    s = ''
    for key, value in sorted(d.items()):
        s += '    ' * indent + str(key)
        if isinstance(value, dict):
             s += '\n' + pretty_dict_string(value, indent+1)
        else:
             s += '=' + str(value) + '\n'
    return s


def markdown_dict_string(d):
    """Markdown output that can be pasted in the examples/README.md.
    """
    # sort experimental data from simulated data
    from collections import OrderedDict
    types = OrderedDict()
    for key, value in sorted(d.items()):
        if 'type' in value:
            if value['type'] not in types:
                types[value['type']] = []
            types[value['type']].append(key)
        else:
            print(key, 'does not define data type!')
    # format output
    s = ''
    for type in ['scRNAseq', 'scqPCR', 'bulk', 'simulated']:
        s += '\nExamples using ' + type + ' data.\n'
        for key in types[type]:
            value = d[key]
            s += '* [' + key + '](#' + key + ')'
            if 'ref' in value:
                if 'doi' in value:
                    link = 'http://dx.doi.org/' + value['doi']
                elif 'url' in value:
                    link = value['url']
                s += (' - [' +  value['ref'].replace('et al.','*et al.*')
                             + '](' + link +  ')')
            if 'title' in value:
                s += '   \n*' + value['title'] + '*'
            s += '\n'
    return s


def merge_dicts(*ds):
    """Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Notes
    -----
    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = ds[0]
    for d in ds[1:]:
        result.update(d)
    return result


def masks(list_of_index_lists, n):
    """Make an array in which rows store 1d mask arrays from list of index lists.

    Parameters
    ----------
    n : int
        Maximal index / number of samples.
    """
    # make a list of mask arrays, it's easier to store
    # as there is a hdf5 equivalent
    for il,l in enumerate(list_of_index_lists):
        mask = np.zeros(n,dtype=bool)
        mask[l] = True
        list_of_index_lists[il] = mask
    # convert to arrays
    masks = np.array(list_of_index_lists)
    return masks


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


def subsample(X, subsample=1, seed=0):
    """Subsample a fraction of 1/subsample samples from the rows of X.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    subsample : int
        1/subsample is the fraction of data sampled, n = X.shape[0]/subsample.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray
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
    if seed > 0:
        n = int(X.shape[0]/subsample)
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X, n=n)
    logg.m('... subsampled to', n, 'of', X.shape[0], 'data points')
    return Xsampled, rows


def subsample_n(X, n=0, seed=0):
    """Subsample n samples from rows of array.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray
        Indices of rows that are stored in Xsampled.
    """
    if n < 0:
        raise ValueError('n must be greater 0')
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0], size=n, replace=False)
    Xsampled = X[rows]
    return Xsampled, rows


def check_presence_download(filename, backup_url):
    """Check if file is present otherwise download."""
    import os
    filename = str(filename) #  Throws error for Path on 3.5
    if not os.path.exists(filename):
        from .readwrite import download_progress
        dr = os.path.dirname(filename)
        try:
            os.makedirs(dr)
        except FileExistsError:
            pass  # ignore if dir already exists
        from urllib.request import urlretrieve
        urlretrieve(backup_url, filename, reporthook=download_progress)


def hierarch_cluster(M):
    """Cluster matrix using hierarchical clustering.

    Parameters
    ----------
    M : np.ndarray
        Matrix, for example, distance matrix.

    Returns
    -------
    Mclus : np.ndarray
        Clustered matrix.
    indices : np.ndarray
        Indices used to cluster the matrix.
    """
    import scipy as sp
    import scipy.cluster
    link = sp.cluster.hierarchy.linkage(M)
    indices = sp.cluster.hierarchy.leaves_list(link)
    Mclus = np.array(M[:, indices])
    Mclus = Mclus[indices, :]
    if False:
        pl.matshow(Mclus)
        pl.colorbar()
    return Mclus, indices
