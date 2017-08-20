# Author: F. Alex Wolf (http://falexwolf.de)
"""Utility functions and classes
"""

from collections import namedtuple
import numpy as np
from natsort import natsorted
from . import settings
from . import logging as logg


# --------------------------------------------------------------------------------
# Deal with stuff
# --------------------------------------------------------------------------------


def compute_minimum_spanning_tree(adjacency):
    from scipy.sparse.csgraph import minimum_spanning_tree
    return minimum_spanning_tree(adjacency)


def get_graph_tool_from_adjacency(adjacency, directed=None):
    """Get graph_tool graph from adjacency matrix."""
    import graph_tool as gt
    adjacency_edge_list = adjacency
    if not directed:
        from scipy.sparse import tril
        adjacency_edge_list = tril(adjacency)
    g = gt.Graph(directed=directed)
    g.add_vertex(adjacency.shape[0])
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
    weights = np.array(weights)[0]
    g = ig.Graph(list(zip(sources, targets)),
                 directed=directed,
                 edge_attrs={'weight': weights})
    if g.vcount() != adjacency.shape[0]:
        logg.warn('The constructed graph has only {} nodes. '
                  'Your adjacency matrix contained redundant nodes.'
                  .format(g.vcount()))
    return g


def compute_association_matrix_of_groups(adata, prediction, reference,
                                         normalization='prediction',
                                         threshold=0.01, max_n_names=2):
    """Identify predicted groups with reference groups.

    Parameters
    ----------
    adata : AnnData
    prediction : str
        Field name of adata.smp.
    reference : str
        Field name of adata.smp.
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
    check_adata(adata)
    asso_names = []
    asso_matrix = []
    for ipred_group, pred_group in enumerate(adata.add[prediction + '_order']):
        if '?' in pred_group: pred_group = str(ipred_group)
        # starting from numpy version 1.13, subtractions of boolean arrays are deprecated
        mask_pred = adata.smp[prediction] == pred_group
        mask_pred_int = mask_pred.astype(np.int8)
        asso_matrix += [[]]
        for ref_group in adata.add[reference + '_order']:
            mask_ref = (adata.smp[reference] == ref_group).astype(np.int8)
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
        name_list_pred = [adata.add[reference + '_order'][i]
                          for i in np.argsort(asso_matrix[-1])[::-1]
                          if asso_matrix[-1][i] > threshold]
        asso_names += ['\n'.join(name_list_pred[:max_n_names])]
    Result = namedtuple('compute_association_matrix_of_groups',
                        ['asso_names', 'asso_matrix'])
    return Result(asso_names=asso_names, asso_matrix=np.array(asso_matrix))


identify_categories = compute_association_matrix_of_groups
"""Deprecated: For backwards compatibility. See compute_association_matrix_of_groups."""


def get_associated_colors_of_groups(reference_colors, asso_matrix):
    asso_colors = [{reference_colors[i_ref]: asso_matrix[i_pred, i_ref]
                    for i_ref in range(asso_matrix.shape[1])}
                   for i_pred in range(asso_matrix.shape[0])]
    return asso_colors


get_associated_colors = get_associated_colors_of_groups
"""Deprecated: For backwards compatibility. See get_associated_colors_of_groups."""


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


def compute_fraction_of_agreeing_paths(adata1, adata2,
                                       adjacency_key='aga_attachedness_absolute'):
    """Compute the fraction of agreeing paths between leafs (nodes with degree
    1), a measure for the topological similarity between graphs.

    Notes
    -----

    By increasing the verbosity to level 4 and 5, the paths that do not agree
    and the paths that agree are written to the output, respectively.
    """
    import networkx as nx
    g1 = nx.Graph(adata1.add[adjacency_key])
    g2 = nx.Graph(adata2.add[adjacency_key])
    leaf_nodes1 = [str(x) for x in g1.nodes_iter() if g1.degree(x) == 1]
    asso_groups1 = identify_groups(adata1.smp['aga_groups'], adata2.smp['aga_groups'])
    asso_groups2 = identify_groups(adata2.smp['aga_groups'], adata1.smp['aga_groups'])

    import itertools
    n_steps = 0
    n_agreeing_steps = 0
    n_paths = 0
    n_agreeing_paths = 0
    # loop over all pairs of leaf nodes in the reference adata1
    for (r, s) in itertools.combinations(leaf_nodes1, r=2):
        r2, s2 = asso_groups1[r][0], asso_groups1[s][0]
        logg.m('compare shortest paths between leafs ({}, {}) in graph1 and ({}, {}) in graph2:'
               .format(r, s, r2, s2), v=4, no_indent=True)
        path1 = [str(x) for x in nx.shortest_path(g1, int(r), int(s))]
        path2 = [str(x) for x in nx.shortest_path(g2, int(r2), int(s2))]
        if len(path1) >= len(path2):
            path_mapped = [asso_groups1[l] for l in path1]
            path_compare = path2
            path_compare_id = 2
        else:
            path_mapped = [asso_groups2[l] for l in path2]
            path_compare = path1
            path_compare_id = 1
        n_agreeing_steps_path = 0
        ip_progress = 0
        for il, l in enumerate(path_compare[:-1]):
            for ip, p in enumerate(path_mapped):
                if ip >= ip_progress and l in p:
                    # check whether we can find the step forward of path_compare in path_mapped
                    if (ip + 1 < len(path_mapped)
                        and
                        path_compare[il + 1] in path_mapped[ip + 1]):
                        # make sure that a step backward leads us to the same value of l
                        # in case we "jumped"
                        logg.m('found matching step ({} -> {}) at position {} in path{} and position {} in path_mapped'
                               .format(l, path_compare[il + 1], il, path_compare_id, ip), v=6)
                        consistent_history = True
                        for iip in range(ip, ip_progress, -1):
                            if l not in path_mapped[iip - 1]:
                                consistent_history = False
                        if consistent_history:
                            # here, we take one step further back (ip_progress - 1); it's implied that this
                            # was ok in the previous step
                            logg.m('    step(s) backward to position(s) {} in path_mapped are fine, too: valid step'
                                   .format(list(range(ip - 1, ip_progress - 2, -1))), v=6)
                            n_agreeing_steps_path += 1
                            ip_progress = ip + 1
                            break
        n_agreeing_steps += n_agreeing_steps_path
        n_steps_path = len(path_compare) - 1
        n_steps += n_steps_path
        n_paths += 1

        if n_agreeing_steps_path == n_steps_path:
            n_agreeing_paths += 1
        logg.m('      path1 = {},\n'
               'path_mapped = {},\n'
               '      path2 = {},\n'
               '-> n_agreeing_steps = {} / n_steps = {}.'
               .format(path1, [list(p) for p in path_mapped], path2,
                       n_agreeing_steps_path, n_steps_path), v=5, no_indent=True)
    Result = namedtuple('compute_fraction_of_agreeing_paths_result',
                        ['frac_steps', 'n_steps', 'frac_paths', 'n_paths'])
    return Result(frac_steps=n_agreeing_steps/n_steps,
                  n_steps=n_steps,
                  frac_paths=n_agreeing_paths/n_paths,
                  n_paths=n_paths)



def remove_repetitions_from_list(l):
    return [l[0]] + [e for (i, e) in enumerate(l[1:]) if l[i] != e]


def plot_category_association(adata, prediction, reference, asso_matrix):
    pl.figure(figsize=(5, 5))
    pl.imshow(np.array(asso_matrix)[:], shape=(12, 4))
    pl.xticks(range(len(adata.add[reference + '_order'])), adata.add[reference + '_order'], rotation='vertical')
    pl.yticks(range(len(adata.add[prediction + '_order'])), adata.add[prediction + '_order'])
    pl.colorbar()


def unique_categories(categories):
    """Pass array-like categories, return sorted cleaned unique categories."""
    categories = np.unique(categories)
    categories = np.setdiff1d(categories, np.array(settings._ignore_categories))
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


_howto_specify_subgroups = '''sample annotation in adata only consists of sample names
--> you can provide additional annotation by setting, for example,
    adata.smp['groups'] = ['A', 'B', 'A', ... ]
    adata.smp['time'] = [0.1, 0.2, 0.7, ... ]'''


def check_adata(adata, verbosity=-3):
    """Do sanity checks on adata object.

    Checks whether adata contains annotation.
    """
    if len(adata.smp_keys()) == 0:
        settings.m(1-verbosity, _howto_specify_subgroups)
    else:
        if len(adata.smp_keys()) > 0 and settings.verbosity > 1-verbosity:
            info = 'sample annotation: '
        for ismp, smp in enumerate(adata.smp_keys()):
            # ordered unique categories for categorical annotation
            if not smp + '_order' in adata.add and adata.smp[smp].dtype.char in {'U', 'S'}:
                adata.add[smp + '_order'] = unique_categories(adata.smp[smp])
            if settings.verbosity > 1-verbosity:
                info += '"' + smp + '" = '
                if adata.smp[smp].dtype.char in {'U', 'S'}:
                    ann_info = str(adata.add[smp + '_order'])
                    if len(adata.add[smp + '_order']) > 7:
                        ann_info = (str(adata.add[smp + '_order'][0:3]).replace(']', '')
                                    + ' ...'
                                    + str(adata.add[smp + '_order'][-2:]).replace('[', ''))
                    info += ann_info
                else:
                    info += 'continuous'
                if ismp < len(adata.smp_keys())-1:
                    info += ', '
        if len(adata.smp_keys()) > 0 and settings.verbosity > 1-verbosity:
            settings.m(1-verbosity, info)
    return adata


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


def select_groups(adata, groups_order_subset='all', smp='groups'):
    """Get subset of groups in adata.smp[smp].
    """
    groups_order = adata.add[smp + '_order']
    if smp + '_masks' in adata.add:
        groups_masks = adata.add[smp + '_masks']
    else:
        groups_masks = np.zeros((len(adata.add[smp + '_order']),
                                 adata.smp[smp].size), dtype=bool)
        for iname, name in enumerate(adata.add[smp + '_order']):
            # if the name is not found, fallback to index retrieval
            if adata.add[smp + '_order'][iname] in adata.smp[smp]:
                mask = adata.add[smp + '_order'][iname] == adata.smp[smp]
            else:
                mask = str(iname) == adata.smp[smp]
            groups_masks[iname] = mask
    groups_ids = list(range(len(groups_order)))
    if groups_order_subset != 'all':
        # get list from string
        if isinstance(groups_order_subset, str):
            groups_order_subset = groups_order_subset.split(',')
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(np.where(adata.add[smp + '_order'] == name)[0][0])
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(np.in1d(np.arange(len(adata.add[smp + '_order'])).astype(str),
                                          np.array(groups_order_subset)))[0]
        if len(groups_ids) == 0:
            logg.m(np.array(groups_order_subset),
                   'invalid! specify valid groups_order (or indices) one of',
                   adata.add[smp + '_order'])
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.add[smp + '_order'][groups_ids]
    else:
        groups_order_subset = groups_order
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

    Note
    ----
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
    import warnings
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


def comp_distance(X, metric='euclidean'):
    """Compute distance matrix for data array X

    Parameters
    ----------
    X : np.ndarray
        Data array (rows store samples, columns store variables).
    metric : string
        For example 'euclidean', 'sqeuclidean', see sp.spatial.distance.pdist.

    Returns
    -------
    D : np.ndarray
        Distance matrix.
    """
    from scipy.spatial import distance
    return distance.squareform(distance.pdist(X, metric=metric))


def comp_sqeuclidean_distance_using_matrix_mult(X, Y):
    """Compute distance matrix for data array X

    Use matrix multiplication as in sklearn.

    Parameters
    ----------
    X : np.ndarray
        Data array (rows store samples, columns store variables).
    metric : string
        For example 'euclidean', 'sqeuclidean', see sp.spatial.distance.pdist.

    Returns
    -------
    D : np.ndarray
        Distance matrix.
    """
    XX = np.einsum('ij,ij->i', X, X)[:, np.newaxis]
    if X is Y:
        YY = XX
    else:
        YY = np.einsum('ij,ij->i', Y, Y)[:, np.newaxis]
    distances = np.dot(X, Y.T)
    distances *= -2
    distances += XX
    distances += YY.T
    np.maximum(distances, 0, out=distances)
    if X is Y:
        distances.flat[::distances.shape[0] + 1] = 0.
    return distances


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
