# Author: F. Alex Wolf (http://falexwolf.de)
"""Utility functions and classes
"""

import numpy as np
from natsort import natsorted
from . import settings as sett
from . import logging as logg

# --------------------------------------------------------------------------------
# Deal with stuff
# --------------------------------------------------------------------------------


def get_igraph_from_adjacency(adjacency, directed=None):
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


def identify_categories(adata, prediction, reference,
                        normalization='prediction',
                        threshold=0.01, max_n_names=2):
    """Identify predicted categories with reference.

    Parameters
    ----------
    adata : AnnData
    prediction : str
        smp_key of adata
    reference : str
        smp_key of adata
    maximum : int or None, optional (default: 2)
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
    for ipred_group, pred_group in enumerate(adata.add[prediction + '_names']):
        if '?' in pred_group: pred_group = str(ipred_group)
        # starting from numpy version 1.13, subtractions of boolean arrays are deprecated
        mask_pred = adata.smp[prediction] == pred_group
        mask_pred_int = mask_pred.astype(np.int8)
        asso_matrix += [[]]
        for ref_group in adata.add[reference + '_names']:
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
        name_list_pred = [adata.add[reference + '_names'][i]
                          for i in np.argsort(asso_matrix[-1])[::-1]
                          if asso_matrix[-1][i] > threshold]
        asso_names += ['\n'.join(name_list_pred[:max_n_names])]
    return asso_names, np.array(asso_matrix)


def get_associated_colors(reference_colors, asso_matrix):
    asso_colors = [{reference_colors[i_ref]: asso_matrix[i_pred, i_ref]
                    for i_ref in range(asso_matrix.shape[1])}
                   for i_pred in range(asso_matrix.shape[0])]
    return asso_colors


def plot_category_association(adata, prediction, reference, asso_matrix):
    pl.figure(figsize=(5, 5))
    pl.imshow(np.array(asso_matrix)[:], shape=(12, 4))
    pl.xticks(range(len(adata.add[reference + '_names'])), adata.add[reference + '_names'], rotation='vertical')
    pl.yticks(range(len(adata.add[prediction + '_names'])), adata.add[prediction + '_names'])
    pl.colorbar()


def unique_categories(categories):
    """Pass array-like categories, return sorted cleaned unique categories."""
    categories = np.unique(categories)
    categories = np.setdiff1d(categories, np.array(sett._ignore_categories))
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
        sett.m(1-verbosity, _howto_specify_subgroups)
    else:
        if len(adata.smp_keys()) > 0 and sett.verbosity > 1-verbosity:
            info = 'sample annotation: '
        for ismp, smp in enumerate(adata.smp_keys()):
            # ordered unique categories for categorical annotation
            if not smp + '_names' in adata.add and adata.smp[smp].dtype.char in {'U', 'S'}:
                adata.add[smp + '_names'] = unique_categories(adata.smp[smp])
            if sett.verbosity > 1-verbosity:
                info += '"' + smp + '" = '
                if adata.smp[smp].dtype.char in {'U', 'S'}:
                    ann_info = str(adata.add[smp + '_names'])
                    if len(adata.add[smp + '_names']) > 7:
                        ann_info = (str(adata.add[smp + '_names'][0:3]).replace(']', '')
                                    + ' ...'
                                    + str(adata.add[smp + '_names'][-2:]).replace('[', ''))
                    info += ann_info
                else:
                    info += 'continuous'
                if ismp < len(adata.smp_keys())-1:
                    info += ', '
        if len(adata.smp_keys()) > 0 and sett.verbosity > 1-verbosity:
            sett.m(1-verbosity, info)
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
    p = default_tool_argparser(sc.help(toolkey, string=True), example_parameters)
    if tool_add_args is None:
        p = add_args(p)
    else:
        p = tool_add_args(p)
    args = vars(p.parse_args())
    args = sett.process_args(args)
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


def select_groups(adata, groups_names_subset='all', smp='groups'):
    """Get subset of groups in adata.smp[smp].
    """
    groups_names = adata.add[smp + '_names']
    if smp + '_masks' in adata.add:
        groups_masks = adata.add[smp + '_masks']
    else:
        groups_masks = np.zeros((len(adata.add[smp + '_names']),
                                 adata.smp[smp].size), dtype=bool)
        for iname, name in enumerate(adata.add[smp + '_names']):
            # if the name is not found, fallback to index retrieval
            if adata.add[smp + '_names'][iname] in adata.smp[smp]:
                mask = adata.add[smp + '_names'][iname] == adata.smp[smp]
            else:
                mask = str(iname) == adata.smp[smp]
            groups_masks[iname] = mask
    groups_ids = list(range(len(groups_names)))
    if groups_names_subset != 'all':
        # get list from string
        if isinstance(groups_names_subset, str):
            groups_names_subset = groups_names_subset.split(',')
        groups_ids = []
        for name in groups_names_subset:
            groups_ids.append(np.where(adata.add[smp + '_names'] == name)[0][0])
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(np.in1d(np.arange(len(adata.add[smp + '_names'])).astype(str),
                                          np.array(groups_names_subset)))[0]
        if len(groups_ids) == 0:
            logg.m(np.array(groups_names_subset),
                   'invalid! specify valid groups_names (or indices) one of',
                   adata.add[smp + '_names'])
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_names_subset = adata.add[smp + '_names'][groups_ids]
    else:
        groups_names_subset = groups_names
    return groups_names_subset, groups_masks


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
    sett.write(warnings.formatwarning(message, category, filename, lineno, line))


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
    sett.mt(0, 'clustered matrix')
    if False:
        pl.matshow(Mclus)
        pl.colorbar()
    return Mclus, indices
