# Author: T. Callies
#         F. Alex Wolf (http://falexwolf.de)
"""Exporting to formats for other software.
"""

import numpy as np
import os
import json
from math import inf
from scipy.sparse import issparse
import logging as logg
from pandas.api.types import is_string_dtype, is_categorical
from .plotting.utils import add_colors_for_categorical_sample_annotation


def spring_project(adata, project_dir, use_genes=None, cell_groupings=None,
                   custom_color_tracks=None):
    """Exports to a SPRING project directory [Weinreb17]_

    See https://github.com/AllonKleinLab/SPRING or Weinreb17_ for details.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix: `adata.uns['data_graph_distance_local']` needs to
        be present.
    project_dir : `str`
        Path to a directory where SPRING readable files will be written. The
        directory does not have to exist before running this function.
    use_genes : `str` or `list`, optional (default: `None`)
        Select a subset of genes. If a `str`, looks for annotation in
        `adata.uns` useful to plot marker genes found with
        :func:`~scanpy.api.tl.rank_gene_groups`.
    cell_groupings : `str`, `list` of `str`, optional (default: `None`)
        Optional list of strings containing `adata.obs` key to grouping.
    custom_color_tracks : `str`, `list` of `str`, optional (default: `None`)
        Optional list of strings containing `adata.obs` key for continuous
        coloring.
    """

    gene_list = adata.var_names
    # We allow to include rank_genes annotation.
    if isinstance(use_genes, str):
        use_genes_list = []
        if use_genes in adata.uns:
            for rank in adata.uns[use_genes]:
                for groups in rank:
                    use_genes_list.append(groups)
            use_genes = use_genes_list
        else:
            # TODO: the following check seems fishy
            if use_genes not in adata.var_names:
                logg.warn(
                    '{} annotation not found. Call `rank_gene_groups` '
                    'or make sure gene names are in `adata.var_names`.'
                    .format(use_genes))
            else:
                # if a single gene contained in AnnData object
                use_genes_list = [use_genes]
        gene_list = use_genes_list
    elif isinstance(use_genes, (list, np.ndarray)):
        gene_list = list(use_genes)

    # File can be safed anywhere. However, for easy access via SPRING, safe it
    # somewhere in the spring directory
    os.system('mkdir ' + project_dir)
    if not project_dir[-1] == '/': project_dir += '/'

    if 'data_graph_distance_local' not in adata.uns:
        raise ValueError(
            'Run any tool that produces a data graph first, '
            'e.g. sc.tl.diffmap or sc.tl.louvain')
    # Note that output here will always be sparse
    D = adata.uns['data_graph_distance_local']
    k = adata.uns['data_graph_distance_local'][0].nonzero()[0].size
    edges = get_knn_edges_sparse(D, k)

    # write custom color tracks
    if isinstance(custom_color_tracks, str):
        custom_color_tracks = [custom_color_tracks]
    if custom_color_tracks is not None:
        custom_colors = {g: adata.obs[g] for g in custom_color_tracks}
    else:
        # write all annotation that's neither categorical or string
        custom_colors = {k: adata.obs[k] for k in adata.obs_keys()
                         if not (is_categorical(adata.obs[k])
                                 or is_string_dtype(adata.obs[k]))}
    if len(custom_colors) > 0:
        write_color_tracks(custom_colors, project_dir + 'color_data_gene_sets.csv')

    all = []

    # save gene colortracks
    os.system('mkdir ' + project_dir + 'gene_colors')
    # The following Split into left right (+ casting) makes sure that every gene
    # is included, no out of bounds
    II = int(len(gene_list) / 50) + 1
    left = 0
    right = II
    for j in range(50):
        fname = project_dir + 'gene_colors/color_data_all_genes-' + repr(j) + '.csv'
        X_writeable_chunk = adata[:, gene_list[left:right]].X
        if issparse(X_writeable_chunk):
            X_writeable_chunk = X_writeable_chunk.toarray()
        if X_writeable_chunk.ndim == 1:
            X_writeable_chunk = X_writeable_chunk[:, None]
        all_gene_colors = {
            g: X_writeable_chunk[:, i]
            for i, g in enumerate(gene_list[left:right])}
        write_color_tracks(all_gene_colors, fname, adata.X.shape[0])
        left += II
        right += II
        if right >= len(gene_list): right = len(gene_list)
        if right <= left: break
        all += all_gene_colors.keys()

    # Create and save a dictionary of color profiles to be used by the visualizer
    # Cast: numpy datatypes as input not json serializable
    # Pre-calculate statistics before writing to speed up calculations
    X = adata.X
    color_stats = {g: (float(np.mean(X[:, i])),
                       float(np.std(X[:, i].todense() if issparse(X) else X[:, i])),
                       float(np.max(X[:, i])),
                       float(np.percentile(
                           X[:, i].todense() if issparse(X) else X[:, i],
                           99)))
                   for i, g in enumerate(gene_list)}
    json.dump(color_stats,
              open(project_dir + '/color_stats.json', 'w'), indent=4, sort_keys=True)

    # save cell labels

    # Categorical coloring data:
    categorical_coloring_data = {}
    # Adapt groupby
    if isinstance(cell_groupings, str):
        cell_groupings = [cell_groupings]

    if cell_groupings is None:
        cell_groupings = [k for k in adata.obs_keys() if
                          is_categorical(adata.obs[k])]
    for j, i in enumerate(cell_groupings):
        if cell_groupings[j] not in adata.obs:
            logg.warn('adata annotation key for cell grouping does not exist. '
                      'Inspect observation annotation.')
        else:
            group_names = []
            groups = adata.obs[cell_groupings[j]]
            group_names = adata.obs[cell_groupings[j]].cat.categories
            add_colors_for_categorical_sample_annotation(adata,
                                                         cell_groupings[j],
                                                         palette=None)
            group_colors = adata.uns[i + '_colors']
            label_colors = {l: group_colors[i] for i, l in enumerate(group_names)}
            labels = list(groups)
            # SPRING expects a Dictionary for label_colors, but a list for labels !
            categorical_coloring_data[cell_groupings[j]] = {
                'label_colors': label_colors, 'label_list': labels}
    json.dump(categorical_coloring_data, open(
              project_dir + '/categorical_coloring_data.json', 'w'), indent=4)

    nodes = [{'name': i, 'number': i} for i in range(X.shape[0])]
    edges = [{'source': int(i), 'target': int(j)} for i, j in edges]
    out = {'nodes': nodes, 'links': edges}
    open(project_dir + 'graph_data.json', 'w').write(
        json.dumps(out, indent=4, separators=(',', ': ')))


# The following method is only used when a full (non-sparse) distance matrix is given as an input parameter
# Depending on input size, this can be very cost-inefficient
def get_knn_edges(dmat, k):
    edge_dict = {}
    for i in range(dmat.shape[0]):
        # Save modified coordinate values, rewrite so that adata_object is not changed!
        l=k
        saved_values={}
        while l>0:
            j = dmat[i, :].argmin()
            saved_values[j]=dmat[i,j]
            if i != j:
                ii, jj = tuple(sorted([i, j]))
                edge_dict[(ii, jj)] = dmat[i, j]
            dmat[i, j] = 0
            l=l-1
        # Rewrite safed values:
        for j, val in enumerate(saved_values):
            dmat[i,j]=val

    return edge_dict.keys()


# This is a (preliminary) alternative to get_knn_edges
# We assume that D is a distance matrix containing only non-zero entries for the (k-1)nn
# (as is the value for data graph distance local)
# This is the version for knn as in graph distance local.
def get_knn_edges_sparse(dmat, k):
    edge_dict = {}
    if not issparse(dmat):
        return get_knn_edges(dmat,k)
    else:
        for i in range(dmat.shape[0]):
            l=1
            saved_values={}
            while l<k:
                row = dmat.getrow(i)
                data_index=row.data.argmin()
                j=row.indices[data_index]
                saved_values[j] = dmat[i, j]
                if i != j:
                    ii, jj = tuple(sorted([i, j]))
                    edge_dict[(ii, jj)] = dmat[i, j]
                dmat[i, j] = inf
                l = l + 1
            # Rewrite safed values:
            for j in saved_values:
                dmat[i, j] = saved_values[j]
    return edge_dict.keys()


def write_color_tracks(ctracks, fname, n_cells=0):
    out = []
    for name, score in ctracks.items():
        line = ','.join([name] + [repr(round(x, 1)) for x in score])
        out += [line]
    out = sorted(out, key=lambda x: x.split(',')[0])
    open(fname, 'w').write('\n'.join(out))
