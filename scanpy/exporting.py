"""Exporting to formats for other software.
"""

import numpy as np
import os
import json
from scipy.sparse import issparse
import logging as logg
from pandas.api.types import is_string_dtype, is_categorical
from .plotting.utils import add_colors_for_categorical_sample_annotation


def spring_project(
        adata, project_dir, use_genes=None, cell_groupings=None,
        custom_color_tracks=None):
    """Exports to a SPRING project directory [Weinreb17]_.

    Visualize annotation present in `adata`. By default, export all categorical
    annotation present in `adata.obs`.

    See `SPRING <https://github.com/AllonKleinLab/SPRING>`_ or [Weinreb17]_ for details.

    Parameters
    ----------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix: `adata.uns['neighbors']` needs to be present.
    project_dir : `str`
        Path to SPRING directory.
    use_genes : `'rank_gene_groups'` or `list`, optional (default: `None`)
        Select a subset of genes. If 'rank_gene_groups', this uses the
        annotation written by :func:`~scanpy.api.tl.rank_gene_groups`.
    cell_groupings : `str`, `list` of `str`, optional (default: `None`)
        Instead of importing all categorical annotation, pass a list of keys for
        `adata.obs`.
    custom_color_tracks : `str`, `list` of `str`, optional (default: `None`)
        Specify specific `adata.obs` keys for continuous coloring.

    Examples
    --------
    See this `tutorial <https://github.com/theislab/scanpy_usage/tree/master/171111_SPRING_export>`_.
    """

    gene_list = adata.var_names
    # We allow to include rank_genes_groups output.
    if isinstance(use_genes, str):
        use_genes_list = []
        if use_genes == 'rank_genes_groups':
            for rank in adata.uns['rank_genes_groups']['names']:
                for groups in rank:
                    use_genes_list.append(groups)
            use_genes = list(set(use_genes_list))
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

    if 'neighbors' not in adata.uns:
        raise ValueError('Run `sc.pp.neighbors` first.')

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
    adata_raw = adata.raw if adata.raw is not None else adata
    II = int(len(gene_list) / 50) + 1
    left = 0
    right = II
    for j in range(50):
        fname = project_dir + 'gene_colors/color_data_all_genes-' + repr(j) + '.csv'
        X_writeable_chunk = adata_raw[:, gene_list[left:right]].X
        if issparse(X_writeable_chunk):
            X_writeable_chunk = X_writeable_chunk.toarray()
        if X_writeable_chunk.ndim == 1:
            X_writeable_chunk = X_writeable_chunk[:, None]
        all_gene_colors = {
            g: X_writeable_chunk[:, i]
            for i, g in enumerate(gene_list[left:right])}
        write_color_tracks(all_gene_colors, fname, adata_raw.X.shape[0])
        left += II
        right += II
        if right >= len(gene_list): right = len(gene_list)
        if right <= left: break
        all += all_gene_colors.keys()

    # Create and save a dictionary of color profiles to be used by the visualizer
    # Cast: numpy datatypes as input not json serializable
    # Pre-calculate statistics before writing to speed up calculations
    X = adata_raw.X
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

    # The actual graph
    nodes = [{'name': i, 'number': i} for i in range(X.shape[0])]
    if 'distances' in adata.uns['neighbors']:  # these are sparse matrices
        matrix = adata.uns['neighbors']['distances']
    else:
        matrix = adata.uns['neighbors']['connectivities']
    matrix = matrix.tocoo()
    edges = [{'source': int(i), 'target': int(j)}
             for i, j in zip(matrix.row, matrix.col)]
    out = {'nodes': nodes, 'links': edges}
    open(project_dir + 'graph_data.json', 'w').write(
        json.dumps(out, indent=4, separators=(',', ': ')))


def write_color_tracks(ctracks, fname, n_cells=0):
    out = []
    for name, score in ctracks.items():
        line = ','.join([name] + [repr(round(x, 1)) for x in score])
        out += [line]
    out = sorted(out, key=lambda x: x.split(',')[0])
    open(fname, 'w').write('\n'.join(out))
