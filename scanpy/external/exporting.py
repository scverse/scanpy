"""\
Exporting to formats for other software.
"""
import json
import logging as logg
from pathlib import Path
from typing import Union, Optional, Iterable, Mapping

import numpy as np
import scipy.sparse
import h5py
import matplotlib.pyplot as plt
from anndata import AnnData
from pandas.api.types import is_categorical_dtype

from ..preprocessing._utils import _get_mean_var
from .._utils import NeighborsView


def spring_project(
    adata: AnnData,
    project_dir: Union[Path, str],
    embedding_method: str,
    subplot_name: Optional[str] = None,
    cell_groupings: Union[str, Iterable[str], None] = None,
    custom_color_tracks: Union[str, Iterable[str], None] = None,
    total_counts_key: str = 'n_counts',
    neighbors_key: Optional[str] = None,
    overwrite: bool = False,
):
    """\
    Exports to a SPRING project directory [Weinreb17]_.

    Visualize annotation present in `adata`. By default, export all gene expression data
    from `adata.raw` and categorical and continuous annotations present in `adata.obs`.

    See `SPRING <https://github.com/AllonKleinLab/SPRING>`__ or [Weinreb17]_ for details.

    Parameters
    ----------
    adata
        Annotated data matrix: `adata.uns['neighbors']` needs to
        be present.
    project_dir
        Path to directory for exported SPRING files.
    embedding_method
        Name of a 2-D embedding in `adata.obsm`
    subplot_name
        Name of subplot folder to be created at `project_dir+"/"+subplot_name`
    cell_groupings
        Instead of importing all categorical annotations when `None`,
        pass a list of keys for `adata.obs`.
    custom_color_tracks
        Specify specific `adata.obs` keys for continuous coloring.
    total_counts_key
        Name of key for total transcript counts in `adata.obs`.
    overwrite
        When `True`, existing counts matrices in `project_dir` are overwritten.

    Examples
    --------
    See this `tutorial <https://github.com/scverse/scanpy_usage/tree/master/171111_SPRING_export>`__.
    """

    # need to get nearest neighbors first
    if neighbors_key is None:
        neighbors_key = 'neighbors'

    if neighbors_key not in adata.uns:
        raise ValueError('Run `sc.pp.neighbors` first.')

    # check that requested 2-D embedding has been generated
    if embedding_method not in adata.obsm_keys():
        if 'X_' + embedding_method in adata.obsm_keys():
            embedding_method = 'X_' + embedding_method
        else:
            if embedding_method in adata.uns:
                embedding_method = (
                    'X_'
                    + embedding_method
                    + '_'
                    + adata.uns[embedding_method]['params']['layout']
                )
            else:
                raise ValueError(
                    'Run the specified embedding method `%s` first.' % embedding_method
                )

    coords = adata.obsm[embedding_method]

    # Make project directory and subplot directory (subplot has same name as project)
    # For now, the subplot is just all cells in adata
    project_dir: Path = Path(project_dir)
    subplot_dir: Path = (
        project_dir.parent if subplot_name is None else project_dir / subplot_name
    )
    subplot_dir.mkdir(parents=True, exist_ok=True)
    print(f'Writing subplot to {subplot_dir}')

    # Write counts matrices as hdf5 files and npz if they do not already exist
    # or if user requires overwrite.
    # To do: check if Alex's h5sparse format will allow fast loading from just
    # one file.
    write_counts_matrices = True
    base_dir_filelist = [
        'counts_norm_sparse_genes.hdf5',
        'counts_norm_sparse_cells.hdf5',
        'counts_norm.npz',
        'total_counts.txt',
        'genes.txt',
    ]
    if all((project_dir / f).is_file() for f in base_dir_filelist):
        if not overwrite:
            logg.warning(
                f'{project_dir} is an existing SPRING folder. A new subplot will be created, but '
                'you must set `overwrite=True` to overwrite counts matrices.'
            )
            write_counts_matrices = False
        else:
            logg.warning(f'Overwriting the files in {project_dir}.')

    # Ideally, all genes will be written from adata.raw
    if adata.raw is not None:
        E = adata.raw.X.tocsc()
        gene_list = list(adata.raw.var_names)
    else:
        E = adata.X.tocsc()
        gene_list = list(adata.var_names)

    # Keep track of total counts per cell if present
    if total_counts_key in adata.obs:
        total_counts = np.array(adata.obs[total_counts_key])
    else:
        total_counts = E.sum(1).A1

    # Write the counts matrices to project directory
    if write_counts_matrices:
        write_hdf5_genes(E, gene_list, project_dir / 'counts_norm_sparse_genes.hdf5')
        write_hdf5_cells(E, project_dir / 'counts_norm_sparse_cells.hdf5')
        write_sparse_npz(E, project_dir / 'counts_norm.npz')
        with (project_dir / 'genes.txt').open('w') as o:
            for g in gene_list:
                o.write(g + '\n')
        np.savetxt(project_dir / 'total_counts.txt', total_counts)

    # Get categorical and continuous metadata
    categorical_extras = {}
    continuous_extras = {}
    if cell_groupings is None:
        for obs_name in adata.obs:
            if is_categorical_dtype(adata.obs[obs_name]):
                categorical_extras[obs_name] = [str(x) for x in adata.obs[obs_name]]
    else:
        if isinstance(cell_groupings, str):
            cell_groupings = [cell_groupings]
        for obs_name in cell_groupings:
            if obs_name not in adata.obs:
                logg.warning(f'Cell grouping {obs_name!r} is not in adata.obs')
            elif is_categorical_dtype(adata.obs[obs_name]):
                categorical_extras[obs_name] = [str(x) for x in adata.obs[obs_name]]
            else:
                logg.warning(
                    f'Cell grouping {obs_name!r} is not a categorical variable'
                )
    if custom_color_tracks is None:
        for obs_name in adata.obs:
            if not is_categorical_dtype(adata.obs[obs_name]):
                continuous_extras[obs_name] = np.array(adata.obs[obs_name])
    else:
        if isinstance(custom_color_tracks, str):
            custom_color_tracks = [custom_color_tracks]
        for obs_name in custom_color_tracks:
            if obs_name not in adata.obs:
                logg.warning(f'Custom color track {obs_name!r} is not in adata.obs')
            elif not is_categorical_dtype(adata.obs[obs_name]):
                continuous_extras[obs_name] = np.array(adata.obs[obs_name])
            else:
                logg.warning(
                    f'Custom color track {obs_name!r} is not a continuous variable'
                )

    # Write continuous colors
    continuous_extras['Uniform'] = np.zeros(E.shape[0])
    _write_color_tracks(continuous_extras, subplot_dir / 'color_data_gene_sets.csv')

    # Create and write a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    color_stats = _get_color_stats_genes(color_stats, E, gene_list)
    color_stats = _get_color_stats_custom(color_stats, continuous_extras)
    _write_color_stats(subplot_dir / 'color_stats.json', color_stats)

    # Write categorical data
    categorical_coloring_data = {}
    categorical_coloring_data = _build_categ_colors(
        categorical_coloring_data, categorical_extras
    )
    _write_cell_groupings(
        subplot_dir / 'categorical_coloring_data.json', categorical_coloring_data
    )

    # Write graph in two formats for backwards compatibility
    edges = _get_edges(adata, neighbors_key)
    _write_graph(subplot_dir / 'graph_data.json', E.shape[0], edges)
    _write_edges(subplot_dir / 'edges.csv', edges)

    # Write cell filter; for now, subplots must be generated from within SPRING,
    # so cell filter includes all cells.
    np.savetxt(subplot_dir / 'cell_filter.txt', np.arange(E.shape[0]), fmt='%i')
    np.save(subplot_dir / 'cell_filter.npy', np.arange(E.shape[0]))

    # Write 2-D coordinates, after adjusting to roughly match SPRING's default d3js force layout parameters
    coords = coords - coords.min(0)[None, :]
    coords = (
        coords * (np.array([1000, 1000]) / coords.ptp(0))[None, :]
        + np.array([200, -200])[None, :]
    )
    np.savetxt(
        subplot_dir / 'coordinates.txt',
        np.hstack((np.arange(E.shape[0])[:, None], coords)),
        fmt='%i,%.6f,%.6f',
    )

    # Write some useful intermediates, if they exist
    if 'X_pca' in adata.obsm_keys():
        np.savez_compressed(
            subplot_dir / 'intermediates.npz',
            Epca=adata.obsm['X_pca'],
            total_counts=total_counts,
        )

    # Write PAGA data, if present
    if 'paga' in adata.uns:
        clusts = np.array(adata.obs[adata.uns['paga']['groups']].cat.codes)
        uniq_clusts = adata.obs[adata.uns['paga']['groups']].cat.categories
        paga_coords = [coords[clusts == i, :].mean(0) for i in range(len(uniq_clusts))]
        _export_PAGA_to_SPRING(adata, paga_coords, subplot_dir / 'PAGA_data.json')


# --------------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------------


def _get_edges(adata, neighbors_key=None):
    neighbors = NeighborsView(adata, neighbors_key)
    if 'distances' in neighbors:  # these are sparse matrices
        matrix = neighbors['distances']
    else:
        matrix = neighbors['connectivities']
    matrix = matrix.tocoo()
    edges = [(i, j) for i, j in zip(matrix.row, matrix.col)]

    return edges


def write_hdf5_genes(E, gene_list, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_genes.hdf5"'''

    E = E.tocsc()

    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    cix_group = hf.create_group('cell_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iG, g in enumerate(gene_list):
        counts = E[:, iG].A.squeeze()
        cell_ix = np.nonzero(counts)[0]
        counts = counts[cell_ix]
        counts_group.create_dataset(g, data=counts)
        cix_group.create_dataset(g, data=cell_ix)

    hf.close()


def write_hdf5_cells(E, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_cells.hdf5"'''

    E = E.tocsr()

    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    gix_group = hf.create_group('gene_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iC in range(E.shape[0]):
        counts = E[iC, :].A.squeeze()
        gene_ix = np.nonzero(counts)[0]
        counts = counts[gene_ix]
        counts_group.create_dataset(str(iC), data=counts)
        gix_group.create_dataset(str(iC), data=gene_ix)

    hf.close()


def write_sparse_npz(E, filename, compressed=False):
    '''SPRING standard: filename = main_spring_dir + "/counts_norm.npz"'''
    E = E.tocsc()
    scipy.sparse.save_npz(filename, E, compressed=compressed)


def _write_graph(filename, n_nodes, edges):
    nodes = [{'name': int(i), 'number': int(i)} for i in range(n_nodes)]
    edges = [{'source': int(i), 'target': int(j), 'distance': 0} for i, j in edges]
    out = {'nodes': nodes, 'links': edges}
    open(filename, 'w').write(json.dumps(out, indent=4, separators=(',', ': ')))


def _write_edges(filename, edges):
    with open(filename, 'w') as f:
        for e in edges:
            f.write('%i;%i\n' % (e[0], e[1]))


def _write_color_tracks(ctracks, fname):
    out = []
    for name, score in ctracks.items():
        line = name + ',' + ','.join(['%.3f' % x for x in score])
        out += [line]
    out = sorted(out, key=lambda x: x.split(',')[0])
    open(fname, 'w').write('\n'.join(out))


def _frac_to_hex(frac):
    rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3]) * 255, dtype=int))
    return '#%02x%02x%02x' % rgb


def _get_color_stats_genes(color_stats, E, gene_list):
    means, variances = _get_mean_var(E)
    stdevs = np.zeros(variances.shape, dtype=float)
    stdevs[variances > 0] = np.sqrt(variances[variances > 0])
    mins = E.min(0).todense().A1
    maxes = E.max(0).todense().A1

    pctl = 99.6
    pctl_n = (100 - pctl) / 100.0 * E.shape[0]
    pctls = np.zeros(E.shape[1], dtype=float)
    for iG in range(E.shape[1]):
        n_nonzero = E.indptr[iG + 1] - E.indptr[iG]
        if n_nonzero > pctl_n:
            pctls[iG] = np.percentile(
                E.data[E.indptr[iG] : E.indptr[iG + 1]], 100 - 100 * pctl_n / n_nonzero
            )
        else:
            pctls[iG] = 0
        color_stats[gene_list[iG]] = tuple(
            map(float, (means[iG], stdevs[iG], mins[iG], maxes[iG], pctls[iG]))
        )
    return color_stats


def _get_color_stats_custom(color_stats, custom_colors):
    for k, v in custom_colors.items():
        color_stats[k] = tuple(
            map(
                float,
                (np.mean(v), np.std(v), np.min(v), np.max(v), np.percentile(v, 99)),
            )
        )
    return color_stats


def _write_color_stats(filename, color_stats):
    with open(filename, 'w') as f:
        f.write(json.dumps(color_stats, indent=4, sort_keys=True))  # .decode('utf-8'))


def _build_categ_colors(categorical_coloring_data, cell_groupings):
    for k, labels in cell_groupings.items():
        label_colors = {
            l: _frac_to_hex(float(i) / len(set(labels)))
            for i, l in enumerate(list(set(labels)))
        }
        categorical_coloring_data[k] = {
            'label_colors': label_colors,
            'label_list': labels,
        }
    return categorical_coloring_data


def _write_cell_groupings(filename, categorical_coloring_data):
    with open(filename, 'w') as f:
        f.write(
            json.dumps(categorical_coloring_data, indent=4, sort_keys=True)
        )  # .decode('utf-8'))


def _export_PAGA_to_SPRING(adata, paga_coords, outpath):
    # retrieve node data
    group_key = adata.uns['paga']['groups']
    names = adata.obs[group_key].cat.categories
    coords = [list(xy) for xy in paga_coords]

    sizes = list(adata.uns[group_key + '_sizes'])
    clus_labels = adata.obs[group_key].cat.codes.values
    cell_groups = [
        [int(j) for j in np.nonzero(clus_labels == i)[0]] for i in range(len(names))
    ]

    if group_key + '_colors' in adata.uns:
        colors = list(adata.uns[group_key + '_colors'])
    else:
        import scanpy.plotting.utils

        scanpy.plotting.utils.add_colors_for_categorical_sample_annotation(
            adata, group_key
        )
        colors = list(adata.uns[group_key + '_colors'])

    # retrieve edge level data
    sources, targets = adata.uns['paga']['connectivities'].nonzero()
    weights = np.sqrt(adata.uns['paga']['connectivities'].data) / 3

    # save a threshold weight for showing edges so that by default,
    # the number of edges shown is 8X the number of nodes
    if len(names) * 8 > len(weights):
        min_edge_weight_view = 0
    else:
        min_edge_weight_view = sorted(weights)[-len(names) * 8]

    # save another threshold for even saving edges at all, with 100 edges per node
    if len(weights) < 100 * len(names):
        min_edge_weight_save = 0
    else:
        min_edge_weight_save = sorted(weights)[-len(names) * 100]

    # make node list
    nodes = []
    for i, name, xy, color, size, cells in zip(
        range(len(names)), names, coords, colors, sizes, cell_groups
    ):
        nodes.append(
            {
                'index': i,
                'size': int(size),
                'color': color,
                'coordinates': xy,
                'cells': cells,
                'name': name,
            }
        )

    # make link list, avoid redundant encoding (graph is undirected)
    links = []
    for source, target, weight in zip(sources, targets, weights):
        if source < target and weight > min_edge_weight_save:
            links.append(
                {'source': int(source), 'target': int(target), 'weight': float(weight)}
            )

    # save data about edge weights
    edge_weight_meta = {
        'min_edge_weight': min_edge_weight_view,
        'max_edge_weight': np.max(weights),
    }

    PAGA_data = {'nodes': nodes, 'links': links, 'edge_weight_meta': edge_weight_meta}

    import json

    json.dump(PAGA_data, open(outpath, 'w'), indent=4)

    return None


def cellbrowser(
    adata: AnnData,
    data_dir: Union[Path, str],
    data_name: str,
    embedding_keys: Union[Iterable[str], Mapping[str, str], str, None] = None,
    annot_keys: Union[Iterable[str], Mapping[str, str], None] = (
        "louvain",
        "percent_mito",
        "n_genes",
        "n_counts",
    ),
    cluster_field: str = "louvain",
    nb_marker: int = 50,
    skip_matrix: bool = False,
    html_dir: Union[Path, str, None] = None,
    port: Optional[int] = None,
    do_debug: bool = False,
):
    """\
    Export adata to a UCSC Cell Browser project directory. If `html_dir` is
    set, subsequently build the html files from the project directory into
    `html_dir`. If `port` is set, start an HTTP server in the background and
    serve `html_dir` on `port`.

    By default, export all gene expression data from `adata.raw`, the
    annotations `louvain`, `percent_mito`, `n_genes` and `n_counts` and the top
    `nb_marker` cluster markers. All existing files in data_dir are
    overwritten, except `cellbrowser.conf`.

    See `UCSC Cellbrowser <https://github.com/maximilianh/cellBrowser>`__ for
    details.

    Parameters
    ----------
    adata
        Annotated data matrix
    data_dir
        Path to directory for exported Cell Browser files.
        Usually these are the files `exprMatrix.tsv.gz`, `meta.tsv`,
        coordinate files like `tsne.coords.tsv`,
        and cluster marker gene lists like `markers.tsv`.
        A file `cellbrowser.conf` is also created with pointers to these files.
        As a result, each adata object should have its own project_dir.
    data_name
        Name of dataset in Cell Browser, a string without special characters.
        This is written to `data_dir/cellbrowser.conf`.
        Ideally this is a short unique name for the dataset,
        like `"pbmc3k"` or `"tabulamuris"`.
    embedding_keys
        2-D embeddings in `adata.obsm` to export.
        The prefix `X_` or `X_draw_graph_` is not necessary.
        Coordinates missing from `adata` are skipped.
        By default (or when specifying `'all'` or `None`), these keys are tried:
        [`"tsne"`, `"umap"`, `"pagaFa"`, `"pagaFr"`, `"pagaUmap"`, `"phate"`,
        `"fa"`, `"fr"`, `"kk"`, `"drl"`, `"rt"`, `"trimap"`].
        For these, default display labels are automatically used.
        For other values, you can specify a mapping from coordinate name to
        display label, e.g. `{"tsne": "t-SNE by Scanpy"}`.
    annot_keys
        Annotations in `adata.obsm` to export.
        Can be a mapping from annotation column name to display label.
        Specify `None` for all available columns in `.obs`.
    skip_matrix
        Do not export the matrix.
        If you had previously exported this adata into the same `data_dir`,
        then there is no need to export the whole matrix again.
        This option will make the export a lot faster,
        e.g. when only coordinates or meta data were changed.
    html_dir
        If this variable is set, the export will build html
        files from `data_dir` to `html_dir`, creating html/js/json files.
        Usually there is one global html output directory for all datasets.
        Often, `html_dir` is located under a webserver's (like Apache)
        htdocs directory or is copied to one.
        A directory `html_dir`/`project_name` will be created and
        an index.html will be created under `html_dir` for all subdirectories.
        Existing files will be overwritten.
        If do not to use html_dir,
        you can use the command line tool `cbBuild` to build the html directory.
    port
        If this variable and `html_dir` are set,
        Python's built-in web server will be spawned as a daemon in the
        background and serve the files under `html_dir`.
        To kill the process, call `cellbrowser.cellbrowser.stop()`.
    do_debug
        Activate debugging output

    Examples
    --------
    See this
    `tutorial <https://github.com/scverse/scanpy_usage/tree/master/181126_Cellbrowser_exports>`__.
    """

    try:
        import cellbrowser.cellbrowser as cb
    except ImportError:
        logg.error(
            "The package cellbrowser is not installed. "
            "Install with 'pip install cellbrowser' and retry."
        )
        raise

    data_dir = str(data_dir)

    cb.setDebug(do_debug)
    cb.scanpyToCellbrowser(
        adata,
        data_dir,
        data_name,
        coordFields=embedding_keys,
        metaFields=annot_keys,
        clusterField=cluster_field,
        nb_marker=nb_marker,
        skipMatrix=skip_matrix,
        doDebug=None,
    )

    if html_dir is not None:
        html_dir = str(html_dir)
        cb.build(data_dir, html_dir, doDebug=None)
        if port is not None:
            cb.serve(html_dir, port)
