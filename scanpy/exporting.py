# Author: F. Alex Wolf (http://falexwolf.de)
#         T. Callies 
"""Exporting to formats for other software.
"""

import numpy as np
import os
import json
import pdb
from math import inf
from .data_structs.data_graph import add_or_update_graph_in_adata
from scipy.sparse import issparse

def save_spring_dir(adata, k, project_directory, D=None,
                    custom_colors={}, cell_groupings=None, use_genes=[]):
    """Builds a SPRING project directory.
    This is based on a preprocessing function by Caleb Weinreb:
    https://github.com/AllonKleinLab/SPRING/
    Parameters
    ----------
    adata : AnnData() object
        Matrix of gene expression. Rows correspond to cells and columns
        correspond to genes.
    k : int
        Number of edges assigned to each node in knn graph
    project_directory : str
        Path to a directory where SPRING readable files will be written. The
        directory does not have to exist before running this function.
    D : np.ndarray (default: None)
        Distance matrix for construction of knn graph. Any distance matrix can
        be used as long as higher values correspond to greater distances.
        If nothing is given, local graph distance is used (if available) or (re-)computed
    cell_groupings : list of str (default: None)
        Optional list of strings containing adata.add key to grouping. Adata.add should contain cell_groupings+'_order'
        and cell_groupings+'colors' as keys with names / colors for groupings and for each cell the corresponding group
        Furthermore. adata.smp[cell_groupings] should return an array of adata.X.shape[0] elements
    custom_colors : dict (default: {})
        Dictionary with one key-value pair for each custom color.  The key is
        the name of the color track and the value is a list of scalar values
        (i.e. color intensities). If there are N cells total (i.e. X.shape[0] ==
        N), then the list of labels should have N entries.
        Currently not used
    use_genes : list, default: []
        Selects certain genes that are written into the coloring files. Default ([]) selects all genes
    """
    # TODO: (LATER:) Make sure that this works for sparse adata objects as well

    X= adata.X
    # gene_list: list, np.ndarry - like
    # An ordered list of gene names with length X.shape[1].
    gene_list=adata.var_names


    # File can be safed anywhere. However, for easy access via SPRING, safe it somewhere in the spring directory
    os.system('mkdir ' + project_directory)
    if not project_directory[-1] == '/': project_directory += '/'

    if D==None:
        if 'data_graph_distance_local' in adata.add:
            D=adata.add['data_graph_distance_local']
            edges= get_knn_edges_sparse(D, k)
        else:
            # if not available and nothing is given, calculate distances and add data_graph_distance_local to dictionary
            add_or_update_graph_in_adata(
                adata,
                n_neighbors=30,
                n_pcs=50,
                n_dcs=15,
                knn=None,
                recompute_pca=False,
                recompute_distances=False,
                recompute_graph=False,
                n_jobs=None)
            # Note that output here will always be sparse
            D = adata.add['data_graph_distance_local']
            edges = get_knn_edges_sparse(D, k)
    else:
        if issparse(D):
            edges = get_knn_edges_sparse(D, k)
        else:
            edges = get_knn_edges(D,k)

    # save genesets
    #TODO: (LATER:) Include when everything else works. Check how to include efficiently
    # custom_colors['Uniform'] = np.zeros(X.shape[0])
    # write_color_tracks(custom_colors, project_directory + 'color_data_gene_sets.csv')
    all = []

    # save gene colortracks
    os.system('mkdir ' + project_directory + 'gene_colors')
    # The following Split into left right (+ casting) makes sure that every gene is included, no out of bounds
    II = int(len(gene_list) / 50) + 1
    left=0
    right=II
    for j in range(50):
        fname = project_directory + '/gene_colors/color_data_all_genes-' + repr(j) + '.csv'
        if len(use_genes) > 0:
            all_gene_colors = {
                # Adapted slicing, so that it won't be OOB
                g: X[:, i + left] for i, g in enumerate(gene_list[left: right]) if g in use_genes}
        else:
            all_gene_colors = {
                # Here, the original control mechansim to included genes in Color if average expression across all cells is above 0.05
                # This doesn't translate to anything meaningful here. Actually, no genes are then selected
                g: X[:, i + left] for i, g in enumerate(
                gene_list[left: right]) }
        write_color_tracks(all_gene_colors, fname)
        left+=II
        right+=II
        # Avoid OOB:
        if right >=len(gene_list):
            right=len(gene_list)
        all += all_gene_colors.keys()

    # Create and save a dictionary of color profiles to be used by the visualizer
    # Cast: numpy datatypes as input not json serializable
    color_stats = {}
    for i in range(X.shape[1]):
        mean = float(np.mean(X[:, i]))
        std = float(np.std(X[:, i]))
        max = float(np.max(X[:, i]))
        centile = float(np.percentile(X[:, i], 99.6))
        color_stats[gene_list[i]] = (mean, std, 0, max, centile)
    json.dump(color_stats,
              open(project_directory + '/color_stats.json', 'w'), indent=4, sort_keys=True)

    # save cell labels
    # Categorical coloring data:
    categorical_coloring_data = {}
    # Adapt groupby
    if cell_groupings is None:
        # In this case, do nothing
        pass
    else:
        for j in range(cell_groupings):
            if (cell_groupings[j]+'_order' not in adata.add) or (cell_groupings[j]+'_colors' not in adata.add) :
                # TODO: Change to logging
                print('Adata annotation does not exist. Check input' )
            else:
                groups=adata.smp[cell_groupings[j]]
                group_names=adata.add[cell_groupings[j]+'_order']
                group_colors=adata.add[cell_groupings[j]+'_colors']
                label_colors = {l: group_colors[i] for i, l in enumerate(group_names)}
                labels = list(groups)
                # SPRING expects a Dictionary for label_colors, but a list for labels !
                categorical_coloring_data[cell_groupings[j]] = {'label_colors': label_colors, 'label_list': labels}
        json.dump(categorical_coloring_data, open(
                project_directory + '/categorical_coloring_data.json', 'w'), indent=4)
    #

    nodes = [{'name': i, 'number': i} for i in range(X.shape[0])]
    edges = [{'source': int(i), 'target': int(j)} for i, j in edges]
    out = {'nodes': nodes, 'links': edges}
    # Possible Error: ' instead of ": For now, it seems to work
    open(project_directory + 'graph_data.json', 'w').write(
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
            row = dmat.getrow(i)
            l=1
            saved_values={}
            while l<k:
                data_index=row.data.argmin()
                j=row.indices[data_index]
                saved_values[j] = dmat[i, j]
                if i != j:
                    ii, jj = tuple(sorted([i, j]))
                    edge_dict[(ii, jj)] = dmat[i, j]
                dmat[i, j] = inf
                l = l + 1
            # Rewrite safed values:
        for j, val in enumerate(saved_values):
            dmat[i, j] = val
    return edge_dict.keys()
    
def write_color_tracks(ctracks, fname):
    out = []
    for name, score in ctracks.items():
        line = ','.join([name] + [repr(round(x, 1)) for x in score])
        out += [line]
    out = sorted(out, key=lambda x: x.split(',')[0])
    open(fname, 'w').write('\n'.join(out))
