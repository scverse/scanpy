# Author: F. Alex Wolf (http://falexwolf.de)
"""Exporting to formats for other software.
"""

def save_spring_dir(X, D, k, gene_list, project_directory,
                    custom_colors={}, cell_groupings={}, use_genes=[]):
    """Builds a SPRING project directory.

    This is based on a preprocessing function by Caleb Weinreb:
    https://github.com/AllonKleinLab/SPRING/

    Parameters
    ----------
    X : np.ndarray
        Matrix of gene expression. Rows correspond to cells and columns
        correspond to genes.
    D : np.ndarray
        Distance matrix for construction of knn graph. Any distance matrix can
        be used as long as higher values correspond to greater distances.
    k : int
        Number of edges assigned to each node in knn graph
    gene_list : list, np.ndarry-like
        An ordered list of gene names with length X.shape[1].
    project_directory : str
        Path to a directory where SPRING readable files will be written. The
        directory does not have to exist before running this function.
    cell_groupings : dict
        Dictionary with one key-value pair for each cell grouping.  The key is
        the name of the grouping (e.g. "SampleID") and the value is a list of
        labels (e.g. ["sample1","sample2"...])  If there are N cells total
        (i.e. X.shape[0] == N), then the list of labels should have N entries.
    custom_colors : dict
        Dictionary with one key-value pair for each custom color.  The key is
        the name of the color track and the value is a list of scalar values
        (i.e. color intensities). If there are N cells total (i.e. X.shape[0] ==
        N), then the list of labels should have N entries.
    """
    os.system('mkdir ' + project_directory)
    if not project_directory[-1] == '/': project_directory += '/'
    # Build graph
    edges = get_knn_edges(D, k)

    # save genesets
    custom_colors['Uniform'] = np.zeros(X.shape[0])
    write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')
    all = []

    # save gene colortracks
    os.system('mkdir '+project_directory+'gene_colors')
    II = len(gene_list) / 50 + 1
    for j in range(50):
        fname = project_directory+'/gene_colors/color_data_all_genes-' + repr(j) + '.csv'
        if len(use_genes) > 0: all_gene_colors = {
            g: X[:, i+II*j] for i, g in enumerate(gene_list[II*j: II*(j+1)]) if g in use_genes}
        else:
            all_gene_colors = {
                g: X[:, i+II*j] for i, g in enumerate(
                    gene_list[II*j: II*(j+1)]) if np.mean(X[:, i+II*j]) > 0.05}
        write_color_tracks(all_gene_colors, fname)
        all += all_gene_colors.keys()

    # Create and save a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    for i in range(X.shape[1]):
        mean = np.mean(X[:, i])
        std = np.std(X[:, i])
        max = np.max(X[:, i])
        centile = np.percentile(X[:, i], 99.6)
        color_stats[gene_list[i]] = (mean, std, 0, max, centile)
    for k, v in custom_colors.items():
        color_stats[k] = (0, 1, np.min(v), np.max(v)+.01, np.percentile(v, 99))
    json.dump(color_stats,
              open(project_directory + '/color_stats.json', 'w'), indent=4, sort_keys=True)

    # save cell labels
    categorical_coloring_data = {}
    from matplotlib.colors import cnames
    for k, labels in cell_groupings.items():
        label_colors = {l: frac_to_hex(float(i)/len(set(labels)))
                                       for i, l in enumerate(list(set(labels)))}
        if k == 'ctpaths':
            label_colors['dontknow'] = cnames['grey']
        if k == 'celltypes':
            label_colors['no_gate'] = cnames['grey']
        categorical_coloring_data[k] = {'label_colors': label_colors, 'label_list': labels}
    json.dump(categorical_coloring_data, open(
        project_directory + '/categorical_coloring_data.json', 'w'), indent=4)

    nodes = [{'name': i, 'number': i} for i in range(X.shape[0])]
    edges = [{'source': i, 'target': j, 'distance': 0} for i, j in edges]
    out = {'nodes': nodes, 'links': edges}
    open(project_directory + 'graph_data.json', 'w').write(
        json.dumps(out, indent=4, separators=(',', ': ')))

    
def get_knn_edges(dmat, k):
    edge_dict = {}
    for i in range(dmat.shape[0]):
        for j in np.nonzero(dmat[i,:] <= sorted(dmat[i, :])[k])[0]:
            if i != j:
                ii,jj = tuple(sorted([i, j]))
                edge_dict[(ii, jj)] = dmat[i,j]
    return edge_dict.keys()
    

def write_color_tracks(ctracks, fname):
    out = []
    for name, score in ctracks.items():
        line = ','.join([name]+[repr(round(x, 1)) for x in score])
        out += [line]
    out = sorted(out, key=lambda x: x.split(',')[0])
    open(fname, 'w').write('\n'.join(out))
