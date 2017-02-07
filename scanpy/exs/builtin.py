"""
Example Data and Use Cases - Builtin Examples

Summarizes example data and use cases.

Provides a set of functions for reading raw data, annotating the raw data and
preprocessing of the raw data.

Attributes
----------
Functions read data and do preprocessing, furthermore

dexdata: dict
    Stores information about example data.
dexamples: dict 
    Stores information about example use cases. The keys in this dictionary also
    exist as a function attribute of this module that performs data reading and
    preprocessing.
"""

# this is necessary to import scanpy from within package
from __future__ import absolute_import, print_function
# standard modules
from collections import OrderedDict as odict
# scientific modules
import numpy as np
# scanpy
import scanpy as sc
from .. import utils
from .. import settings as sett

#--------------------------------------------------------------------------------
# The 'dexdata dictionary' stores information about example data.
# - please respect formatting of the 'addedby' entry as 
#   "Initials Surname (github_name), 2016-12-15"
#--------------------------------------------------------------------------------

dexdata = {
'burczynski06': {
    'ref': 'Burczynski et al., J Mol Diagn 8, 51 (2006)',
    'title': 'Molecular classification of Crohn\'s disease and ulcerative colitis '
             'patients using transcriptional profiles in peripheral blood '
             'mononuclear cells',
    'doi': '10.2353/jmoldx.2006.050079',
    'type': 'bulk',
    'addedby': 'FA Wolf (falexwolf), 2016-12-15' },
'krumsiek11': {
    'ref': 'Krumsiek et al., PLoS ONE 6, e22649 (2011)',
    'title': 'Hierarchical Differentiation of Myeloid Progenitors Is Encoded in '
             'the Transcription Factor Network',
    'doi': '10.1371/journal.pone.0022649',
    'type': 'simulated',
    'addedby': 'FA Wolf (falexwolf), 2016-12-15' },
'moignard15': {
    'ref': 'Moignard et al., Nature Biotechnology 33, 269 (2015)',
    'title': 'Decoding the regulatory network of early blood development from '
             'single-cell gene expression measurements',
    'type': 'scqPCR',
    'doi': '10.1038/nbt.3154',
    'addedby': 'FA Wolf (falexwolf), 2016-12-15' },
'paul15': {
    'ref': 'Paul et al., Cell 163, 1663 (2015)',
    'title': 'Transcriptional Heterogeneity and Lineage Commitment in Myeloid '
             'Progenitors',
    'type': 'scRNAseq',
    'doi': '10.1016/j.cell.2015.11.013',
    'addedby': 'FA Wolf (falexwolf), 2016-12-15' },
'toggleswitch': {
    'title': 'Simple toggle switch model.',
    'type': 'simulated',
    'addedby': 'FA Wolf (falexwolf), 2016-12-15' },
} 

#--------------------------------------------------------------------------------
# The 'example dictionary' provides keys that are used to select a preprocessing
# method with the same name. Also, it provides information about tool parameters 
# that deviate from default parameters.
#
# By default, any 'examplekey' ('exkey') is used as 'datakey'. If the
# 'examplekey' does not exist as a datakey (e.g. "paul15_alternative"), a
# datakey has to be specified. It will be used to complete dexamples with
# entries from dexdata during runtime.
#
# If you specified an example data file in dexdata above, you can also be lazy
# and omit to specify the example here. An entry will be generated
# automatically from the dexdata, assuming default settings everywhere.
#--------------------------------------------------------------------------------

dexamples = {
'krumsiek11': {
    'dpt': { 
        'nr_branchings': 2, # detect two branching points (default 1)
        'allow_branching_at_root': True }, # allow branching directly at root
    'ctpaths': { 
        'k': 5,
        'num_fates': 4, # detect two branching points (default 1)
        }
    },
'moignard15': {
    'ctpaths': { 
        'k': 5,
        'num_fates': 2, # detect two branching points (default 1)
        }
    },
'paul15': {
    'paths': { 
        'fates': {0: 877, 1: 2156},
        'k': 20, # increase number of neighbors (default 5)
        'knn': True }, # set a hard threshold on number of neighbors
    'dpt/diffmap': {'k': 20, 'knn': True},
    'difftest': {'log': False, 'groups': 'GMP,MEP'},
    'tgdyn': {'groups': 'GMP,MEP'}
    },
'paul15pca': {
    'datakey': 'paul15',
    'paths': { 
        'fates': {0: 193, 1: 2201},
        'num_fates': 2,
        'k': 4, # increase number of neighbors (default 5)
        'knn': True}, # set a hard threshold on number of neighbors
    'dpt/diffmap': {'k': 20, 'knn': True},
    'difftest': {'log': False, 'groups': 'GMP,MEP'},
    'tgdyn': {'groups': 'GMP,MEP'}
    },
'toggleswitch': {
    'ctpaths': {'fates': {0: 95, 1: 189}},
    'difftest': {'log': False}
    }
}

#--------------------------------------------------------------------------------
# One function per example that reads, annotates and preprocesses data
# - one function 'exkey()' per 'exkey' (key in dexamples)
#--------------------------------------------------------------------------------

def burczynski06():
    """
    Bulk data with conditions ulcerative colitis (UC) and Crohn's disease (CD).

    The study assesses transcriptional profiles in peripheral blood mononuclear
    cells from 42 healthy individuals, 59 CD patients, and 26 UC patients by
    hybridization to microarrays interrogating more than 22,000 sequences.

    Reference
    ---------
    Burczynski ME, Peterson RL, Twine NC, Zuberek KA et al. 
    "Molecular classification of Crohn's disease and ulcerative colitis patients
    using transcriptional profiles in peripheral blood mononuclear cells"
    J Mol Diagn 8, 51 (2006). PMID:16436634.

    https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1615

    Note
    ----
    The function is based on a script by Kerby Shedden.
    http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html
    """
    filename = 'data/burczynski06/GDS1615_full.soft.gz'
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz'
    ddata = sc.read(filename, backup_url=url)
    groupnames_n = ddata['groupnames_n']
    # locations (indices) of samples for the ulcerative colitis group
    locs_UC = [i for i, x in enumerate(groupnames_n) if x == 'ulcerative colitis']
    # locations (indices) of samples for the Crohn's disease group
    locs_CD = [i for i, x in enumerate(groupnames_n) if x == 'Crohn\'s disease']
    grouplocs = [locs_UC, locs_CD]
    # this is just a label that distinguishes the sets
    ddata['grouplabels'] = ['ulcerative colitis','Crohn\'s disease']
    ddata['groupmasks'] = sc.utils.masks(grouplocs, ddata['X'].shape[0])
    # this is not actually needed
    ddata['STP'] = groupnames_n
    ddata['UC'] = locs_UC
    ddata['CD'] = locs_CD
    return ddata

def krumsiek11():
    """
    Simulated myeloid progenitor data.

    Uses a literature-curated boolean network from the reference below.

    Reference
    ---------
    Krumsiek et al.
    "Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the
    Transcription Factor Network"
    PLoS ONE 6, e22649 (2011).

    Returns
    -------
    See paul15().
    """ 
    filename = 'write/krumsiek11_sim/sim_000000.txt'
    ddata = sc.read(filename, first_column_names=True)
    ddata['xroot'] = ddata['X'][0]
    return ddata

def moignard15():
    """
    Hematopoiesis in early mouse embryos.

    Reference
    ---------
    Moignard et al.,
    "Decoding the regulatory network of early blood development from single-cell
    gene expression measurements"
    Nature Biotechnology 33, 269 (2015)

    Returns
    -------
    See paul15. 
    """
    ddata = moignard15_raw() 
    return ddata

def paul15():
    """ 
    Get preprocessed data matrix, gene names, cell names, and root cell.

    This largely follows an R tutorial by Maren Buttner.
    https://github.com/theislab/scAnalysisTutorial

    Reference
    ---------
    Paul et al.,
    "Transcriptional Heterogeneity and Lineage Commitment in Myeloid
    Progenitors",
    Cell 163, 1663 (2015)

    Returns
    -------
    ddata: dict containing
        X: np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        rownames: np.ndarray
            Array storing the experimental labels of samples.
        colnames: np.ndarray
            Array storing the names of genes.
        xroot: np.ndarray
            Expression vector of root cell.
    """
    ddata = paul15_raw()
    ddata['X'] = sc.pp('log', ddata['X'])
    # adjust expression vector of root cell
    ddata['xroot'] = ddata['X'][ddata['iroot']]
    return ddata

def paul15pca():
    ddata = paul15_raw()
    ddata['X'] = sc.pp('log', ddata['X'])
    #ddata['X'] = sc.pp.log(ddata['X'])
    # reduce to 50 components
    ddata['Xpca'] = sc.pca(ddata['X'], nr_comps=50)
    # adjust expression vector of root cell
    ddata['xroot'] = ddata['Xpca'][ddata['iroot']]
    return ddata    

def toggleswitch():
    """ 
    Returns
    -------
    See paul15.
    """
    filename = 'write/toggleswitch_sim/sim_000000.txt'
    ddata = sc.read(filename, first_column_names=True)
    ddata['xroot'] = ddata['X'][0]
    return ddata

#--------------------------------------------------------------------------------
# Optional functions for Raw Data, Annotation, Postprocessing, respectively
# - instead of having just one function per example in the section above, one we
#   split this in parts for 'Raw Data', and, if wanted, tool-specific
#   post-processing (e.g. annotation of groups identified by tools)
# - this is useful, if one wants to experiment with different preprocessing
#   steps, all of which require the same raw data, annotation, and
#   postprocessing steps
#--------------------------------------------------------------------------------

def moignard15_raw():
    """
    1. Filter out a few genes.
    2. Choose 'root cell'.
    3. Define groupnames by inspecting cellnames.
    """
    filename = 'data/moignard15/nbt.3154-S3.xlsx'
    url = 'http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3154-S3.xlsx'
    ddata = sc.read(filename, sheet='dCt_values.txt', backup_url=url)
    X = ddata['X'] # data matrix
    genenames = ddata['colnames'] 
    cellnames = ddata['rownames'] 
    # filter genes
    # filter out the 4th column (Eif2b1), the 31nd (Mrpl19), the 36th
    # (Polr2a) and the 45th (last,UBC), as done by Haghverdi et al. (2016)
    genes = np.r_[np.arange(0, 4), np.arange(5, 31), 
                  np.arange(32, 36), np.arange(37, 45)]
    print('selected', len(genes), 'genes')
    ddata['X'] = X[:, genes] # filter data matrix
    ddata['colnames'] = genenames[genes] # filter genenames
    # choose root cell as in Haghverdi et al. (2016)
    ddata['iroot'] = 532 # note that in Matlab/R, counting starts at 1
    ddata['xroot'] = ddata['X'][ddata['iroot']] 
    # annotate rows of X with Moignard et al. (2015) experimental cell groups
    groups_names = ['HF', 'NP', 'PS', '4SG', '4SFG']
    groups = [] # a list with n entries (one for each sample)
    for name in cellnames:
        for groupname in groups_names:
            if name.startswith(groupname):
                groups.append(groupname)
    # groups are categorical annotation for rows of X
    ddata['rowcat'] = {'groups': groups}
    # add additional metadata for the groups, we want them to be ordered
    # and we want custom colors for each group
    ddata['groups_names'] = groups_names    
    ddata['groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    return ddata

def moignard15_dpt(ddpt):
    # switch on annotation by uncommenting the following
    groups_names = ['trunk', 'undecided/endothelial', 
                  'endothelial', 'erythrocytes']
    ddpt['groups_names'] = [str(i) + ': ' + n for i, n in enumerate(groups_names)]
    return ddpt

def paul15_raw():
    filename = 'data/paul15/paul15.h5'
    url = 'http://falexwolf.de/data/paul15.h5'
    ddata = sc.read(filename, 'data.debatched', backup_url=url)
    # the data has to be transposed (in the hdf5 and R files, each row
    # corresponds to one gene, we use the opposite convention)
    ddata = utils.transpose_ddata(ddata)
    # define local variables to manipulate
    X = ddata['X']
    genenames = ddata['colnames']
    # cluster assocations identified by Paul et al.
    # groups = sc.read(filename,'cluster.id')['X']
    infogenenames = sc.read(filename, 'info.genes_strings')['X']
    # print('the first 10 informative gene names are \n',infogenenames[:10])
    # just keep the first of the equivalent names for each gene
    genenames = np.array([gn.split(';')[0] for gn in genenames])
    # print('the first 10 trunkated gene names are \n',genenames[:10])
    # mask array for the informative genes
    infogenes_idcs = np.array([(True if gn in infogenenames else False)
                                for gn in genenames])
    # restrict data array to the 3451 informative genes
    X = X[:, infogenes_idcs]
    genenames = genenames[infogenes_idcs]
    # print('after selecting info genes, the first 10 gene names are \n',
    #       genenames[:10])
    # write to dict
    ddata['X'] = X
    ddata['colnames'] = genenames
    # set root cell as in Haghverdi et al. (2016)
    ddata['iroot'] = 840 # note that in Matlab/R, counting starts at 1
    ddata['xroot'] = X[ddata['iroot']] 
    return ddata

def paul15_dpt(ddpt):
    ddpt['groups_names'] = ['','GMP','','MEP']
    return ddpt

def paul15pca_dpt(ddpt):
    ddpt['groups_names'] = ['','','GMP','MEP']
    return ddpt

def paul15_paths(dtool):
    dtool['groups_names'] = ['GMP','MEP']
    return dtool

def paul15pca_paths(dtool):
    dtool['groups_names'] = ['GMP','MEP']
    return dtool

