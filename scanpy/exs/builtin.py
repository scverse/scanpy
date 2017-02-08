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
from ..ann_data import AnnData

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

    Available from https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1615.

    Reference
    ---------
    Burczynski ME, Peterson RL, Twine NC, Zuberek KA et al. 
    "Molecular classification of Crohn's disease and ulcerative colitis patients
    using transcriptional profiles in peripheral blood mononuclear cells"
    J Mol Diagn 8, 51 (2006). PMID:16436634.
    """
    filename = 'data/burczynski06/GDS1615_full.soft.gz'
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz'
    ddata = sc.read(filename, backup_url=url)
    adata = AnnData(ddata)
    return adata

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
    """ 
    filename = 'write/krumsiek11_sim/sim_000000.txt'
    ddata = sc.read(filename, first_column_names=True)
    adata = AnnData(ddata)
    adata['xroot'] = adata.X[0]
    return adata

def moignard15():
    """
    Hematopoiesis in early mouse embryos.

    Reference
    ---------
    Moignard et al.,
    "Decoding the regulatory network of early blood development from single-cell
    gene expression measurements"
    Nature Biotechnology 33, 269 (2015)
    """
    adata = moignard15_raw()
    return adata

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
    adata: AnnData
    """
    adata = paul15_raw()
    adata.X = sc.pp.log(adata.X)
    # adjust expression vector of root cell
    adata['xroot'] = adata.X[adata['iroot']]
    return adata

def paul15pca():
    adata = paul15_raw()
    adata.X = sc.pp.log(adata.X)
    # reduce to 50 components
    adata['Xpca'] = sc.pca(adata.X, nr_comps=50)
    # adjust expression vector of root cell
    adata['xroot'] = adata['Xpca'][adata['iroot']]
    return adata    

def toggleswitch():
    """ 
    """
    filename = 'write/toggleswitch_sim/sim_000000.txt'
    ddata = sc.read(filename, first_column_names=True)
    adata = AnnData(ddata)
    adata['xroot'] = adata.X[0]
    return adata

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
    adata = AnnData(ddata)
    # filter genes
    # filter out the 4th column (Eif2b1), the 31nd (Mrpl19), the 36th
    # (Polr2a) and the 45th (last,UBC), as done by Haghverdi et al. (2016)
    genes = np.r_[np.arange(0, 4), np.arange(5, 31), 
                  np.arange(32, 36), np.arange(37, 45)]
    print('selected', len(genes), 'genes')
    adata = adata[:, genes] # filter data matrix
    # choose root cell as in Haghverdi et al. (2016)
    adata['iroot'] = iroot = 532 # note that in Matlab/R, counting starts at 1
    adata['xroot'] = adata.X[iroot]
    # annotate with Moignard et al. (2015) experimental cell groups
    adata['groups_names'] = groups_names = ['HF', 'NP', 'PS', '4SG', '4SFG']
    adata['groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    # get name for each cell
    adata.smp['groups'] = [
        next(gname for gname in groups_names if sname.startswith(gname))
        for sname in adata.smp_names]
    return adata

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
    adata = AnnData(ddata)
    # the data has to be transposed (in the hdf5 and R files, each row
    # corresponds to one gene, we use the opposite convention)
    adata = adata.transpose()
    # cluster assocations identified by Paul et al.
    # groups = sc.read(filename,'cluster.id')['X']
    infogenes_names = sc.read(filename, 'info.genes_strings')['X']
    # just keep the first of the two equivalent names per gene
    adata.var_names = np.array([gn.split(';')[0] for gn in adata.var_names])
    # index array for the informative genes
    infogenes_idcs = np.array([gn in infogenes_names for gn in adata.var_names])
    # restrict data array to the 3451 informative genes
    adata = adata[:, infogenes_idcs]
    # set root cell as in Haghverdi et al. (2016)
    adata['iroot'] = iroot = 840 # note that in Matlab/R, counting starts at 1
    adata['xroot'] = adata.X[iroot]
    return adata

def paul15_dpt(ddpt):
    ddpt['groups_names'] = ['', 'GMP', '', 'MEP']
    return ddpt

def paul15pca_dpt(ddpt):
    ddpt['groups_names'] = ['', '', 'GMP', 'MEP']
    return ddpt

def paul15_paths(dtool):
    dtool['groups_names'] = ['GMP', 'MEP']
    return dtool

def paul15pca_paths(dtool):
    dtool['groups_names'] = ['GMP', 'MEP']
    return dtool

