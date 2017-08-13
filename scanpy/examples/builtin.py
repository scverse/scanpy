"""Builtin Examples.

Provides functions for preprocessing data and default parameters.

Defines command-line "runs" for builtin examples.
"""

import numpy as np
from . import api_without_examples as sc


def blobs(centers=5, cluster_std=1.0):
    """Make Gaussian Blobs.

    Same sample number as in krumsiek11, to compare with the latter.
    """
    import sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=640,
                                       n_features=11,
                                       centers=centers,
                                       cluster_std=cluster_std,
                                       random_state=0)
    return sc.AnnData(X, smp={'clusters': y.astype(str)})


def burczynski06():
    """Bulk data with conditions ulcerative colitis (UC) and Crohn's disease (CD).

    The study assesses transcriptional profiles in peripheral blood mononuclear
    cells from 42 healthy individuals, 59 CD patients, and 26 UC patients by
    hybridization to microarrays interrogating more than 22,000 sequences.

    Reference
    ---------
    Burczynski et al., "Molecular classification of Crohn's disease and
    ulcerative colitis patients using transcriptional profiles in peripheral
    blood mononuclear cells"
    J Mol Diagn 8, 51 (2006). PMID:16436634.
    """
    filename = 'data/burczynski06/GDS1615_full.soft.gz'
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz'
    adata = sc.read(filename, backup_url=url, cache=True)
    return adata


def krumsiek11():
    """Simulated myeloid progenitor data.

    A literature-curated boolean network from the reference below was
    used to simulate this data with four cell fates.

    Simulate the data by running "scanpy krumsiek11 sim" on the command line
    or running `sc.tl.sim('krumsiek11')`.

    Reference
    ---------
    Krumsiek et al., "Hierarchical Differentiation of Myeloid Progenitors Is
    Encoded in the Transcription Factor Network"
    PLoS ONE 6, e22649 (2011).
    """
    import os
    from .. import logging as logg
    filename = 'write/krumsiek11_sim/sim_000000.txt'
    if not os.path.exists(filename):
        filename = os.path.dirname(__file__) + '/krumsiek11.txt'
        logg.hint('you can reproduce the data file {} '
                  'by running `sc.tl.sim("krumsiek11")`'
                  .format(filename))
    adata = sc.read(filename, first_column_names=True, cache=True)
    adata.add['iroot'] = 0
    return adata


krumsiek11_diffmap_params = {'n_neighbors': 5, 'knn': False}
krumsiek11_dpt_params = {'n_neighbors': 5, 'knn': False, 'n_branchings': 2}


def moignard15():
    """Hematopoiesis in early mouse embryos.

    1. Filter out a few genes.
    2. Choose 'root cell'.
    3. Define experimental groups by cropping cell names.

    Reference
    ---------
    Moignard et al., "Decoding the regulatory network of early blood development
    from single-cell gene expression measurements"
    Nature Biotechnology 33, 269 (2015)
    """
    filename = 'data/moignard15/nbt.3154-S3.xlsx'
    backup_url = 'http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3154-S3.xlsx'
    adata = sc.read(filename, sheet='dCt_values.txt', cache=True, backup_url=backup_url)
    # filter out 4 genes as in Haghverdi et al. (2016)
    gene_subset = ~np.in1d(adata.var_names, ['Eif2b1', 'Mrpl19', 'Polr2a', 'Ubc'])
    adata = adata[:, gene_subset]  # retain non-removed genes
    # choose root cell for DPT analysis as in Haghverdi et al. (2016)
    adata.add['iroot'] = 532  # note that in Matlab/R, counting starts at 1
    # annotate with Moignard et al. (2015) experimental cell groups
    groups_names = ['HF', 'NP', 'PS', '4SG', '4SFG']
    # annotate each sample/cell
    adata.smp['exp_groups'] = [
        next(gname for gname in groups_names if sname.startswith(gname))
        for sname in adata.smp_names]
    # fix the order and colors of names in "groups"
    adata.add['exp_groups_names'] = groups_names
    adata.add['exp_groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    return adata


moignard15_diffmap_params = {'n_neighbors': 5, 'knn': False}
moignard15_dpt_params = {'n_neighbors': 5, 'knn': False}


def moignard15_dpt(adata):
    """Add some labeling information to DPT result.
    """
    sc.logg.m('... adding annotation for DPT groups')
    if len(adata.add['dpt_groups_names']) > 1:
        groups_names = ['undecided', 'endothelial',
                        'erythrocytes', 'trunk']
        adata.add['dpt_groups_names'] = ['{}: {}'.format(i, n) for i, n in enumerate(groups_names)]
    return adata


def paul15():
    """Myeloid Progenitors.

    This largely follows an R tutorial by Maren Buttner.
    https://github.com/theislab/scAnalysisTutorial

    Reference
    ---------
    Paul et al., "Transcriptional Heterogeneity and Lineage Commitment in
    Myeloid Progenitors",
    Cell 163, 1663 (2015)
    """
    adata = paul15_raw()
    sc.pp.log1p(adata)
    adata.add['xroot'] = adata.X[adata.add['iroot']]  # adjust expression vector of root cell
    return adata

paul15_diffmap_params = {'n_neighbors': 20, 'n_pcs': 0}
paul15_dpt_params = {'n_neighbors': 20, 'n_pcs': 0}


def paul15pca():
    """Same as paul15.
    """
    adata = paul15_raw()
    sc.pp.log1p(adata)
    adata.add['xroot'] = adata.X[adata.add['iroot']]  # adjust expression vector of root cell
    return adata


paul15pca_dpt_params = {'n_neighbors': 20, 'knn': True}


def paul15_raw():
    filename = 'data/paul15/paul15.h5'
    backup_url = 'http://falexwolf.de/data/paul15.h5'
    adata = sc.read(filename, 'data.debatched', backup_url=backup_url)
    # each row has to correspond to a sample, therefore transpose
    adata = adata.transpose()
    # cluster assocations identified by Paul et al.
    # groups = sc.read(filename, 'cluster.id', return_dict=True)['X']
    infogenes_names = sc.read(filename, 'info.genes_strings', return_dict=True)['X']
    # just keep the first of the two equivalent names per gene
    adata.var_names = np.array([gn.split(';')[0] for gn in adata.var_names])
    # remove 10 corrupted gene names
    infogenes_names = np.intersect1d(infogenes_names, adata.var_names)
    # restrict data array to the 3461 informative genes
    adata = adata[:, infogenes_names]
    # set root cell as in Haghverdi et al. (2016)
    adata.add['iroot'] = iroot = 840  # note that other than in Matlab/R, counting starts at 1
    adata.add['xroot'] = adata.X[iroot]
    return adata


def paul15_dpt(adata):
    adata.add['dpt_groups_names'] = ['', 'GMP', '', 'MEP']


def paul15pca_dpt(adata):
    adata.add['dpt_groups_names'] = ['', '', 'GMP', 'MEP']


def toggleswitch():
    """Simple toggleswitch from simulated data.

    Simulate the data by running "scanpy sim toggleswitch" on the command line
    or `sc.tl.sim("toggleswitch")`.
    """
    import os
    from .. import logging as logg
    filename = 'write/toggleswitch_sim/sim_000000.txt'
    if not os.path.exists(filename):
        filename = os.path.dirname(__file__) + '/toggleswitch.txt'
        logg.hint('you can reproduce the data file {} '
                  'by running `sc.tl.sim("toggleswitch")`'
                  .format(filename))
    adata = sc.read(filename, first_column_names=True, cache=True)
    adata.add['xroot'] = adata.X[0]
    return adata


toggleswitch_diffmap_params = {'n_neighbors': 5, 'knn': False}
toggleswitch_dpt_params = {'n_neighbors': 5, 'knn': False}
