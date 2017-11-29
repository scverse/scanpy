# Author: Alex Wolf (http://falexwolf.de)
"""Builtin Datasets.
"""

import numpy as np
from . import api_without_examples as sc


def blobs(n_centers=5, cluster_std=1.0, n_samples=640):
    """Gaussian Blobs.

    Parameters
    ----------
    n_centers : `int`, optional (default: 5)
        Number of cluster centers.
    cluster_std : `float`, optional (default: 1.0)
        Standard deviation of clusters.
    n_samples : `int`, optional (default: 640)
        Number of samples. By default, this is the same sample number as in
        ``sc.examples.krumsiek11()``.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix containing a sample annotation 'blobs' that
        indicates cluster identity.
    """
    import sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=n_samples,
                                       n_features=11,
                                       centers=n_centers,
                                       cluster_std=cluster_std,
                                       random_state=0)
    return sc.AnnData(X, smp={'blobs': y.astype(str)})


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
    """Simulated myeloid progenitors [Krumsiek11]_.

    The literature-curated boolean network from [Krumsiek11]_ was used to
    simulate the data. It describes development to four cell fates: 'monocyte',
    'erythrocyte', 'megakaryocyte' and 'neutrophil'.

    Simulate via :func:`~scanpy.api.sim`.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
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
    adata.uns['iroot'] = 0
    fate_labels = {0: 'progenitor', 159: 'monocyte', 319: 'erythrocyte',
                   459: 'megakaryocyte', 619: 'neutrophil'}
    adata.uns['highlights'] = fate_labels
    cell_type = np.array(['progenitor' for i in range(adata.n_smps)])
    cell_type[80:160] = 'monocyte'
    cell_type[240:320] = 'erythrocyte'
    cell_type[400:480] = 'megakaryocyte'
    cell_type[560:640] = 'neutrophil'
    adata.smp['cell_type'] = cell_type
    return adata


krumsiek11_diffmap_params = {'n_neighbors': 5, 'knn': False}
krumsiek11_dpt_params = {'n_neighbors': 5, 'knn': False, 'n_branchings': 2}


def moignard15():
    """Hematopoiesis in early mouse embryos [Moignard15]_.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    """
    filename = 'data/moignard15/nbt.3154-S3.xlsx'
    backup_url = 'http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3154-S3.xlsx'
    adata = sc.read(filename, sheet='dCt_values.txt', cache=True, backup_url=backup_url)
    # filter out 4 genes as in Haghverdi et al. (2016)
    gene_subset = ~np.in1d(adata.var_names, ['Eif2b1', 'Mrpl19', 'Polr2a', 'Ubc'])
    adata = adata[:, gene_subset]  # retain non-removed genes
    # choose root cell for DPT analysis as in Haghverdi et al. (2016)
    adata.uns['iroot'] = 532  # note that in Matlab/R, counting starts at 1
    # annotate with Moignard et al. (2015) experimental cell groups
    groups_order = ['HF', 'NP', 'PS', '4SG', '4SFG']
    # annotate each sample/cell
    adata.smp['exp_groups'] = [
        next(gname for gname in groups_order if sname.startswith(gname))
        for sname in adata.smp_names]
    # fix the order and colors of names in "groups"
    adata.uns['exp_groups_order'] = groups_order
    adata.uns['exp_groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    return adata


moignard15_diffmap_params = {'n_neighbors': 5, 'knn': False}
moignard15_dpt_params = {'n_neighbors': 5, 'knn': False}


def moignard15_dpt(adata):
    """Add some labeling information to DPT result.
    """
    sc.logg.m('... adding annotation for DPT groups')
    if len(adata.uns['dpt_groups_order']) > 1:
        groups_order = ['undecided', 'endothelial',
                        'erythrocytes', 'trunk']
        adata.uns['dpt_groups_order'] = ['{}: {}'.format(i, n) for i, n in enumerate(groups_order)]
    return adata


def paul15():
    """Development of Myeloid Progenitors [Paul15]_.

    Logarithmized raw data.

    The data has been sent out by Email from the Amit Lab. An R version for
    loading the data can be found here
    https://github.com/theislab/scAnalysisTutorial

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    """
    adata = paul15_raw()
    sc.pp.log1p(adata)
    return adata

paul15_diffmap_params = {'n_neighbors': 20, 'n_pcs': 0}
paul15_dpt_params = {'n_neighbors': 20, 'n_pcs': 0}


def paul15_raw():
    """Development of Myeloid Progenitors [Paul15]_.

    Non-logarithmized raw data.

    The data has been sent out by Email from the Amit Lab. An R version for
    loading the data can be found here
    https://github.com/theislab/scAnalysisTutorial

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    """
    import h5py
    filename = 'data/paul15/paul15.h5'
    backup_url = 'http://falexwolf.de/data/paul15.h5'
    sc.utils.check_presence_download(filename, backup_url)
    with h5py.File(filename, 'r') as f:
        X = f['data.debatched'][()]
        gene_names = f['data.debatched_rownames'][()].astype(str)
        cell_names = f['data.debatched_colnames'][()].astype(str)
        clusters = f['cluster.id'][()].flatten()
        infogenes_names = f['info.genes_strings'][()].astype(str)
    # each row has to correspond to a sample, therefore transpose
    adata = sc.AnnData(X.transpose())
    adata.var_names = gene_names
    adata.row_names = cell_names
    # names reflecting the cell type identifications from the paper
    cell_type = {7: 'MEP', 8: 'Mk', 9: 'GMP', 10: 'GMP', 11: 'DC',
                 12: 'Baso', 13: 'Baso', 14: 'Mo', 15: 'Mo',
                 16: 'Neu', 17: 'Neu', 18: 'Eos', 19: 'Lymph'}
    cell_type.update({i: 'Ery' for i in range(1, 7)})
    adata.smp['paul15_clusters'] = [
        str(i) + cell_type[i] for i in clusters.astype(int)]
    # make string annotations categorical (optional)
    sc.utils.sanitize_anndata(adata)
    # just keep the first of the two equivalent names per gene
    adata.var_names = [gn.split(';')[0] for gn in adata.var_names]
    # remove 10 corrupted gene names
    infogenes_names = np.intersect1d(infogenes_names, adata.var_names)
    # restrict data array to the 3461 informative genes
    adata = adata[:, infogenes_names]
    # usually we'd set the root cell to an arbitrary cell in the MEP cluster
    # adata.uns['iroot': np.flatnonzero(adata.smp['paul15_clusters']  == '7MEP')[0]
    # here, set the root cell as in Haghverdi et al. (2016)
    adata.uns['iroot'] = 840  # note that other than in Matlab/R, counting starts at 1
    return adata


def paul15_dpt(adata):
    """Post-processing for DPT."""
    adata.uns['dpt_groups_order'] = ['', 'GMP', '', 'MEP']


def toggleswitch():
    """Simulated toggleswitch.

    Data obtained simulating a simple toggleswitch `Gardner *et al.*, Nature
    (2000) <https://doi.org/10.1038/35002131>`_.

    Simulate via :func:`~scanpy.api.sim`.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    """
    import os
    from .. import logging as logg
    filename = 'write/toggleswitch_sim/sim_000000.txt'
    if not os.path.exists(filename):
        filename = os.path.dirname(__file__) + '/toggleswitch.txt'
        logg.hint('You can reproduce the data file {} '
                  'by running `sc.tl.sim("toggleswitch")`.'
                  .format(filename))
    adata = sc.read(filename, first_column_names=True, cache=True)
    adata.uns['iroot'] = 0
    return adata


toggleswitch_diffmap_params = {'n_neighbors': 5, 'knn': False}
toggleswitch_dpt_params = {'n_neighbors': 5, 'knn': False}
