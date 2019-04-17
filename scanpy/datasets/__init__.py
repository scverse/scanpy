"""Builtin Datasets.
"""

import os
import numpy as np
import pandas as pd
from anndata import AnnData

from .. import logging as logg
from .._settings import settings
import scanpy as sc
from ._ebi_expression_atlas import ebi_expression_atlas


def blobs(n_variables=11, n_centers=5, cluster_std=1.0, n_observations=640) -> AnnData:
    """Gaussian Blobs.

    Parameters
    ----------
    n_variables : `int`, optional (default: 11)
        Dimension of feature space.
    n_centers : `int`, optional (default: 5)
        Number of cluster centers.
    cluster_std : `float`, optional (default: 1.0)
        Standard deviation of clusters.
    n_observations : `int`, optional (default: 640)
        Number of observations. By default, this is the same observation number as in
        ``sc.datasets.krumsiek11()``.

    Returns
    -------
    Annotated data matrix containing a observation annotation 'blobs' that
    indicates cluster identity.
    """
    import sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=n_observations,
                                       n_features=n_variables,
                                       centers=n_centers,
                                       cluster_std=cluster_std,
                                       random_state=0)
    return AnnData(X, obs={'blobs': y.astype(str)})


def burczynski06() -> AnnData:
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
    filename = settings.datasetdir / 'burczynski06/GDS1615_full.soft.gz'
    url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz'
    adata = sc.read(filename, backup_url=url)
    return adata


def krumsiek11() -> AnnData:
    """Simulated myeloid progenitors [Krumsiek11]_.

    The literature-curated boolean network from [Krumsiek11]_ was used to
    simulate the data. It describes development to four cell fates: 'monocyte',
    'erythrocyte', 'megakaryocyte' and 'neutrophil'.

    See also the discussion of this data in [Wolf19]_.

    Simulate via :func:`~scanpy.api.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = os.path.dirname(__file__) + '/krumsiek11.txt'
    verbosity_save = sc.settings.verbosity
    sc.settings.verbosity = 'error'  # suppress output...
    adata = sc.read(filename, first_column_names=True)
    sc.settings.verbosity = verbosity_save
    adata.uns['iroot'] = 0
    fate_labels = {0: 'Stem', 159: 'Mo', 319: 'Ery',
                   459: 'Mk', 619: 'Neu'}
    adata.uns['highlights'] = fate_labels
    cell_type = np.array(['progenitor' for i in range(adata.n_obs)])
    cell_type[80:160] = 'Mo'
    cell_type[240:320] = 'Ery'
    cell_type[400:480] = 'Mk'
    cell_type[560:640] = 'Neu'
    adata.obs['cell_type'] = cell_type
    sc.utils.sanitize_anndata(adata)
    return adata


def moignard15() -> AnnData:
    """Hematopoiesis in early mouse embryos [Moignard15]_.

    Returns
    -------
    Annotated data matrix.
    """
    filename = settings.datasetdir / 'moignard15/nbt.3154-S3.xlsx'
    backup_url = 'http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3154-S3.xlsx'
    adata = sc.read(filename, sheet='dCt_values.txt', backup_url=backup_url)
    # filter out 4 genes as in Haghverdi et al. (2016)
    gene_subset = ~np.in1d(adata.var_names, ['Eif2b1', 'Mrpl19', 'Polr2a', 'Ubc'])
    adata = adata[:, gene_subset]  # retain non-removed genes
    # choose root cell for DPT analysis as in Haghverdi et al. (2016)
    adata.uns['iroot'] = 532  # note that in Matlab/R, counting starts at 1
    # annotate with Moignard et al. (2015) experimental cell groups
    groups_order = ['HF', 'NP', 'PS', '4SG', '4SFG']
    # annotate each observation/cell
    adata.obs['exp_groups'] = [
        next(gname for gname in groups_order if sname.startswith(gname))
        for sname in adata.obs_names]
    # fix the order and colors of names in "groups"
    adata.obs['exp_groups'] = pd.Categorical(adata.obs['exp_groups'],
                                             categories=groups_order)
    adata.uns['exp_groups_colors'] = ['#D7A83E', '#7AAE5D', '#497ABC', '#AF353A', '#765099']
    return adata


def paul15() -> AnnData:
    """Development of Myeloid Progenitors [Paul15]_.

    Non-logarithmized raw data.

    The data has been sent out by Email from the Amit Lab. An R version for
    loading the data can be found here
    https://github.com/theislab/scAnalysisTutorial

    Returns
    -------
    Annotated data matrix.
    """
    logg.warn('In Scanpy 0.*, this returned logarithmized data. '
              'Now it returns non-logarithmized data.')
    import h5py
    filename = settings.datasetdir / 'paul15/paul15.h5'
    backup_url = 'http://falexwolf.de/data/paul15.h5'
    sc.utils.check_presence_download(filename, backup_url)
    with h5py.File(filename, 'r') as f:
        X = f['data.debatched'][()]
        gene_names = f['data.debatched_rownames'][()].astype(str)
        cell_names = f['data.debatched_colnames'][()].astype(str)
        clusters = f['cluster.id'][()].flatten()
        infogenes_names = f['info.genes_strings'][()].astype(str)
    # each row has to correspond to a observation, therefore transpose
    adata = AnnData(X.transpose())
    adata.var_names = gene_names
    adata.row_names = cell_names
    # names reflecting the cell type identifications from the paper
    cell_type = {7: 'MEP', 8: 'Mk', 9: 'GMP', 10: 'GMP', 11: 'DC',
                 12: 'Baso', 13: 'Baso', 14: 'Mo', 15: 'Mo',
                 16: 'Neu', 17: 'Neu', 18: 'Eos', 19: 'Lymph'}
    cell_type.update({i: 'Ery' for i in range(1, 7)})
    adata.obs['paul15_clusters'] = [
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
    # adata.uns['iroot': np.flatnonzero(adata.obs['paul15_clusters']  == '7MEP')[0]
    # here, set the root cell as in Haghverdi et al. (2016)
    adata.uns['iroot'] = 840  # note that other than in Matlab/R, counting starts at 1
    return adata


def toggleswitch() -> AnnData:
    """Simulated toggleswitch.

    Data obtained simulating a simple toggleswitch `Gardner *et al.*, Nature
    (2000) <https://doi.org/10.1038/35002131>`__.

    Simulate via :func:`~scanpy.api.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = os.path.dirname(__file__) + '/toggleswitch.txt'
    adata = sc.read(filename, first_column_names=True)
    adata.uns['iroot'] = 0
    return adata


def pbmc68k_reduced() -> AnnData:
    """Subsampled and processed 68k PBMCs.

    10x PBMC 68k dataset from
    https://support.10xgenomics.com/single-cell-gene-expression/datasets

    The original PBMC 68k dataset was preprocessed using scanpy and was saved
    keeping only 724 cells and 221 highly variable genes.

    The saved file contains the annotation of cell types (key: 'bulk_labels'), UMAP coordinates,
    louvain clustering and gene rankings based on the bulk_labels.

    Returns
    -------
    Annotated data matrix.
    """

    filename = os.path.dirname(__file__) + '/10x_pbmc68k_reduced.h5ad'
    return sc.read(filename)


def pbmc3k() -> AnnData:
    """3k PBMCs from 10x Genomics.

    The data consists in 3k PBMCs from a Healthy Donor and is freely available
    from 10x Genomics (`here
    <http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>`__
    from this `webpage
    <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>`__).

    The exact same data is also used in Seurat's
    `basic clustering tutorial <https://satijalab.org/seurat/pbmc3k_tutorial.html>`__.

    .. note::

        This downloads 5.9 MB of data upon the first call of the function and stores it in `./data/pbmc3k_raw.h5ad`.

    The following code was run to produce the file.

    .. code:: python

        adata = sc.read_10x_mtx(
        './data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
        var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
        cache=True)                                # write a cache file for faster subsequent reading

        adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
        adata.write('write/pbmc3k_raw.h5ad', compression='gzip')

    Returns
    -------
    Annotated data matrix.
    """
    adata = sc.read(settings.datasetdir / 'pbmc3k_raw.h5ad', backup_url='http://falexwolf.de/data/pbmc3k_raw.h5ad')
    return adata
