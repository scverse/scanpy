from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from tqdm.auto import tqdm

from .. import logging as logg, _utils
from .._compat import Literal
from .._settings import settings
from ..readwrite import read, read_visium

HERE = Path(__file__).parent


def blobs(
    n_variables: int = 11,
    n_centers: int = 5,
    cluster_std: float = 1.0,
    n_observations: int = 640,
) -> AnnData:
    """\
    Gaussian Blobs.

    Parameters
    ----------
    n_variables
        Dimension of feature space.
    n_centers
        Number of cluster centers.
    cluster_std
        Standard deviation of clusters.
    n_observations
        Number of observations. By default, this is the same observation number
        as in :func:`scanpy.datasets.krumsiek11`.

    Returns
    -------
    Annotated data matrix containing a observation annotation 'blobs' that
    indicates cluster identity.
    """
    import sklearn.datasets

    X, y = sklearn.datasets.make_blobs(
        n_samples=n_observations,
        n_features=n_variables,
        centers=n_centers,
        cluster_std=cluster_std,
        random_state=0,
    )
    return AnnData(X, obs=dict(blobs=y.astype(str)))


def burczynski06() -> AnnData:
    """\
    Bulk data with conditions ulcerative colitis (UC) and Crohn's disease (CD).

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
    adata = read(filename, backup_url=url)
    return adata


def krumsiek11() -> AnnData:
    """\
    Simulated myeloid progenitors [Krumsiek11]_.

    The literature-curated boolean network from [Krumsiek11]_ was used to
    simulate the data. It describes development to four cell fates: 'monocyte',
    'erythrocyte', 'megakaryocyte' and 'neutrophil'.

    See also the discussion of this data in [Wolf19]_.

    Simulate via :func:`~scanpy.tl.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = HERE / 'krumsiek11.txt'
    verbosity_save = settings.verbosity
    settings.verbosity = 'error'  # suppress output...
    adata = read(filename, first_column_names=True)
    settings.verbosity = verbosity_save
    adata.uns['iroot'] = 0
    fate_labels = {0: 'Stem', 159: 'Mo', 319: 'Ery', 459: 'Mk', 619: 'Neu'}
    adata.uns['highlights'] = fate_labels
    cell_type = np.array(['progenitor' for i in range(adata.n_obs)])
    cell_type[80:160] = 'Mo'
    cell_type[240:320] = 'Ery'
    cell_type[400:480] = 'Mk'
    cell_type[560:640] = 'Neu'
    adata.obs['cell_type'] = cell_type
    _utils.sanitize_anndata(adata)
    return adata


def moignard15() -> AnnData:
    """\
    Hematopoiesis in early mouse embryos [Moignard15]_.

    Returns
    -------
    Annotated data matrix.
    """
    filename = settings.datasetdir / 'moignard15/nbt.3154-S3.xlsx'
    backup_url = 'http://www.nature.com/nbt/journal/v33/n3/extref/nbt.3154-S3.xlsx'
    adata = read(filename, sheet='dCt_values.txt', backup_url=backup_url)
    # filter out 4 genes as in Haghverdi et al. (2016)
    gene_subset = ~np.in1d(adata.var_names, ['Eif2b1', 'Mrpl19', 'Polr2a', 'Ubc'])
    adata = adata[:, gene_subset]  # retain non-removed genes
    # choose root cell for DPT analysis as in Haghverdi et al. (2016)
    adata.uns["iroot"] = 532  # note that in Matlab/R, counting starts at 1
    # annotate with Moignard et al. (2015) experimental cell groups
    groups = {
        'HF': '#D7A83E',
        'NP': '#7AAE5D',
        'PS': '#497ABC',
        '4SG': '#AF353A',
        '4SFG': '#765099',
    }
    # annotate each observation/cell
    adata.obs['exp_groups'] = [
        next(gname for gname in groups.keys() if sname.startswith(gname))
        for sname in adata.obs_names
    ]
    # fix the order and colors of names in "groups"
    adata.obs['exp_groups'] = pd.Categorical(
        adata.obs['exp_groups'], categories=list(groups.keys())
    )
    adata.uns['exp_groups_colors'] = list(groups.values())
    return adata


def paul15() -> AnnData:
    """\
    Development of Myeloid Progenitors [Paul15]_.

    Non-logarithmized raw data.

    The data has been sent out by Email from the Amit Lab. An R version for
    loading the data can be found here
    https://github.com/theislab/scAnalysisTutorial

    Returns
    -------
    Annotated data matrix.
    """
    logg.warning(
        'In Scanpy 0.*, this returned logarithmized data. '
        'Now it returns non-logarithmized data.'
    )
    import h5py

    filename = settings.datasetdir / 'paul15/paul15.h5'
    backup_url = 'http://falexwolf.de/data/paul15.h5'
    _utils.check_presence_download(filename, backup_url)
    with h5py.File(filename, 'r') as f:
        X = f['data.debatched'][()]
        gene_names = f['data.debatched_rownames'][()].astype(str)
        cell_names = f['data.debatched_colnames'][()].astype(str)
        clusters = f['cluster.id'][()].flatten().astype(int)
        infogenes_names = f['info.genes_strings'][()].astype(str)
    # each row has to correspond to a observation, therefore transpose
    adata = AnnData(X.transpose())
    adata.var_names = gene_names
    adata.row_names = cell_names
    # names reflecting the cell type identifications from the paper
    cell_type = 6 * ['Ery']
    cell_type += 'MEP Mk GMP GMP DC Baso Baso Mo Mo Neu Neu Eos Lymph'.split()
    adata.obs['paul15_clusters'] = [f'{i}{cell_type[i-1]}' for i in clusters]
    # make string annotations categorical (optional)
    _utils.sanitize_anndata(adata)
    # just keep the first of the two equivalent names per gene
    adata.var_names = [gn.split(';')[0] for gn in adata.var_names]
    # remove 10 corrupted gene names
    infogenes_names = np.intersect1d(infogenes_names, adata.var_names)
    # restrict data array to the 3461 informative genes
    adata = adata[:, infogenes_names]
    # usually we'd set the root cell to an arbitrary cell in the MEP cluster
    # adata.uns['iroot'] = np.flatnonzero(adata.obs['paul15_clusters'] == '7MEP')[0]
    # here, set the root cell as in Haghverdi et al. (2016)
    # note that other than in Matlab/R, counting starts at 0
    adata.uns['iroot'] = 840
    return adata


def toggleswitch() -> AnnData:
    """\
    Simulated toggleswitch.

    Data obtained simulating a simple toggleswitch [Gardner00]_

    Simulate via :func:`~scanpy.tl.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = HERE / 'toggleswitch.txt'
    adata = read(filename, first_column_names=True)
    adata.uns['iroot'] = 0
    return adata


def pbmc68k_reduced() -> AnnData:
    """\
    Subsampled and processed 68k PBMCs.

    10x PBMC 68k dataset from
    https://support.10xgenomics.com/single-cell-gene-expression/datasets

    The original PBMC 68k dataset was preprocessed using scanpy and was saved
    keeping only 724 cells and 221 highly variable genes.

    The saved file contains the annotation of cell types (key: `'bulk_labels'`),
    UMAP coordinates, louvain clustering and gene rankings based on the
    `bulk_labels`.

    Returns
    -------
    Annotated data matrix.
    """

    filename = HERE / '10x_pbmc68k_reduced.h5ad'
    return read(filename)


def pbmc3k() -> AnnData:
    """\
    3k PBMCs from 10x Genomics.

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
            # the directory with the `.mtx` file
            './data/filtered_gene_bc_matrices/hg19/',
            # use gene symbols for the variable names (variables-axis index)
            var_names='gene_symbols',
            # write a cache file for faster subsequent reading
            cache=True,
        )

        adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
        adata.write('write/pbmc3k_raw.h5ad', compression='gzip')

    Returns
    -------
    Annotated data matrix.
    """
    url = 'http://falexwolf.de/data/pbmc3k_raw.h5ad'
    adata = read(settings.datasetdir / 'pbmc3k_raw.h5ad', backup_url=url)
    return adata


def pbmc3k_processed() -> AnnData:
    """Processed 3k PBMCs from 10x Genomics.

    Processed using the `basic tutorial <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>`__.

    Returns
    -------
    Annotated data matrix.
    """
    adata = read(
        settings.datasetdir / 'pbmc3k_processed.h5ad',
        backup_url='https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad',
    )
    return adata


def _download_visium_dataset(sample_id: str, base_dir: Optional[Path] = None):
    """
    Params
    ------
    sample_id
        String name of example visium dataset.
    base_dir
        Where to download the dataset to.
    """
    import tarfile

    if base_dir is None:
        base_dir = settings.datasetdir

    url_prefix = f'http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/{sample_id}/'

    sample_dir = base_dir / sample_id
    sample_dir.mkdir(exist_ok=True)

    # Download spatial data
    tar_filename = f"{sample_id}_spatial.tar.gz"
    tar_pth = sample_dir / tar_filename
    _utils.check_presence_download(
        filename=tar_pth, backup_url=url_prefix + tar_filename
    )
    with tarfile.open(tar_pth) as f:
        for el in f:
            if not (sample_dir / el.name).exists():
                f.extract(el, sample_dir)

    # Download counts
    _utils.check_presence_download(
        filename=sample_dir / "filtered_feature_bc_matrix.h5",
        backup_url=url_prefix + f"{sample_id}_filtered_feature_bc_matrix.h5",
    )


def visium_sge(
    sample_id: Literal[
        'V1_Breast_Cancer_Block_A_Section_1',
        'V1_Breast_Cancer_Block_A_Section_2',
        'V1_Human_Heart',
        'V1_Human_Lymph_Node',
        'V1_Mouse_Kidney',
        'V1_Adult_Mouse_Brain',
        'V1_Mouse_Brain_Sagittal_Posterior',
        'V1_Mouse_Brain_Sagittal_Posterior_Section_2',
        'V1_Mouse_Brain_Sagittal_Anterior',
        'V1_Mouse_Brain_Sagittal_Anterior_Section_2',
    ] = 'V1_Breast_Cancer_Block_A_Section_1',
) -> AnnData:
    """\
    Processed Visium Spatial Gene Expression data from 10x Genomics.
    Database: https://support.10xgenomics.com/spatial-gene-expression/datasets

    Parameters
    ----------
    sample_id
        The ID of the data sample in 10x’s spatial database.

    Returns
    -------
    Annotated data matrix.
    """
    _download_visium_dataset(sample_id)
    return read_visium(settings.datasetdir / sample_id)
