from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from anndata import AnnData

from .. import _utils
from .. import logging as logg
from .._compat import old_positionals
from .._settings import settings
from ..readwrite import read, read_visium
from ._utils import check_datasetdir_exists, filter_oldformatwarning

if TYPE_CHECKING:
    from .._utils import AnyRandom

HERE = Path(__file__).parent


@old_positionals(
    "n_variables", "n_centers", "cluster_std", "n_observations", "random_state"
)
def blobs(
    *,
    n_variables: int = 11,
    n_centers: int = 5,
    cluster_std: float = 1.0,
    n_observations: int = 640,
    random_state: AnyRandom = 0,
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
    random_state
        Determines random number generation for dataset creation.

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
        random_state=random_state,
    )
    return AnnData(X, obs=dict(blobs=y.astype(str)))


@check_datasetdir_exists
def burczynski06() -> AnnData:
    """\
    Bulk data with conditions ulcerative colitis (UC) and Crohn’s disease (CD) :cite:p:`Burczynski2006`.

    The study assesses transcriptional profiles in peripheral blood mononuclear
    cells from 42 healthy individuals, 59 CD patients, and 26 UC patients by
    hybridization to microarrays interrogating more than 22,000 sequences.
    """
    filename = settings.datasetdir / "burczynski06/GDS1615_full.soft.gz"
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz"
    return read(filename, backup_url=url)


def krumsiek11() -> AnnData:
    """\
    Simulated myeloid progenitors :cite:p:`Krumsiek2011`.

    The literature-curated boolean network from :cite:t:`Krumsiek2011` was used to
    simulate the data. It describes development to four cell fates: 'monocyte',
    'erythrocyte', 'megakaryocyte' and 'neutrophil'.

    See also the discussion of this data in :cite:t:`Wolf2019`.

    Simulate via :func:`~scanpy.tl.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    with settings.verbosity.override("error"):  # suppress output...
        adata = read(HERE / "krumsiek11.txt", first_column_names=True)
    adata.uns["iroot"] = 0
    fate_labels = {0: "Stem", 159: "Mo", 319: "Ery", 459: "Mk", 619: "Neu"}
    adata.uns["highlights"] = fate_labels
    cell_type = pd.array(["progenitor"]).repeat(adata.n_obs)
    cell_type[80:160] = "Mo"
    cell_type[240:320] = "Ery"
    cell_type[400:480] = "Mk"
    cell_type[560:640] = "Neu"
    adata.obs["cell_type"] = cell_type
    _utils.sanitize_anndata(adata)
    return adata


@check_datasetdir_exists
def moignard15() -> AnnData:
    """\
    Hematopoiesis in early mouse embryos :cite:p:`Moignard2015`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = settings.datasetdir / "moignard15/nbt.3154-S3.xlsx"
    backup_url = "https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3154/MediaObjects/41587_2015_BFnbt3154_MOESM4_ESM.xlsx"
    adata = read(filename, sheet="dCt_values.txt", backup_url=backup_url)
    # filter out 4 genes as in Haghverdi et al. (2016)
    gene_subset = ~np.in1d(adata.var_names, ["Eif2b1", "Mrpl19", "Polr2a", "Ubc"])
    adata = adata[:, gene_subset].copy()  # retain non-removed genes
    # choose root cell for DPT analysis as in Haghverdi et al. (2016)
    adata.uns["iroot"] = 532  # note that in Matlab/R, counting starts at 1
    # annotate with Moignard et al. (2015) experimental cell groups
    groups = {
        "HF": "#D7A83E",
        "NP": "#7AAE5D",
        "PS": "#497ABC",
        "4SG": "#AF353A",
        "4SFG": "#765099",
    }
    # annotate each observation/cell
    adata.obs["exp_groups"] = [
        next(gname for gname in groups.keys() if sname.startswith(gname))
        for sname in adata.obs_names
    ]
    # fix the order and colors of names in "groups"
    adata.obs["exp_groups"] = pd.Categorical(
        adata.obs["exp_groups"], categories=list(groups.keys())
    )
    adata.uns["exp_groups_colors"] = list(groups.values())
    return adata


@check_datasetdir_exists
def paul15() -> AnnData:
    """\
    Development of Myeloid Progenitors :cite:p:`Paul2015`.

    Non-logarithmized raw data.

    The data has been sent out by Email from the Amit Lab. An R version for
    loading the data can be found here
    https://github.com/theislab/scAnalysisTutorial

    Returns
    -------
    Annotated data matrix.
    """
    logg.warning(
        "In Scanpy 0.*, this returned logarithmized data. "
        "Now it returns non-logarithmized data."
    )
    import h5py

    filename = settings.datasetdir / "paul15/paul15.h5"
    filename.parent.mkdir(exist_ok=True)
    backup_url = "https://falexwolf.de/data/paul15.h5"
    _utils.check_presence_download(filename, backup_url)
    with h5py.File(filename, "r") as f:
        # Coercing to float32 for backwards compatibility
        X = f["data.debatched"][()].astype(np.float32)
        gene_names = f["data.debatched_rownames"][()].astype(str)
        cell_names = f["data.debatched_colnames"][()].astype(str)
        clusters = f["cluster.id"][()].flatten().astype(int)
        infogenes_names = f["info.genes_strings"][()].astype(str)
    # each row has to correspond to a observation, therefore transpose
    adata = AnnData(X.transpose())
    adata.var_names = gene_names
    adata.obs_names = cell_names
    # names reflecting the cell type identifications from the paper
    cell_type = 6 * ["Ery"]
    cell_type += "MEP Mk GMP GMP DC Baso Baso Mo Mo Neu Neu Eos Lymph".split()
    adata.obs["paul15_clusters"] = [f"{i}{cell_type[i - 1]}" for i in clusters]
    # make string annotations categorical (optional)
    _utils.sanitize_anndata(adata)
    # just keep the first of the two equivalent names per gene
    adata.var_names = [gn.split(";")[0] for gn in adata.var_names]
    # remove 10 corrupted gene names
    infogenes_names = np.intersect1d(infogenes_names, adata.var_names)
    # restrict data array to the 3461 informative genes
    adata = adata[:, infogenes_names].copy()
    # usually we'd set the root cell to an arbitrary cell in the MEP cluster
    # adata.uns['iroot'] = np.flatnonzero(adata.obs['paul15_clusters'] == '7MEP')[0]
    # here, set the root cell as in Haghverdi et al. (2016)
    # note that other than in Matlab/R, counting starts at 0
    adata.uns["iroot"] = 840
    return adata


def toggleswitch() -> AnnData:
    """\
    Simulated toggleswitch.

    Data obtained simulating a simple toggleswitch :cite:p:`Gardner2000`

    Simulate via :func:`~scanpy.tl.sim`.

    Returns
    -------
    Annotated data matrix.
    """
    filename = HERE / "toggleswitch.txt"
    adata = read(filename, first_column_names=True)
    adata.uns["iroot"] = 0
    return adata


@filter_oldformatwarning
def pbmc68k_reduced() -> AnnData:
    """\
    Subsampled and processed 68k PBMCs.

    `PBMC 68k dataset`_ from 10x Genomics.

    The original PBMC 68k dataset was preprocessed using scanpy and was saved
    keeping only 724 cells and 221 highly variable genes.

    The saved file contains the annotation of cell types (key: `'bulk_labels'`),
    UMAP coordinates, louvain clustering and gene rankings based on the
    `bulk_labels`.

    .. _PBMC 68k dataset: https://www.10xgenomics.com/datasets/fresh-68-k-pbm-cs-donor-a-1-standard-1-1-0

    Returns
    -------
    Annotated data matrix.
    """

    filename = HERE / "10x_pbmc68k_reduced.h5ad"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")
        return read(filename)


@filter_oldformatwarning
@check_datasetdir_exists
def pbmc3k() -> AnnData:
    r"""\
    3k PBMCs from 10x Genomics.

    The data consists in 3k PBMCs from a Healthy Donor and is freely available
    from 10x Genomics (file_ from this webpage_).

    The exact same data is also used in Seurat’s `basic clustering tutorial`_.

    .. _file: https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
    .. _webpage: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
    .. _basic clustering tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

    .. note::
       This downloads 5.9 MB of data upon the first call of the function and stores it in
       :attr:`~scanpy._settings.ScanpyConfig.datasetdir`\\ `/pbmc3k_raw.h5ad`.

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
    url = "https://falexwolf.de/data/pbmc3k_raw.h5ad"
    adata = read(settings.datasetdir / "pbmc3k_raw.h5ad", backup_url=url)
    return adata


@filter_oldformatwarning
@check_datasetdir_exists
def pbmc3k_processed() -> AnnData:
    """\
    Processed 3k PBMCs from 10x Genomics.

    Processed using the basic tutorial :doc:`/tutorials/basics/clustering-2017`.

    Returns
    -------
    Annotated data matrix.
    """
    url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")
        return read(settings.datasetdir / "pbmc3k_processed.h5ad", backup_url=url)


VisiumSampleID = Literal[
    "V1_Breast_Cancer_Block_A_Section_1",
    "V1_Breast_Cancer_Block_A_Section_2",
    "V1_Human_Heart",
    "V1_Human_Lymph_Node",
    "V1_Mouse_Kidney",
    "V1_Adult_Mouse_Brain",
    "V1_Mouse_Brain_Sagittal_Posterior",
    "V1_Mouse_Brain_Sagittal_Posterior_Section_2",
    "V1_Mouse_Brain_Sagittal_Anterior",
    "V1_Mouse_Brain_Sagittal_Anterior_Section_2",
    "V1_Human_Brain_Section_1",
    "V1_Human_Brain_Section_2",
    "V1_Adult_Mouse_Brain_Coronal_Section_1",
    "V1_Adult_Mouse_Brain_Coronal_Section_2",
    # spaceranger version 1.2.0
    "Targeted_Visium_Human_Cerebellum_Neuroscience",
    "Parent_Visium_Human_Cerebellum",
    "Targeted_Visium_Human_SpinalCord_Neuroscience",
    "Parent_Visium_Human_SpinalCord",
    "Targeted_Visium_Human_Glioblastoma_Pan_Cancer",
    "Parent_Visium_Human_Glioblastoma",
    "Targeted_Visium_Human_BreastCancer_Immunology",
    "Parent_Visium_Human_BreastCancer",
    "Targeted_Visium_Human_OvarianCancer_Pan_Cancer",
    "Targeted_Visium_Human_OvarianCancer_Immunology",
    "Parent_Visium_Human_OvarianCancer",
    "Targeted_Visium_Human_ColorectalCancer_GeneSignature",
    "Parent_Visium_Human_ColorectalCancer",
]


def _download_visium_dataset(
    sample_id: VisiumSampleID,
    spaceranger_version: Literal["1.1.0", "1.2.0"],
    *,
    base_dir: Path | None = None,
    download_image: bool = False,
) -> Path:
    """\
    Download Visium spatial data from 10x Genomics’ database.

    Params
    ------
    sample_id
        String name of example visium dataset.
    base_dir
        Where to download the dataset to.
    download_image
        Whether to download the high-resolution tissue section.
    """
    import tarfile

    if base_dir is None:
        base_dir = settings.datasetdir

    url_prefix = f"https://cf.10xgenomics.com/samples/spatial-exp/{spaceranger_version}/{sample_id}"

    sample_dir = base_dir / sample_id
    sample_dir.mkdir(exist_ok=True)

    # Download spatial data
    tar_filename = f"{sample_id}_spatial.tar.gz"
    tar_pth = sample_dir / tar_filename
    _utils.check_presence_download(
        filename=tar_pth, backup_url=f"{url_prefix}/{tar_filename}"
    )
    with tarfile.open(tar_pth) as f:
        f.extraction_filter = tarfile.data_filter
        for el in f:
            if not (sample_dir / el.name).exists():
                f.extract(el, sample_dir)

    # Download counts
    _utils.check_presence_download(
        filename=sample_dir / "filtered_feature_bc_matrix.h5",
        backup_url=f"{url_prefix}/{sample_id}_filtered_feature_bc_matrix.h5",
    )

    # Download image
    if download_image:
        _utils.check_presence_download(
            filename=sample_dir / "image.tif",
            backup_url=f"{url_prefix}/{sample_id}_image.tif",
        )

    return sample_dir


@check_datasetdir_exists
def visium_sge(
    sample_id: VisiumSampleID = "V1_Breast_Cancer_Block_A_Section_1",
    *,
    include_hires_tiff: bool = False,
) -> AnnData:
    """\
    Processed Visium Spatial Gene Expression data from 10x Genomics’ database_.

    .. _database: https://support.10xgenomics.com/spatial-gene-expression/datasets

    Parameters
    ----------
    sample_id
        The ID of the data sample in 10x’s spatial database.
    include_hires_tiff
        Download and include the high-resolution tissue image (tiff) in
        `adata.uns["spatial"][sample_id]["metadata"]["source_image_path"]`.

    Returns
    -------
    Annotated data matrix.
    """
    spaceranger_version = "1.1.0" if "V1_" in sample_id else "1.2.0"
    sample_dir = _download_visium_dataset(
        sample_id, spaceranger_version, download_image=include_hires_tiff
    )
    source_image_path = sample_dir / "image.tif" if include_hires_tiff else None
    return read_visium(sample_dir, source_image_path=source_image_path)
