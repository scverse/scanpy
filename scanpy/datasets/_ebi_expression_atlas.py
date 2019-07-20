from urllib.request import urlopen
from urllib.error import HTTPError
from zipfile import ZipFile
from typing import BinaryIO

import anndata
import pandas as pd
import numpy as np
from scipy import sparse

from ..readwrite import download
from .._settings import settings
from .. import logging as logg


def _filter_boring(dataframe: pd.DataFrame) -> pd.DataFrame:
    unique_vals = dataframe.apply(lambda x: len(x.unique()))
    is_boring = (unique_vals == 1) | (unique_vals == len(dataframe))
    return dataframe.loc[:, ~is_boring]


def sniff_url(accession: str):
    # Note that data is downloaded from gxa/sc/experiment, not experiments
    base_url = f"https://www.ebi.ac.uk/gxa/sc/experiments/{accession}/"
    try:
        with urlopen(base_url):  # Check if server up/ dataset exists
            pass
    except HTTPError as e:
        e.msg = e.msg + f" ({base_url})"  # Report failed url
        raise


def download_experiment(accession: str):
    sniff_url(accession)

    base_url = f"https://www.ebi.ac.uk/gxa/sc/experiment/{accession}/"
    quantification_path = "download/zip?fileType=quantification-filtered&accessKey="
    sampledata_path = "download?fileType=experiment-design&accessKey="

    experiment_dir = settings.datasetdir / accession
    experiment_dir.mkdir(parents=True, exist_ok=True)

    download(
        base_url + sampledata_path,
        experiment_dir / "experimental_design.tsv",
    )
    download(
        base_url + quantification_path,
        experiment_dir / "expression_archive.zip",
    )


def read_mtx_from_stream(stream: BinaryIO) -> sparse.csr_matrix:
    stream.readline()
    n, m, _ = (int(x) for x in stream.readline()[:-1].split(b" "))
    data = pd.read_csv(
        stream,
        sep=r"\s+",
        header=None,
        dtype={0: np.integer, 1: np.integer, 2: np.float32},
    )
    mtx = sparse.csr_matrix((data[2], (data[1] - 1, data[0] - 1)), shape=(m, n))
    return mtx


def read_expression_from_archive(archive: ZipFile) -> anndata.AnnData:
    info = archive.infolist()
    assert len(info) == 3
    mtx_data_info = next(i for i in info if i.filename.endswith(".mtx"))
    mtx_rows_info = next(i for i in info if i.filename.endswith(".mtx_rows"))
    mtx_cols_info = next(i for i in info if i.filename.endswith(".mtx_cols"))
    with archive.open(mtx_data_info, "r") as f:
        expr = read_mtx_from_stream(f)
    with archive.open(mtx_rows_info, "r") as f:
        varname = pd.read_csv(f, sep="\t", header=None)[1]  # TODO: Check what other value could be
    with archive.open(mtx_cols_info, "r") as f:
        obsname = pd.read_csv(f, sep="\t", header=None).iloc[:, 0]
    adata = anndata.AnnData(expr)
    adata.var_names = varname
    adata.obs_names = obsname
    return adata


def ebi_expression_atlas(accession: str, *, filter_boring: bool = False) -> anndata.AnnData:
    """Load a dataset from the `EBI Single Cell Expression Atlas <https://www.ebi.ac.uk/gxa/sc/experiments>`__.

    Downloaded datasets are saved in directory specified by `sc.settings.datasetdir`.

    Params
    ------
    accession
        Dataset accession. Like ``E-GEOD-98816`` or ``E-MTAB-4888``. This can
        be found in the url on the datasets page. For example:
        ``https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-98816/results/tsne``
    filter_boring
        Whether boring labels in `.obs` should be automatically removed.

    Example
    -------
    >>> adata = sc.datasets.ebi_expression_atlas("E-MTAB-4888")
    """
    experiment_dir = settings.datasetdir / accession
    dataset_path = experiment_dir / f"{accession}.h5ad"
    try:
        adata = anndata.read(dataset_path)
        if filter_boring:
            adata.obs = _filter_boring(adata.obs)
        return adata
    except OSError:
        # Dataset couldn't be read for whatever reason
        pass

    download_experiment(accession)

    logg.info(f"Downloaded {accession} to {experiment_dir.absolute()}")

    with ZipFile(experiment_dir / "expression_archive.zip", "r") as f:
        adata = read_expression_from_archive(f)
    obs = pd.read_csv(
        experiment_dir / "experimental_design.tsv", sep="\t", index_col=0
    )

    adata.obs[obs.columns] = obs
    adata.write(dataset_path, compression="gzip")  # To be kind to disk space

    if filter_boring:
        adata.obs = _filter_boring(adata.obs)

    return adata
