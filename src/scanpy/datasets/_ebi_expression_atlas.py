from __future__ import annotations

from typing import TYPE_CHECKING
from urllib.error import HTTPError
from urllib.request import urlopen
from zipfile import ZipFile

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

from .. import logging as logg
from .._compat import add_note
from .._settings import settings
from .._utils._doctests import doctest_internet
from ..readwrite import _download
from ._utils import check_datasetdir_exists

if TYPE_CHECKING:
    from pandas._typing import ReadCsvBuffer

    from scanpy._compat import CSRBase


CHUNK_SIZE = int(1e7)


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
        add_note(e, base_url)
        raise


@check_datasetdir_exists
def download_experiment(accession: str):
    sniff_url(accession)

    base_url = f"https://www.ebi.ac.uk/gxa/sc/experiment/{accession}"
    design_url = f"{base_url}/download?accessKey=&fileType="
    mtx_url = f"{base_url}/download/zip?accessKey=&fileType="

    experiment_dir = settings.datasetdir / accession
    experiment_dir.mkdir(parents=True, exist_ok=True)

    _download(
        design_url + "experiment-design",
        experiment_dir / "experimental_design.tsv",
    )
    _download(
        mtx_url + "quantification-raw",
        experiment_dir / "expression_archive.zip",
    )


def read_mtx_from_stream(stream: ReadCsvBuffer[bytes]) -> CSRBase:
    curline = stream.readline()
    while curline.startswith(b"%"):
        curline = stream.readline()
    n, m, e = map(int, curline[:-1].split(b" "))

    dtype_data = np.float32
    max_int32 = np.iinfo(np.int32).max
    dtype_coord = np.int64 if n > max_int32 or m > max_int32 else np.int32

    data = np.ndarray((e,), dtype=dtype_data)
    i = np.ndarray((e,), dtype=dtype_coord)
    j = np.ndarray((e,), dtype=dtype_coord)
    start = 0
    with pd.read_csv(
        stream,
        sep=r"\s+",
        header=None,
        dtype={0: dtype_coord, 1: dtype_coord, 2: dtype_data},
        chunksize=CHUNK_SIZE,
    ) as reader:
        chunk: pd.DataFrame
        for chunk in reader:
            l = chunk.shape[0]
            data[start : start + l] = chunk[2]
            i[start : start + l] = chunk[1] - 1
            j[start : start + l] = chunk[0] - 1
            start += l
    return sparse.csr_matrix((data, (i, j)), shape=(m, n))  # noqa: TID251


def read_expression_from_archive(archive: ZipFile) -> anndata.AnnData:
    info = archive.infolist()
    assert len(info) == 3
    mtx_data_info = next(i for i in info if i.filename.endswith(".mtx"))
    mtx_rows_info = next(i for i in info if i.filename.endswith(".mtx_rows"))
    mtx_cols_info = next(i for i in info if i.filename.endswith(".mtx_cols"))
    with archive.open(mtx_data_info, "r") as f:
        expr = read_mtx_from_stream(f)
    with archive.open(mtx_rows_info, "r") as f:
        # TODO: Check what other value could be
        varname = pd.read_csv(f, sep="\t", header=None)[1]
    with archive.open(mtx_cols_info, "r") as f:
        obsname = pd.read_csv(f, sep="\t", header=None).iloc[:, 0]
    adata = anndata.AnnData(expr)
    adata.var_names = varname
    adata.obs_names = obsname
    return adata


@doctest_internet
def ebi_expression_atlas(
    accession: str, *, filter_boring: bool = False
) -> anndata.AnnData:
    """Load a dataset from the EBI Single Cell Expression Atlas.

    The atlas_ can be browsed online to find the ``accession`` you want.
    Downloaded datasets are saved in the directory specified by
    :attr:`~scanpy.settings.datasetdir`.

    .. _atlas: https://www.ebi.ac.uk/gxa/sc/experiments

    Params
    ------
    accession
        Dataset accession. Like ``E-GEOD-98816`` or ``E-MTAB-4888``.
        This can be found in the url on the datasets page, for example E-GEOD-98816_.

        .. _E-GEOD-98816: https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-98816/results/tsne
    filter_boring
        Whether boring labels in `.obs` should be automatically removed, such as
        labels with a single or :attr:`~anndata.AnnData.n_obs` distinct values.

    Returns
    -------
    Annotated data matrix.

    Example
    -------
    >>> import scanpy as sc
    >>> sc.datasets.ebi_expression_atlas("E-MTAB-4888")  # doctest: +ELLIPSIS
    AnnData object with n_obs × n_vars = 2261 × 23899
        obs: 'Sample Characteristic[organism]', 'Sample Characteristic Ontology Term[organism]', ..., 'Factor Value[cell type]', 'Factor Value Ontology Term[cell type]'

    """
    experiment_dir = settings.datasetdir / accession
    dataset_path = experiment_dir / f"{accession}.h5ad"
    try:
        adata = anndata.read_h5ad(dataset_path)
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
    obs = pd.read_csv(experiment_dir / "experimental_design.tsv", sep="\t", index_col=0)

    adata.obs[obs.columns] = obs
    adata.write(dataset_path, compression="gzip")  # To be kind to disk space

    if filter_boring:
        adata.obs = _filter_boring(adata.obs)

    return adata
