from urllib.request import urlretrieve, urlopen
from urllib.error import HTTPError
from zipfile import ZipFile

from scipy import sparse
import pandas as pd
import numpy as np
from tqdm import tqdm

import anndata
from .._settings import settings


def _filter_boring(dataframe):
    unique_vals = dataframe.apply(lambda x: len(x.unique()))
    is_boring = (unique_vals == 1) | (unique_vals == len(dataframe))
    return dataframe.loc[:, ~is_boring]


# Copied from tqdm examples
def tqdm_hook(t):
    """
    Wraps tqdm instance.

    Don't forget to close() or __exit__()
    the tqdm instance once you're done with it (easiest using `with` syntax).
    Example
    -------
    >>> with tqdm(...) as t:
    ...     reporthook = my_hook(t)
    ...     urllib.urlretrieve(..., reporthook=reporthook)
    """
    last_b = [0]

    def update_to(b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b

    return update_to


def sniff_url(accession):
    # Note that data is downloaded from gxa/sc/experiment, not experiments
    base_url = "https://www.ebi.ac.uk/gxa/sc/experiments/{}/".format(accession)
    try:
        with urlopen(base_url) as req:  # Check if server up/ dataset exists
            pass
    except HTTPError as e:
        e.msg = e.msg + " ({})".format(base_url)  # Report failed url
        raise


def download_experiment(accession):
    sniff_url(accession)

    base_url = "https://www.ebi.ac.uk/gxa/sc/experiment/{}/".format(accession)
    quantification_path = "download/zip?fileType=quantification-filtered&accessKey="
    sampledata_path = "download?fileType=experiment-design&accessKey="

    experiment_dir = settings.datasetdir / accession
    experiment_dir.mkdir(parents=True, exist_ok=True)

    with tqdm(
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        miniters=1,
        desc="experimental_design.tsv",
    ) as t:
        urlretrieve(
            base_url + sampledata_path,
            experiment_dir / "experimental_design.tsv",
            reporthook=tqdm_hook(t),
        )
    with tqdm(
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        miniters=1,
        desc="expression_archive.zip",
    ) as t:
        urlretrieve(
            base_url + quantification_path,
            experiment_dir / "expression_archive.zip",
            reporthook=tqdm_hook(t),
        )


def read_mtx_from_stream(stream):
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


def read_expression_from_archive(archive: ZipFile):
    info = archive.infolist()
    assert len(info) == 3
    mtx_data_info = [i for i in info if i.filename.endswith(".mtx")][0]
    mtx_rows_info = [i for i in info if i.filename.endswith(".mtx_rows")][0]
    mtx_cols_info = [i for i in info if i.filename.endswith(".mtx_cols")][0]
    with archive.open(mtx_data_info, "r") as f:
        expr = read_mtx_from_stream(f)
    with archive.open(mtx_rows_info, "r") as f:
        varname = pd.read_csv(f, sep="\t", header=None)[1]
    with archive.open(mtx_cols_info, "r") as f:
        obsname = pd.read_csv(f, sep="\t", header=None)[1]
    adata = anndata.AnnData(expr)
    adata.var_names = varname
    adata.obs_names = obsname
    return adata


def ebi_expression_atlas(accession: str, *, filter_boring: bool = False):
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
    dataset_path = experiment_dir / "{}.h5ad".format(accession)
    try:
        adata = anndata.read(dataset_path)
        if filter_boring:
            adata.obs = _filter_boring(adata.obs)
        return adata
    except OSError:
        # Dataset couldn't be read for whatever reason
        pass

    download_experiment(accession)

    print("Downloaded {} to {}".format(accession, experiment_dir.absolute()))

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
