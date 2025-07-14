"""Reading and Writing."""

from __future__ import annotations

import json
from functools import partial
from pathlib import Path, PurePath
from typing import TYPE_CHECKING, overload

import anndata.utils
import h5py
import numpy as np
import pandas as pd
from packaging.version import Version

if Version(anndata.__version__) >= Version("0.11.0rc2"):
    from anndata.io import (
        read_csv,
        read_excel,
        read_h5ad,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
    )
else:
    from anndata import (
        read_csv,
        read_excel,
        read_h5ad,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
    )

from anndata import AnnData
from matplotlib.image import imread

from . import logging as logg
from ._compat import add_note, deprecated, old_positionals
from ._settings import settings
from ._utils import _empty

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import BinaryIO, Literal

    from ._utils import Empty

# .gz and .bz2 suffixes are also allowed for text formats
text_exts = {
    "csv",
    "tsv",
    "tab",
    "data",
    "txt",  # these four are all equivalent
}
avail_exts = {
    "anndata",
    "xlsx",
    "h5",
    "h5ad",
    "mtx",
    "mtx.gz",
    "soft.gz",
    "loom",
} | text_exts
"""Available file formats for reading data. """


# --------------------------------------------------------------------------------
# Reading and Writing data files and AnnData objects
# --------------------------------------------------------------------------------


@old_positionals(
    "sheet",
    "ext",
    "delimiter",
    "first_column_names",
    "backup_url",
    "cache",
    "cache_compression",
)
def read(
    filename: Path | str,
    backed: Literal["r", "r+"] | None = None,
    *,
    sheet: str | None = None,
    ext: str | None = None,
    delimiter: str | None = None,
    first_column_names: bool = False,
    backup_url: str | None = None,
    cache: bool = False,
    cache_compression: Literal["gzip", "lzf"] | None | Empty = _empty,
    **kwargs,
) -> AnnData:
    """Read file and return :class:`~anndata.AnnData` object.

    To speed up reading, consider passing ``cache=True``, which creates an hdf5
    cache file.

    Parameters
    ----------
    filename
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via ``sc.settings.writedir / (filename +
        sc.settings.file_format_data)``.  This is the same behavior as in
        ``sc.read(filename, ...)``.
    backed
        If ``'r'``, load :class:`~anndata.AnnData` in ``backed`` mode instead
        of fully loading it into memory (`memory` mode). If you want to modify
        backed attributes of the AnnData object, you need to choose ``'r+'``.
    sheet
        Name of sheet/table in hdf5 or Excel file.
    ext
        Extension that indicates the file type. If ``None``, uses extension of
        filename.
    delimiter
        Delimiter that separates data within text file. If ``None``, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at any single white space ``' '``.
    first_column_names
        Assume the first column stores row names. This is only necessary if
        these are not strings: strings in the first column are automatically
        assumed to be row names.
    backup_url
        Retrieve the file from an URL if not present on disk.
    cache
        If `False`, read from source, if `True`, read from fast 'h5ad' cache.
    cache_compression
        See the h5py :ref:`dataset_compression`.
        (Default: `settings.cache_compression`)
    kwargs
        Parameters passed to :func:`~anndata.io.read_loom`.

    Returns
    -------
    An :class:`~anndata.AnnData` object

    """
    filename = Path(filename)  # allow passing strings
    if is_valid_filename(filename, ext=ext):
        return _read(
            filename,
            backed=backed,
            sheet=sheet,
            ext=ext,
            delimiter=delimiter,
            first_column_names=first_column_names,
            backup_url=backup_url,
            cache=cache,
            cache_compression=cache_compression,
            **kwargs,
        )
    # generate filename and read to dict
    filekey = str(filename)
    filename = settings.writedir / (filekey + "." + settings.file_format_data)
    if not filename.exists():
        msg = (
            f"Reading with filekey {filekey!r} failed, "
            f"the inferred filename {filename!r} does not exist. "
            "If you intended to provide a filename, either use a filename "
            f"ending on one of the available extensions {avail_exts} "
            "or pass the parameter `ext`."
        )
        raise ValueError(msg)
    return read_h5ad(filename, backed=backed)


@old_positionals("genome", "gex_only", "backup_url")
def read_10x_h5(
    filename: Path | str,
    *,
    genome: str | None = None,
    gex_only: bool = True,
    backup_url: str | None = None,
) -> AnnData:
    r"""Read 10x-Genomics-formatted hdf5 file.

    Parameters
    ----------
    filename
        Path to a 10x hdf5 file.
    genome
        Filter expression to genes within this genome. For legacy 10x h5
        files, this must be provided if the data contains more than one genome.
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    backup_url
        Retrieve the file from an URL if not present on disk.

    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:

    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names for a feature barcode matrix, probe names for a probe bc matrix
    :attr:`~anndata.AnnData.var`\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.obs`\ `[filtered_barcodes]`
        filtered barcodes if present in the matrix
    :attr:`~anndata.AnnData.var`
        Any additional metadata present in /matrix/features is read in.

    """
    path = Path(filename)
    start = logg.info(f"reading {path}")
    is_present = _check_datafile_present_and_download(path, backup_url=backup_url)
    if not is_present:
        logg.debug(f"... did not find original file {path}")
    with h5py.File(str(path), "r") as f:
        v3 = "/matrix" in f
    if v3:
        adata = _read_10x_h5(path, _read_v3_10x_h5)
        if genome:
            if genome not in adata.var["genome"].values:
                msg = (
                    f"Could not find data corresponding to genome {genome!r} in {path}. "
                    f"Available genomes are: {list(adata.var['genome'].unique())}."
                )
                raise ValueError(msg)
            adata = adata[:, adata.var["genome"] == genome]
        if gex_only:
            adata = adata[:, adata.var["feature_types"] == "Gene Expression"]
        if adata.is_view:
            adata = adata.copy()
    else:
        adata = _read_10x_h5(path, partial(_read_legacy_10x_h5, genome=genome))
    logg.info("", time=start)
    return adata


def _read_10x_h5(path: Path, cb: Callable[[h5py.File], AnnData]) -> AnnData:
    """Read hdf5 file from Cell Ranger v3 or later versions."""
    with h5py.File(str(path), "r") as f:
        try:
            return cb(f)
        except KeyError as e:
            msg = "File is missing one or more required datasets."
            raise Exception(msg) from e


def _collect_datasets(dsets: dict, group: h5py.Group) -> None:
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            dsets[k] = v[()]
        else:
            _collect_datasets(dsets, v)


def _read_v3_10x_h5(f: h5py.File) -> AnnData:
    dsets = {}
    _collect_datasets(dsets, f["matrix"])

    from scipy.sparse import csr_matrix  # noqa: TID251

    M, N = dsets["shape"]
    data = dsets["data"]
    if dsets["data"].dtype == np.dtype("int32"):
        data = dsets["data"].view("float32")
        data[:] = dsets["data"]
    matrix = csr_matrix(
        (data, dsets["indices"], dsets["indptr"]),
        shape=(N, M),
    )
    obs_dict = {"obs_names": dsets["barcodes"].astype(str)}
    var_dict = {"var_names": dsets["name"].astype(str)}

    if "gene_id" not in dsets:
        # Read metadata specific to a feature-barcode matrix
        var_dict["gene_ids"] = dsets["id"].astype(str)
    else:
        # Read metadata specific to a probe-barcode matrix
        var_dict.update(
            {
                "gene_ids": dsets["gene_id"].astype(str),
                "probe_ids": dsets["id"].astype(str),
            }
        )
    var_dict["feature_types"] = dsets["feature_type"].astype(str)
    if "filtered_barcodes" in f["matrix"]:
        obs_dict["filtered_barcodes"] = dsets["filtered_barcodes"].astype(bool)

    if "features" in f["matrix"]:
        var_dict.update(
            (
                feature_metadata_name,
                dsets[feature_metadata_name].astype(
                    bool if feature_metadata_item.dtype.kind == "b" else str
                ),
            )
            for feature_metadata_name, feature_metadata_item in f["matrix"][
                "features"
            ].items()
            if isinstance(feature_metadata_item, h5py.Dataset)
            and feature_metadata_name
            not in ["name", "feature_type", "id", "gene_id", "_all_tag_keys"]
        )
    else:
        msg = "10x h5 has no features group"
        raise ValueError(msg)
    return AnnData(matrix, obs=obs_dict, var=var_dict)


def _read_legacy_10x_h5(f: h5py.File, genome: str | None) -> AnnData:
    children = list(f.keys())
    if not genome:
        if len(children) > 1:
            msg = (
                f"{f.filename} contains more than one genome. "
                "For legacy 10x h5 files you must specify the genome "
                "if more than one is present. "
                f"Available genomes are: {children}"
            )
            raise ValueError(msg)
        genome = children[0]
    elif genome not in children:
        msg = (
            f"Could not find genome {genome!r} in {f.filename}. "
            f"Available genomes are: {children}"
        )
        raise ValueError(msg)

    dsets = {}
    _collect_datasets(dsets, f[genome])

    # AnnData works with csr matrices
    # 10x stores the transposed data, so we do the transposition right away
    from scipy.sparse import csr_matrix  # noqa: TID251

    M, N = dsets["shape"]
    data = dsets["data"]
    if dsets["data"].dtype == np.dtype("int32"):
        data = dsets["data"].view("float32")
        data[:] = dsets["data"]
    matrix = csr_matrix(
        (data, dsets["indices"], dsets["indptr"]),
        shape=(N, M),
    )
    # the csc matrix is automatically the transposed csr matrix
    # as scanpy expects it, so, no need for a further transpostion
    adata = AnnData(
        matrix,
        obs=dict(obs_names=dsets["barcodes"].astype(str)),
        var=dict(
            var_names=dsets["gene_names"].astype(str),
            gene_ids=dsets["genes"].astype(str),
        ),
    )
    return adata


@deprecated("Use `squidpy.read.visium` instead.")
def read_visium(
    path: Path | str,
    genome: str | None = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str | None = None,
    load_images: bool | None = True,
    source_image_path: Path | str | None = None,
) -> AnnData:
    r"""Read 10x-Genomics-formatted visum dataset.

    .. deprecated:: 1.11.0
       Use :func:`squidpy.read.visium` instead.

    In addition to reading regular 10x output,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    Based on the `Space Ranger output docs`_.

    See :func:`~scanpy.pl.spatial` for a compatible plotting function.

    .. _Space Ranger output docs: <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview>

    Parameters
    ----------
    path
        Path to directory for visium datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.

    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:

    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names for a feature barcode matrix, probe names for a probe bc matrix
    :attr:`~anndata.AnnData.var`\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.obs`\ `[filtered_barcodes]`
        filtered barcodes if present in the matrix
    :attr:`~anndata.AnnData.var`
        Any additional metadata present in /matrix/features is read in.
    :attr:`~anndata.AnnData.uns`\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.

    """
    path = Path(path)
    adata = read_10x_h5(path / count_file, genome=genome)

    adata.uns["spatial"] = dict()

    from h5py import File

    with File(path / count_file, mode="r") as f:
        attrs = dict(f.attrs)
    if library_id is None:
        library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        tissue_positions_file = (
            path / "spatial/tissue_positions.csv"
            if (path / "spatial/tissue_positions.csv").exists()
            else path / "spatial/tissue_positions_list.csv"
        )
        files = dict(
            tissue_positions_file=tissue_positions_file,
            scalefactors_json_file=path / "spatial/scalefactors_json.json",
            hires_image=path / "spatial/tissue_hires_image.png",
            lowres_image=path / "spatial/tissue_lowres_image.png",
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    logg.warning(
                        f"You seem to be missing an image file.\nCould not find {f}."
                    )
                else:
                    msg = f"Could not find {f}"
                    raise OSError(msg)

        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["hires", "lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(
                    str(files[f"{res}_image"])
                )
            except Exception as e:  # noqa: PERF203
                msg = f"Could not find '{res}_image'"
                raise OSError(msg) from e

        # read json scalefactors
        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(
            files["scalefactors_json_file"].read_bytes()
        )

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        # read coordinates
        positions = pd.read_csv(
            files["tissue_positions_file"],
            header=0 if tissue_positions_file.name == "tissue_positions.csv" else None,
            index_col=0,
        )
        positions.columns = [
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm["spatial"] = adata.obs[
            ["pxl_row_in_fullres", "pxl_col_in_fullres"]
        ].to_numpy()
        adata.obs.drop(
            columns=["pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata


@old_positionals("var_names", "make_unique", "cache", "cache_compression", "gex_only")
def read_10x_mtx(
    path: Path | str,
    *,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = False,
    cache_compression: Literal["gzip", "lzf"] | None | Empty = _empty,
    gex_only: bool = True,
    prefix: str | None = None,
) -> AnnData:
    """Read 10x-Genomics-formatted mtx directory.

    Parameters
    ----------
    path
        Path to directory for `.mtx` and `.tsv` files,
        e.g. './filtered_gene_bc_matrices/hg19/'.
    var_names
        The variables index.
    make_unique
        Whether to make the variables index unique by appending '-1',
        '-2' etc. or not.
    cache
        If `False`, read from source, if `True`, read from fast 'h5ad' cache.
    cache_compression
        See the h5py :ref:`dataset_compression`.
        (Default: `settings.cache_compression`)
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    prefix
        Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
        if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
        `patientA_barcodes.tsv` the prefix is `patientA_`.
        (Default: no prefix)

    Returns
    -------
    An :class:`~anndata.AnnData` object

    """
    path = Path(path)
    prefix = "" if prefix is None else prefix
    is_legacy = (path / f"{prefix}genes.tsv").is_file()
    adata = _read_10x_mtx(
        path,
        var_names=var_names,
        make_unique=make_unique,
        cache=cache,
        cache_compression=cache_compression,
        prefix=prefix,
        is_legacy=is_legacy,
    )
    if is_legacy or not gex_only:
        return adata
    gex_rows = adata.var["feature_types"] == "Gene Expression"
    return adata[:, gex_rows].copy()


def _read_10x_mtx(
    path: Path,
    *,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = False,
    cache_compression: Literal["gzip", "lzf"] | None | Empty = _empty,
    prefix: str = "",
    is_legacy: bool,
) -> AnnData:
    """Read mex from output from Cell Ranger v2- or v3+."""
    suffix = "" if is_legacy else ".gz"
    adata = read(
        path / f"{prefix}matrix.mtx{suffix}",
        cache=cache,
        cache_compression=cache_compression,
    ).T  # transpose the data
    genes = pd.read_csv(
        path / f"{prefix}{'genes' if is_legacy else 'features'}.tsv{suffix}",
        header=None,
        sep="\t",
    )
    if var_names == "gene_symbols":
        var_names_idx = pd.Index(genes[1].values)
        if make_unique:
            var_names_idx = anndata.utils.make_index_unique(var_names_idx)
        adata.var_names = var_names_idx
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    else:
        msg = "`var_names` needs to be 'gene_symbols' or 'gene_ids'"
        raise ValueError(msg)
    if not is_legacy:
        adata.var["feature_types"] = genes[2].values
    barcodes = pd.read_csv(path / f"{prefix}barcodes.tsv{suffix}", header=None)
    adata.obs_names = barcodes[0].values
    return adata


@old_positionals("ext", "compression", "compression_opts")
def write(
    filename: Path | str,
    adata: AnnData,
    *,
    ext: Literal["h5", "csv", "txt", "npz"] | None = None,
    compression: Literal["gzip", "lzf"] | None = "gzip",
    compression_opts: int | None = None,
):
    """Write :class:`~anndata.AnnData` objects to file.

    Parameters
    ----------
    filename
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via `sc.settings.writedir / (filename +
        sc.settings.file_format_data)`. This is the same behavior as in
        :func:`~scanpy.read`.
    adata
        Annotated data matrix.
    ext
        File extension from wich to infer file format. If `None`, defaults to
        `sc.settings.file_format_data`.
    compression
        See https://docs.h5py.org/en/latest/high/dataset.html.
    compression_opts
        See https://docs.h5py.org/en/latest/high/dataset.html.

    """
    filename = Path(filename)  # allow passing strings
    if is_valid_filename(filename):
        ext_ = is_valid_filename(filename, return_ext=True)
        if ext is None:
            ext = ext_
        elif ext != ext_:
            msg = (
                "It suffices to provide the file type by "
                "providing a proper extension to the filename."
                'One of "txt", "csv", "h5" or "npz".'
            )
            raise ValueError(msg)
    else:
        key = filename
        ext = settings.file_format_data if ext is None else ext
        filename = _get_filename_from_key(key, ext)
    if ext == "csv":
        adata.write_csvs(filename)
    else:
        adata.write(
            filename, compression=compression, compression_opts=compression_opts
        )


# -------------------------------------------------------------------------------
# Reading and writing parameter files
# -------------------------------------------------------------------------------


@old_positionals("as_header")
def read_params(
    filename: Path | str, *, as_header: bool = False
) -> dict[str, int | float | bool | str | None]:
    """Read parameter dictionary from text file.

    Assumes that parameters are specified in the format::

        par1 = value1
        par2 = value2

    Comments that start with '#' are allowed.

    Parameters
    ----------
    filename
        Filename of data file.
    asheader
        Read the dictionary from the header (comment section) of a file.

    Returns
    -------
    Dictionary that stores parameters.

    """
    filename = Path(filename)  # allow passing str objects
    from collections import OrderedDict

    params = OrderedDict([])
    for line_raw in filename.open():
        if "=" in line_raw and (not as_header or line_raw.startswith("#")):
            line = line_raw[1:] if line_raw.startswith("#") else line_raw
            key, val = line.split("=")
            key = key.strip()
            val = val.strip()
            params[key] = convert_string(val)
    return params


def write_params(path: Path | str, *args, **maps):
    """Write parameters to file, so that it's readable by read_params.

    Uses INI file format.
    """
    path = Path(path)
    if not path.parent.is_dir():
        path.parent.mkdir(parents=True)
    if len(args) == 1:
        maps[None] = args[0]
    with path.open("w") as f:
        for header, map in maps.items():
            if header is not None:
                f.write(f"[{header}]\n")
            for key, val in map.items():
                f.write(f"{key} = {val}\n")


# -------------------------------------------------------------------------------
# Reading and Writing data files
# -------------------------------------------------------------------------------


def _read(  # noqa: PLR0912, PLR0915
    filename: Path,
    *,
    backed=None,
    sheet=None,
    ext=None,
    delimiter=None,
    first_column_names=None,
    backup_url=None,
    cache=False,
    cache_compression=None,
    suppress_cache_warning=False,
    **kwargs,
):
    if ext is not None and ext not in avail_exts:
        msg = f"Please provide one of the available extensions.\n{avail_exts}"
        raise ValueError(msg)
    else:
        ext = is_valid_filename(filename, return_ext=True, ext=ext)
    is_present = _check_datafile_present_and_download(filename, backup_url=backup_url)
    if not is_present:
        logg.debug(f"... did not find original file {filename}")
    # read hdf5 files
    if ext in {"h5", "h5ad"}:
        if sheet is None:
            return read_h5ad(filename, backed=backed)
        else:
            logg.debug(f"reading sheet {sheet} from file {filename}")
            return read_hdf(filename, sheet)
    # read other file types
    path_cache: Path = settings.cachedir / _slugify(filename).replace(
        f".{ext}", ".h5ad"
    )
    if path_cache.suffix in {".gz", ".bz2"}:
        path_cache = path_cache.with_suffix("")
    if cache and path_cache.is_file():
        logg.info(f"... reading from cache file {path_cache}")
        return read_h5ad(path_cache)

    if not is_present:
        msg = f"Did not find file {filename}."
        raise FileNotFoundError(msg)
    logg.debug(f"reading {filename}")
    if not cache and not suppress_cache_warning:
        logg.hint(
            "This might be very slow. Consider passing `cache=True`, "
            "which enables much faster reading from a cache file."
        )
    # do the actual reading
    if ext in {"xlsx", "xls"}:
        if sheet is None:
            msg = "Provide `sheet` parameter when reading '.xlsx' files."
            raise ValueError(msg)
        else:
            adata = read_excel(filename, sheet)
    elif ext in {"mtx", "mtx.gz"}:
        adata = read_mtx(filename)
    elif ext == "csv":
        if delimiter is None:
            delimiter = ","
        adata = read_csv(
            filename, first_column_names=first_column_names, delimiter=delimiter
        )
    elif ext in {"txt", "tab", "data", "tsv"}:
        if ext == "data":
            logg.hint(
                "... assuming '.data' means tab or white-space separated text file"
            )
            logg.hint("change this by passing `ext` to sc.read")
        adata = read_text(filename, delimiter, first_column_names=first_column_names)
    elif ext == "soft.gz":
        adata = _read_softgz(filename)
    elif ext == "loom":
        adata = read_loom(filename=filename, **kwargs)
    else:
        msg = f"Unknown extension {ext}."
        raise ValueError(msg)
    if cache:
        logg.info(
            f"... writing an {settings.file_format_data} "
            "cache file to speedup reading next time"
        )
        if cache_compression is _empty:
            cache_compression = settings.cache_compression
        if not path_cache.parent.is_dir():
            path_cache.parent.mkdir(parents=True)
        # write for faster reading when calling the next time
        adata.write(path_cache, compression=cache_compression)
    return adata


def _slugify(path: str | PurePath) -> str:
    """Make a path into a filename."""
    if not isinstance(path, PurePath):
        path = PurePath(path)
    parts = list(path.parts)
    if parts[0] == "/":
        parts.pop(0)
    elif len(parts[0]) == 3 and parts[0][1:] == ":\\":
        parts[0] = parts[0][0]  # C:\ → C
    filename = "-".join(parts)
    assert "/" not in filename, filename
    assert not filename[1:].startswith(":"), filename
    return filename


def _read_softgz(filename: str | bytes | Path | BinaryIO) -> AnnData:
    """Read a SOFT format data file.

    The SOFT format is documented here
    https://www.ncbi.nlm.nih.gov/geo/info/soft.html.

    Notes
    -----
    The function is based on a script by Kerby Shedden.
    https://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html

    """
    import gzip

    with gzip.open(filename, mode="rt") as file:
        # The header part of the file contains information about the
        # samples. Read that information first.
        samples_info = {}
        for line in file:
            if line.startswith("!dataset_table_begin"):
                break
            elif line.startswith("!subset_description"):
                subset_description = line.split("=")[1].strip()
            elif line.startswith("!subset_sample_id"):
                subset_ids = line.split("=")[1].split(",")
                subset_ids = [x.strip() for x in subset_ids]
                for k in subset_ids:
                    samples_info[k] = subset_description
        # Next line is the column headers (sample id's)
        sample_names = file.readline().strip().split("\t")
        # The column indices that contain gene expression data
        indices = [i for i, x in enumerate(sample_names) if x.startswith("GSM")]
        # Restrict the column headers to those that we keep
        sample_names = [sample_names[i] for i in indices]
        # Get a list of sample labels
        groups = [samples_info[k] for k in sample_names]
        # Read the gene expression data as a list of lists, also get the gene
        # identifiers
        gene_names, X = [], []
        for line in file:
            # This is what signals the end of the gene expression data
            # section in the file
            if line.startswith("!dataset_table_end"):
                break
            V = line.split("\t")
            # Extract the values that correspond to gene expression measures
            # and convert the strings to numbers
            x = [float(V[i]) for i in indices]
            X.append(x)
            gene_names.append(V[1])
    # Convert the Python list of lists to a Numpy array and transpose to match
    # the Scanpy convention of storing samples in rows and variables in colums.
    X = np.array(X).T
    obs = pd.DataFrame({"groups": groups}, index=sample_names)
    var = pd.DataFrame(index=gene_names)
    return AnnData(X=X, obs=obs, var=var)


# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


def is_float(string: str) -> float:
    """Check whether string is float.

    See Also
    --------
    https://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python

    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string: str) -> bool:
    """Check whether string is integer."""
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string: str) -> tuple[bool, bool]:
    """Check whether string is boolean."""
    if string == "True":
        return True, True
    elif string == "False":
        return True, False
    else:
        return False, False


def convert_string(string: str) -> int | float | bool | str | None:
    """Convert string to int, float or bool."""
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == "None":
        return None
    else:
        return string


# -------------------------------------------------------------------------------
# Helper functions for reading and writing
# -------------------------------------------------------------------------------


def get_used_files():
    """Get files used by processes with name scanpy."""
    import psutil

    loop_over_scanpy_processes = (
        proc for proc in psutil.process_iter() if proc.name() == "scanpy"
    )
    filenames = []
    for proc in loop_over_scanpy_processes:
        try:
            flist = proc.open_files()
            filenames.extend(nt.path for nt in flist)
        # This catches a race condition where a process ends
        # before we can examine its files
        except psutil.NoSuchProcess:  # noqa: PERF203
            pass
    return set(filenames)


def _get_filename_from_key(key, ext=None) -> Path:
    ext = settings.file_format_data if ext is None else ext
    return settings.writedir / f"{key}.{ext}"


def _download(url: str, path: Path):
    from urllib.error import URLError
    from urllib.request import Request, urlopen

    from tqdm.auto import tqdm

    blocksize = 1024 * 8
    blocknum = 0

    try:
        req = Request(url, headers={"User-agent": "scanpy-user"})

        try:
            open_url = urlopen(req)
        except URLError:
            if not url.startswith("https://"):
                raise  # No need to try using certifi

            msg = "Failed to open the url with default certificates."
            try:
                from certifi import where
            except ImportError as e:
                add_note(e, f"{msg} Please install `certifi` and try again.")
                raise
            else:
                logg.warning(f"{msg} Trying to use certifi.")

            from ssl import create_default_context

            open_url = urlopen(req, context=create_default_context(cafile=where()))

        with open_url as resp:
            total = resp.info().get("content-length", None)
            with (
                tqdm(
                    unit="B",
                    unit_scale=True,
                    miniters=1,
                    unit_divisor=1024,
                    total=total if total is None else int(total),
                ) as t,
                path.open("wb") as f,
            ):
                block = resp.read(blocksize)
                while block:
                    f.write(block)
                    blocknum += 1
                    t.update(len(block))
                    block = resp.read(blocksize)

    except (KeyboardInterrupt, Exception):
        # Make sure file doesn’t exist half-downloaded
        if path.is_file():
            path.unlink()
        raise


def _check_datafile_present_and_download(path: Path, backup_url=None):
    """Check whether the file is present, otherwise download."""
    path = Path(path)
    if path.is_file():
        return True
    if backup_url is None:
        return False
    logg.info(
        f"try downloading from url\n{backup_url}\n"
        "... this may take a while but only happens once"
    )
    if not path.parent.is_dir():
        logg.info(f"creating directory {path.parent}/ for saving data")
        path.parent.mkdir(parents=True)

    _download(backup_url, path)
    return True


@overload
def is_valid_filename(
    filename: Path, *, return_ext: Literal[False] = False, ext: str | None = None
) -> bool: ...
@overload
def is_valid_filename(
    filename: Path, *, return_ext: Literal[True], ext: str | None = None
) -> str: ...
def is_valid_filename(
    filename: Path, *, return_ext: bool = False, ext: str | None = None
) -> str | bool:
    """Check whether the argument is a filename."""
    ext_from_file = filename.suffixes
    if ext is not None:
        if not (joined_file_ext := ".".join(ext_from_file)).endswith(ext):
            msg = f"{joined_file_ext} does not end in expected extension {ext}"
            raise ValueError(msg)
        return ext if return_ext else True
    if len(ext_from_file) > 2:
        logg.warning(
            f"Your filename has more than two extensions: {ext_from_file}.\n"
            f"Only considering the two last: {ext_from_file[-2:]}."
        )
        ext_from_file = ext_from_file[-2:]

    # cases for gzipped/bzipped text files
    if (
        len(ext_from_file) == 2
        and ext_from_file[0][1:] in text_exts
        and ext_from_file[1][1:] in ("gz", "bz2")
    ):
        return ext_from_file[0][1:] if return_ext else True
    elif ext_from_file and ext_from_file[-1][1:] in avail_exts:
        return ext_from_file[-1][1:] if return_ext else True
    elif "".join(ext_from_file) == ".soft.gz":
        return "soft.gz" if return_ext else True
    elif "".join(ext_from_file) == ".mtx.gz":
        return "mtx.gz" if return_ext else True
    elif not return_ext:
        return False
    msg = f"""\
{filename!r} does not end on a valid extension.
Please, provide one of the available extensions.
{avail_exts}
Text files with .gz and .bz2 extensions are also supported.\
"""
    raise ValueError(msg)
