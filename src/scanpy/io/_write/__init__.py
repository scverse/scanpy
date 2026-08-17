"""Writing."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, cast, get_args

from anndata.io import write_h5ad, write_zarr

from ..._compat import warn
from ..._settings import AnnDataFileFormat, settings
from .exporting import write_cellbrowser, write_spring_project

if TYPE_CHECKING:
    from os import PathLike
    from typing import Literal

    from anndata import AnnData


__all__ = ["write", "write_cellbrowser", "write_spring_project"]


def write(
    filename: PathLike[str] | str,
    adata: AnnData,
    *,
    ext: AnnDataFileFormat | Literal["csv"] | None = None,
    convert_strings_to_categoricals: bool = True,
    compression: Literal["gzip", "lzf"] | None = "gzip",
    compression_opts: int | None = None,
) -> None:
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
        File extension from which to infer file format.
        If `None`, defaults to `sc.settings.file_format_data`.
    convert_strings_to_categoricals
        If anndata supports it, setting this to `False` will avoid
        converting string columns to categorical arrays when writing.
    compression
        See https://docs.h5py.org/en/latest/high/dataset.html.
    compression_opts
        See https://docs.h5py.org/en/latest/high/dataset.html.

    """
    filename = Path(filename)  # allow passing strings
    valid_exts = cast(
        "set[Literal['csv'] | AnnDataFileFormat]", {"csv", *get_args(AnnDataFileFormat)}
    )
    if filename.suffix and (ext_from_name := filename.suffix[1:]) in valid_exts:
        if ext is None:
            ext = ext_from_name
        elif ext != ext_from_name:
            msg = (
                "It suffices to provide the file type by "
                "providing a proper extension to the filename."
                f"One of {valid_exts}."
            )
            raise ValueError(msg)
    else:
        key = filename
        ext = settings.file_format_data if ext is None else ext
        filename = _get_filename_from_key(key, ext)

    if ext == "csv":
        msg = (
            "'csv' is not a good choice for anything, especially storing AnnData, "
            "and will be removed from this function. Use 'h5ad' or 'zarr' instead."
        )
        warn(msg, FutureWarning)
        adata.write_csvs(filename)
        return
    elif ext not in {"h5ad", "h5", "zarr"}:
        msg = f"Unknown file format: {ext} (not in {valid_exts})"
        raise ValueError(msg)

    if ext == "zarr":
        write_zarr(
            filename,
            adata,
            convert_strings_to_categoricals=convert_strings_to_categoricals,
        )
    else:
        write_h5ad(
            filename,
            adata,
            convert_strings_to_categoricals=convert_strings_to_categoricals,
            compression=compression,
            compression_opts=compression_opts,
        )


def _get_filename_from_key(key, ext=None) -> Path:
    ext = settings.file_format_data if ext is None else ext
    return settings.writedir / f"{key}.{ext}"
