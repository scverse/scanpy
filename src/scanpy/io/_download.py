"""Helpers for reading and writing."""

from __future__ import annotations

from pathlib import Path, PurePath

from .. import logging as logg

__all__ = ["check_datafile_present_and_download", "download", "slugify"]


def slugify(path: str | PurePath) -> str:
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


def download(url: str, path: Path) -> None:
    from ssl import create_default_context
    from tempfile import NamedTemporaryFile
    from urllib.request import Request, urlopen

    from certifi import contents
    from tqdm.auto import tqdm

    blocksize = 1024 * 8
    blocknum = 0

    # Write to a temp file and rename so readers never see a partial file (#4097).
    tmp_path: Path | None = None
    try:
        req = Request(url, headers={"User-agent": "scanpy-user"})

        with urlopen(req, context=create_default_context(cadata=contents())) as resp:
            total = resp.info().get("content-length", None)
            with (
                tqdm(
                    unit="B",
                    unit_scale=True,
                    miniters=1,
                    unit_divisor=1024,
                    total=total if total is None else int(total),
                ) as t,
                NamedTemporaryFile(
                    dir=path.parent,
                    prefix=f"{path.name}.",
                    suffix=".part",
                    delete=False,
                ) as f,
            ):
                tmp_path = Path(f.name)
                block = resp.read(blocksize)
                while block:
                    f.write(block)
                    blocknum += 1
                    t.update(len(block))
                    block = resp.read(blocksize)

        tmp_path.replace(path)
        tmp_path = None

    except (KeyboardInterrupt, Exception):
        # Only remove our own temp file; leave path, which may be another process's.
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)
        raise


def check_datafile_present_and_download(
    path: Path, backup_url: str | None = None
) -> bool:
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

    download(backup_url, path)
    return True
