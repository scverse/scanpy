from __future__ import annotations

from dataclasses import replace
from functools import cache, wraps
from importlib.resources import as_file, files
from pathlib import Path
from typing import TYPE_CHECKING

from scverse_misc.datasets import fetch, parse_registry, register_loader

from .. import logging as logg
from .._settings import settings

if TYPE_CHECKING:
    from collections.abc import Callable

    from scverse_misc.datasets import DatasetEntry, DownloadCB


def check_datasetdir_exists[**P, R](f: Callable[P, R]) -> Callable[P, R]:
    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        settings.datasetdir.mkdir(exist_ok=True)
        return f(*args, **kwargs)

    return wrapper


@register_loader("scanpy")
def _load_file(entry: DatasetEntry, target: Path, download: DownloadCB, /) -> Path:
    """Download the single file making up `entry` and return its path.

    Falls back to the entry’s ``fallback_urls`` (the pre-S3 upstream locations)
    if our bucket is unreachable or serves something that fails the hash check.
    """
    (file,) = entry.files
    # `target` is `datasetdir/scanpy`; its parent keeps scanpy’s flat layout.
    dest = target.parent
    try:
        return Path(download(file, dest=dest))
    except (OSError, ValueError) as e:
        if (url := entry.metadata.get("fallback_urls", {}).get(file.name)) is None:
            raise
        logg.warning(f"Downloading {file.name} failed ({e}), falling back to {url}")
        return Path(download(replace(file, url=url), dest=dest))


@cache
def _registry() -> tuple[str | None, dict[str, DatasetEntry]]:
    with as_file(files(__package__) / "registry.yaml") as path:
        return parse_registry(path)


def fetch_dataset(name: str) -> Path:
    """Download `name` into :attr:`~scanpy.settings.datasetdir` if needed, return its path."""
    base_url, datasets = _registry()
    return fetch(datasets[name], settings.datasetdir, base_url=base_url)
