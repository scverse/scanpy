from __future__ import annotations

import enum
import sys
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING, TypedDict, cast

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableSet
    from typing import NotRequired


class TunaColor(enum.IntEnum):
    Func = 0
    Builtin = 1
    Deprecated = 2


class TunaProf(TypedDict):
    text: list[str]
    value: float
    color: TunaColor
    children: NotRequired[list[TunaProf]]


def descend(
    profile: TunaProf, modules: MutableSet[str], path: Iterable[str] = ()
) -> Iterable[str]:
    [module] = profile["text"]
    path = [*path, module]
    if module in modules:
        yield " → ".join(e for e in path if e is not None)
        modules.remove(module)
    for child in profile.get("children", []):
        yield from descend(child, modules, path)


def get_import_paths(modules: Iterable[str]) -> Iterable[str]:
    from tuna import read_import_profile

    proc = run(
        [sys.executable, "-X", "importtime", "-c", "import scanpy"],
        capture_output=True,
        check=True,
    )
    with NamedTemporaryFile() as f:
        Path(f.name).write_bytes(proc.stderr)
        data = cast("TunaProf", read_import_profile(f.name))
    return descend(data, set(modules))


def test_deferred_imports(imported_modules: frozenset[str]) -> None:
    slow_to_import = {
        "umap",  # neighbors, tl.umap
        "seaborn",  # plotting
        "sklearn.metrics",  # neighbors
        "pynndescent",  # neighbors
        "networkx",  # diffmap, paga, plotting._utils
        # TODO: "matplotlib.pyplot",
        # TODO (maybe): "numba",
    }
    falsely_imported = slow_to_import & imported_modules

    assert not falsely_imported, "\n".join(get_import_paths(falsely_imported))
