"""Extension to warn about numpydoc-style parameter types in docstrings."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.ext import autosummary

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from sphinx.application import Sphinx


import_by_name = autosummary.import_by_name


def _patched_import_by_name(
    name: str,
    prefixes: Sequence[str | None] = (None,),
) -> tuple[str, Any, Any, str]:
    if name.startswith("scanpy.settings"):
        prefixed_name, obj, parent, modname = import_by_name(
            name.replace("scanpy.settings", "scanpy._settings.Settings"), prefixes
        )
        prefixed_name = prefixed_name.replace(
            "scanpy._settings.Settings", "scanpy.settings"
        )
        if parent.__name__ == "scanpy._settings":
            parent = import_by_name("scanpy")[1]

        return prefixed_name, obj, parent, "scanpy"

    return import_by_name(name, prefixes)


def setup(app: Sphinx) -> None:
    """App setup hook."""
    autosummary.import_by_name = _patched_import_by_name
