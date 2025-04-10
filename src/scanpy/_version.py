"""Get version from VCS in a dev environment or from package metadata in production.

See <https://github.com/maresb/hatch-vcs-footgun-example>.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["__version__"]


def _get_version_from_vcs() -> str:  # pragma: no cover
    from hatchling.metadata.core import ProjectMetadata
    from hatchling.plugin.exceptions import UnknownPluginError
    from hatchling.plugin.manager import PluginManager
    from hatchling.utils.fs import locate_file

    if (pyproject_toml := locate_file(__file__, "pyproject.toml")) is None:
        msg = "pyproject.toml not found although hatchling is installed"
        raise LookupError(msg)
    root = Path(pyproject_toml).parent
    metadata = ProjectMetadata(root=str(root), plugin_manager=PluginManager())
    try:
        # Version can be either statically set in pyproject.toml or computed dynamically:
        return metadata.core.version or metadata.hatch.version.cached
    except UnknownPluginError as e:
        msg = "Unable to import hatch plugin."
        raise ImportError(msg) from e


try:
    __version__ = _get_version_from_vcs()
except (ImportError, LookupError):
    import importlib.metadata

    __version__ = importlib.metadata.version("scanpy")
