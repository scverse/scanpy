"""Get version from VCS in a dev environment or from package metadata in production.

See <https://github.com/maresb/hatch-vcs-footgun-example>.
"""

from __future__ import annotations

import warnings
from pathlib import Path

__all__ = ["__version__"]

_PROJECT_NAME = "scanpy"


class GetVersionError(Exception):
    pass


def _get_version_from_vcs(project_name: str) -> str:  # pragma: no cover
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
        version = metadata.core.version or metadata.hatch.version.cached
    except UnknownPluginError as e:
        msg = "Unable to import hatch plugin."
        raise ImportError(msg) from e
    except ValueError as e:
        msg = f"Could not find hatchling project data in TOML file, {pyproject_toml}"
        raise GetVersionError(msg) from e
    except TypeError as e:
        msg = "Could not parse build configuration."
        raise GetVersionError(msg) from e
    except Exception as e:
        msg = (
            f"Unknown error getting version from hatchling config for '{project_name}'."
        )
        warnings.warn(f"{msg}: {e}", stacklevel=1)
        raise GetVersionError(msg) from e

    # We found a hatchling environment, but is it ours?
    if metadata.core.name != project_name:
        msg = f"Data in pyproject.toml is not related to {project_name}."
        raise GetVersionError(msg)
    return version


try:
    __version__ = _get_version_from_vcs(_PROJECT_NAME)
except (ImportError, LookupError, GetVersionError):
    import importlib.metadata

    __version__ = importlib.metadata.version(_PROJECT_NAME)
