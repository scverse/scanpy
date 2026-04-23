"""Extension to patch ignored mime types."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from importlib.abc import MetaPathFinder
from importlib.metadata import Distribution, EntryPoint, EntryPoints
from types import MappingProxyType
from typing import TYPE_CHECKING, override

from myst_nb.core.render import MimeRenderPlugin
from sphinx.util.typing import ExtensionMetadata

if TYPE_CHECKING:
    import os
    from collections.abc import Iterable, Sequence
    from importlib.machinery import ModuleSpec
    from importlib.metadata import DistributionFinder, SimplePath
    from types import ModuleType

    from docutils import nodes
    from myst_nb.core.render import MimeData, NbElementRenderer
    from sphinx.application import Sphinx


ignore: set[str] = set()


class _Ignore(MimeRenderPlugin):
    @override
    @staticmethod
    def handle_mime(
        renderer: NbElementRenderer, data: MimeData, inline: bool
    ) -> None | list[nodes.Element]:
        if data.mime_type in ignore:
            return []  # returning a list instead of `None` means “we handled it”
        return None


class _IgnoreMimeDist(Distribution):
    metadata = MappingProxyType(dict(Name=__name__, Version="0.0.0"))

    @override
    def read_text(self, filename: str) -> str | None:
        return None

    @override
    def locate_file(self, path: str | os.PathLike[str]) -> SimplePath:
        raise RuntimeError

    @property
    @override
    def entry_points(self) -> EntryPoints:
        ep = EntryPoint("ignore", f"{__name__}:_Ignore", "myst_nb.mime_renderers")
        return EntryPoints([ep])


@dataclass
class _IgnoreMimeFinder(MetaPathFinder):
    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None:
        return None

    def find_distributions(
        self, context: DistributionFinder.Context | None
    ) -> Iterable[Distribution]:
        """Find fake distribution."""
        yield _IgnoreMimeDist()


def setup(app: Sphinx) -> ExtensionMetadata:
    """App setup hook."""
    global ignore  # noqa: PLW0603

    app.add_config_value("myst_ignore_mime_types", [], "env")
    ignore |= set(app.config.myst_ignore_mime_types)

    sys.meta_path.append(_IgnoreMimeFinder())

    return ExtensionMetadata(parallel_read_safe=True)
