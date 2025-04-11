"""Images for plot functions."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def insert_function_images(  # noqa: PLR0917
    app: Sphinx, what: str, name: str, obj: Any, options: Options, lines: list[str]
) -> None:
    """Insert images for plot functions."""
    path = app.config.api_dir / f"{name}.png"
    if what != "function" or not path.is_file():
        return
    lines[0:0] = [
        f".. image:: {path.name}",
        "   :width: 200",
        "   :align: right",
        "",
    ]


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_config_value("api_dir", Path(), "env")
    app.connect("autodoc-process-docstring", insert_function_images)
