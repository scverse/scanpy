"""Extension to skip deprecated methods and properties in autosummary."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def skip_deprecated(  # noqa: PLR0917
    app: Sphinx,
    what: Literal[
        "module", "class", "exception", "function", "method", "attribute", "property"
    ],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options | dict[str, object],
) -> bool | None:
    """Skip deprecated members."""
    if hasattr(obj, "__deprecated__"):
        return True
    return None


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.connect("autodoc-skip-member", skip_deprecated)
