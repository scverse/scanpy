"""Extension to skip inherited methods and properties in autosummary."""

from __future__ import annotations

from traceback import walk_stack
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def skip_inherited(  # noqa: PLR0917
    app: Sphinx,
    what: Literal["module", "class", "exception", "function", "method", "attribute"],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options | dict[str, object],
) -> bool | None:
    """Skip inherited members."""
    # Skip `getdoc` property
    if what == "method" and name == "getdoc":
        return True

    # find parent class
    for frame, _ in walk_stack(None):
        if frame.f_code.co_name == "_get_members" and frame.f_code.co_filename.endswith(
            "/generate.py"
        ):
            parent = frame.f_locals["obj"]
            if not isinstance(parent, type):
                return None
            break
    else:
        return None

    # skip if not a direct member of parent class
    if name not in parent.__dict__:
        return True
    return None


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.connect("autodoc-skip-member", skip_inherited)
