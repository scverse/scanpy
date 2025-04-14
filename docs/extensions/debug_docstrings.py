"""Extension for debugging docstrings."""

# Just do the following to see the rst of a function:
# rm ./_build/doctrees/api/generated/scanpy.<what you want>.doctree; DEBUG=1 make html
from __future__ import annotations

import os
from typing import TYPE_CHECKING

import sphinx.ext.napoleon

if TYPE_CHECKING:
    from sphinx.application import Sphinx

_pd_orig = sphinx.ext.napoleon._process_docstring


def pd_new(app, what, name, obj, options, lines) -> None:  # noqa: PLR0917
    """Wrap ``sphinx.ext.napoleon._process_docstring``."""
    _pd_orig(app, what, name, obj, options, lines)
    print(*lines, sep="\n")


def setup(app: Sphinx) -> None:
    """App setup hook."""
    if os.environ.get("DEBUG") is not None:
        sphinx.ext.napoleon._process_docstring = pd_new
