# Just do the following to see the rst of a function:
# rm -f _build/doctrees/api/scanpy.<what_you_want>.doctree; DEBUG=1 make html
from __future__ import annotations

import os
from typing import TYPE_CHECKING

import sphinx.ext.napoleon

if TYPE_CHECKING:
    from sphinx.application import Sphinx

_pd_orig = sphinx.ext.napoleon._process_docstring


def pd_new(app, what, name, obj, options, lines):
    _pd_orig(app, what, name, obj, options, lines)
    print(*lines, sep="\n")


def setup(app: Sphinx):
    if os.environ.get("DEBUG") is not None:
        sphinx.ext.napoleon._process_docstring = pd_new
