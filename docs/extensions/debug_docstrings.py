# Just do the following to see the rst of a function:
# rm -f _build/doctrees/api/scanpy.<what_you_want>.doctree; DEBUG=1 make html
import os

from sphinx.application import Sphinx
import sphinx.ext.napoleon


_pd_orig = sphinx.ext.napoleon._process_docstring


def pd_new(app, what, name, obj, options, lines):
    _pd_orig(app, what, name, obj, options, lines)
    print(*lines, sep='\n')


def setup(app: Sphinx):
    if os.environ.get('DEBUG') is not None:
        sphinx.ext.napoleon._process_docstring = pd_new
