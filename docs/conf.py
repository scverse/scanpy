import sys
from pathlib import Path
from datetime import datetime

import matplotlib  # noqa
# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import scanpy  # noqa


# -- General configuration ------------------------------------------------


needs_sphinx = '1.7'  # autosummary bugfix

# General information
project = 'Scanpy'
author = scanpy.__author__
copyright = f'{datetime.now():%Y}, {author}.'
version = scanpy.__version__.replace('.dirty', '')
release = version

# default settings
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
default_role = 'literal'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    # 'plot_generator',
    # 'plot_directive',
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    # 'ipython_directive',
    # 'ipython_console_highlighting',
    'scanpydoc',
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False

intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/latest/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/', None),
    networkx=('https://networkx.github.io/documentation/networkx-1.10/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
)


# -- Options for HTML output ----------------------------------------------


html_theme = 'sphinx_rtd_theme'
html_theme_options = dict(
    navigation_depth=4,
    logo_only=True,           # Only show the logo
)
html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='theislab',   # Username
    github_repo='scanpy',     # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_logo = '_static/img/Scanpy_Logo_RGB.png'


def setup(app):
    app.add_stylesheet('css/custom.css')
    app.connect('autodoc-process-docstring', insert_function_images)


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [
    (master_doc, f'{project}.tex', doc_title, author, 'manual'),
]
man_pages = [
    (master_doc, project, doc_title, [author], 1)
]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, 'One line description of project.', 'Miscellaneous'),
]


# -- Images for plot functions -------------------------------------------------


def insert_function_images(app, what, name, obj, options, lines):
    path = Path(__file__).parent / 'api' / f'{name}.png'
    if what != 'function' or not path.is_file(): return
    lines[0:0] = [f'.. image:: {path.name}', '   :width: 200', '   :align: right', '']


# -- Test for new scanpydoc functionality --------------------------------------


import re
from sphinx.ext.napoleon import NumpyDocstring


def process_return(lines):
    for line in lines:
        m = re.fullmatch(r'(?P<param>\w+)\s+:\s+(?P<type>[\w.]+)', line)
        if m:
            # Once this is in scanpydoc, we can use the fancy hover stuff
            yield f'**{m["param"]}** : :class:`~{m["type"]}`'
        else:
            yield line


def scanpy_parse_returns_section(self, section):
    lines_raw = list(process_return(self._dedent(self._consume_to_next_section())))
    lines = self._format_block(':returns: ', lines_raw)
    if lines and lines[-1]:
        lines.append('')
    return lines


NumpyDocstring._parse_returns_section = scanpy_parse_returns_section


# -- Debug code ----------------------------------------------------------------


# Just do the following to see the rst of a function:
# rm -f _build/doctrees/api/scanpy.<what_you_want>.doctree; DEBUG=1 make html
import os
if os.environ.get('DEBUG') is not None:
    import sphinx.ext.napoleon
    pd = sphinx.ext.napoleon._process_docstring
    def pd_new(app, what, name, obj, options, lines):
        pd(app, what, name, obj, options, lines)
        print(*lines, sep='\n')
    sphinx.ext.napoleon._process_docstring = pd_new
