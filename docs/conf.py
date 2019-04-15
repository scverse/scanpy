import sys
import logging
from pathlib import Path
from datetime import datetime

from jinja2.defaults import DEFAULT_FILTERS

import matplotlib  # noqa
# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import scanpy  # noqa


logger = logging.getLogger(__name__)


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


def api_image(qualname: str) -> str:
    # I’d like to make this a contextfilter, but the jinja context doesn’t contain the path,
    # so no chance to not hardcode “api/” here.
    path = Path(__file__).parent / 'api' / f'{qualname}.png'
    print(path, path.is_file())
    return f'.. image:: {path.name}\n   :width: 200\n   :align: right' if path.is_file() else ''


# html_context doesn’t apply to autosummary templates ☹
# and there’s no way to insert filters into those templates
# so we have to modify the default filters
DEFAULT_FILTERS['api_image'] = api_image


# -- Test for new scanpydoc functionality --------------------------------------

from itertools import repeat, chain
from sphinx.ext.napoleon import NumpyDocstring

# # allow "prose" sections...
# def _consume_returns_section(self):
#     # type: () -> List[Tuple[unicode, unicode, List[unicode]]]
#     # this is the full original function
#     fields = self._consume_fields(prefer_type=True)
#     # Let us postprocess the output of this function.
#     # 
#     # fields with empty descriptions are prose fields, the "actual descriptions"
#     # are stored in the types, hence:
#     #
#     # concat fields with empty descriptions
#     #
#     new_fields = []
#     concat_with_old = False
#     for field in fields:
#         name, type, descr = field
#         if (not descr  # empty description (empty list)
#             or len(descr) == 1 and descr[0] == ''):  # empty description (empty string)
#             new_descr = ''
#             if name != '':
#                 new_descr = name + ': '
#             new_descr += type + '\n'
#             # deal with escaped *
#             new_descr = new_descr.replace('\* ', '* ')
#             if concat_with_old:
#                 # concat to the description section
#                 new_fields[-1][2].append(new_descr)
#             else:
#                 new_fields.append(('', '', [new_descr]))
#                 concat_with_old = True
#         else:
#             new_fields.append(field)
#             concat_with_old = False
#     return new_fields


# This is essentially entirely copied, the only change here is removing the bullets.
def _parse_returns_section(self, section):
    # type: (unicode) -> List[unicode]
    fields = self._consume_returns_section()
    multi = len(fields) > 1
    if multi:
        use_rtype = False
    else:
        use_rtype = self._config.napoleon_use_rtype

    lines = []  # type: List[unicode]
    for _name, _type, _desc in fields:
        if use_rtype:
            field = self._format_field(_name, '', _desc)
        else:
            field = self._format_field(_name, _type, _desc)

        if multi:
            if lines:
                lines.extend(self._format_block('          ', field))
            else:
                lines.extend(self._format_block(':returns: ', field))
        else:
            lines.extend(self._format_block(':returns: ', field))
            if _type and use_rtype:
                lines.extend([':rtype: %s' % _type, ''])
    if lines and lines[-1]:
        lines.append('')
    return lines


NumpyDocstring._consume_returns_section = _consume_returns_section
NumpyDocstring._parse_returns_section = _parse_returns_section
