<<<<<<< HEAD
import sys
import logging
from pathlib import Path
from datetime import datetime

from jinja2.defaults import DEFAULT_FILTERS

import matplotlib  # noqa
=======
import os
import sys
import warnings
from pathlib import Path
from datetime import datetime

import matplotlib  # noqa

>>>>>>> upstream/master
# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
<<<<<<< HEAD
sys.path.insert(0, str(HERE.parent))
import scanpy  # noqa


logger = logging.getLogger(__name__)

=======
sys.path[:0] = [str(HERE.parent), str(HERE / 'extensions')]
import scanpy  # noqa

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=FutureWarning)
    import scanpy.api

on_rtd = os.environ.get('READTHEDOCS') == 'True'
>>>>>>> upstream/master

# -- General configuration ------------------------------------------------


<<<<<<< HEAD
needs_sphinx = '1.7'  # autosummary bugfix
=======
nitpicky = True  # Warn about broken links
needs_sphinx = '2.0'  # Nicer param docs
suppress_warnings = ['ref.citation']
>>>>>>> upstream/master

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
<<<<<<< HEAD

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
)
html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='theislab',   # Username
    github_repo='scanpy',     # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
)
html_static_path = ['_static']


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
=======

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
    *[p.stem for p in (HERE / 'extensions').glob('*.py')],
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
api_dir = HERE / 'api'  # function_images

intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    cycler=('https://matplotlib.org/cycler/', None),
    h5py=('http://docs.h5py.org/en/stable/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/', None),
    networkx=('https://networkx.github.io/documentation/networkx-1.10/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    scvelo=('https://scvelo.readthedocs.io/', None),
    seaborn=('https://seaborn.pydata.org/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
    scanpy_tutorials=(
        'https://scanpy-tutorials.readthedocs.io/en/latest',
        None,
    ),
)


# -- Options for HTML output ----------------------------------------------


html_theme = 'sphinx_rtd_theme'
html_theme_options = dict(
    navigation_depth=4, logo_only=True  # Only show the logo
)
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user='theislab',  # Username
    github_repo='scanpy',  # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',  # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_show_sphinx = False
html_logo = '_static/img/Scanpy_Logo_BrightFG.svg'


def setup(app):
    app.warningiserror = on_rtd
    app.add_stylesheet('css/custom.css')


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [(master_doc, f'{project}.tex', doc_title, author, 'manual')]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        'One line description of project.',
        'Miscellaneous',
    )
]


# -- Suppress link warnings ----------------------------------------------------

qualname_overrides = {
    "sklearn.neighbors.dist_metrics.DistanceMetric": "sklearn.neighbors.DistanceMetric"
}

nitpick_ignore = [
    # Will probably be documented
    ('py:class', 'scanpy._settings.Verbosity'),
    # Currently undocumented: https://github.com/mwaskom/seaborn/issues/1810
    ('py:class', 'seaborn.ClusterGrid'),
    # Won’t be documented
    ('py:class', 'scanpy.readwrite.Empty'),
]

for mod_name in [
    'pp',
    'tl',
    'pl',
    'queries',
    'logging',
    'datasets',
    'export_to',
    None,
]:
    if mod_name is None:
        mod = scanpy.api
        mod_name = 'scanpy.api'
    else:
        mod = getattr(scanpy.api, mod_name)
        mod_name = f'scanpy.api.{mod_name}'
    for name, item in vars(mod).items():
        if not callable(item):
            continue
        for kind in ['func', 'obj']:
            nitpick_ignore.append((f'py:{kind}', f'{mod_name}.{name}'))
>>>>>>> upstream/master
