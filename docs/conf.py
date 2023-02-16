import os
import sys
from pathlib import Path
from datetime import datetime

import matplotlib  # noqa
from packaging.version import parse as parse_version

# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / 'extensions')]
import scanpy  # noqa

on_rtd = os.environ.get('READTHEDOCS') == 'True'

# -- General configuration ------------------------------------------------


nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = '2.0'  # Nicer param docs
suppress_warnings = [
    'ref.citation',
    'myst.header',  # https://github.com/executablebooks/MyST-Parser/issues/262
]

# General information
project = 'Scanpy'
author = 'Scanpy development team'
copyright = f'{datetime.now():%Y}, the Scanpy development team.'
version = scanpy.__version__.replace('.dirty', '')

# Bumping the version updates all docs, so don't do that
if parse_version(version).is_devrelease:
    parsed = parse_version(version)
    version = f"{parsed.major}.{parsed.minor}.{parsed.micro}.dev"

release = version

# default settings
templates_path = ['_templates']
master_doc = 'index'
default_role = 'literal'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    # 'plot_generator',
    'matplotlib.sphinxext.plot_directive',
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    # 'ipython_directive',
    # 'ipython_console_highlighting',
    # 'scanpydoc',  # scanpydoc.elegant_typehints causes full doc rebuilds
    'scanpydoc.rtd_github_links',
    'scanpydoc.theme',
    'scanpydoc.definition_list_typed_field',
    'scanpydoc.autosummary_generate_imported',
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

typehints_defaults = 'braces'

scanpy_tutorials_url = 'https://scanpy-tutorials.readthedocs.io/en/latest/'

intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    cycler=('https://matplotlib.org/cycler/', None),
    h5py=('https://docs.h5py.org/en/stable/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/stable/', None),
    networkx=('https://networkx.org/documentation/stable/', None),
    numpy=('https://numpy.org/doc/stable/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    pytest=('https://docs.pytest.org/en/latest/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/', None),
    seaborn=('https://seaborn.pydata.org/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
    scanpy_tutorials=(scanpy_tutorials_url, None),
)


# -- Options for HTML output ----------------------------------------------


html_theme = 'scanpydoc'
html_theme_options = dict(
    navigation_depth=4,
    logo_only=True,
    docsearch_index='scanpy',
    docsearch_key='fa4304eb95d2134997e3729553a674b2',
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
    app.warningiserror = (
        False  # Temporarily disable warnings as errors to get 1.9.2 out
    )


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
    "sklearn.neighbors._dist_metrics.DistanceMetric": "sklearn.neighbors.DistanceMetric",
    # If the docs are built with an old version of numpy, this will make it work:
    "numpy.random.RandomState": "numpy.random.mtrand.RandomState",
    "scanpy.plotting._matrixplot.MatrixPlot": "scanpy.pl.MatrixPlot",
    "scanpy.plotting._dotplot.DotPlot": "scanpy.pl.DotPlot",
    "scanpy.plotting._stacked_violin.StackedViolin": "scanpy.pl.StackedViolin",
    "pandas.core.series.Series": "pandas.Series",
}

nitpick_ignore = [
    # Will probably be documented
    ('py:class', 'scanpy._settings.Verbosity'),
    # Currently undocumented: https://github.com/mwaskom/seaborn/issues/1810
    ('py:class', 'seaborn.ClusterGrid'),
    # Won’t be documented
    ('py:class', 'scanpy.plotting._utils._AxesSubplot'),
    ('py:class', 'scanpy._utils.Empty'),
    ('py:class', 'numpy.random.mtrand.RandomState'),
    # Will work once scipy 1.8 is released
    ('py:class', 'scipy.sparse.base.spmatrix'),
    ('py:class', 'scipy.sparse.csr.csr_matrix'),
]

# Options for plot examples

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root
