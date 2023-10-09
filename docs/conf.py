import os
import sys
from pathlib import Path
from datetime import datetime
from typing import Any

import matplotlib  # noqa
from sphinx.application import Sphinx
from packaging.version import parse as parse_version

# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / 'extensions')]
import scanpy  # noqa


# -- General configuration ------------------------------------------------


nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = '4.0'  # Nicer param docs
suppress_warnings = [
    'ref.citation',
    'myst.header',  # https://github.com/executablebooks/MyST-Parser/issues/262
]

# General information
project = 'Scanpy'
author = 'Scanpy development team'
repository_url = "https://github.com/scverse/scanpy"
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
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

extensions = [
    'myst_nb',
    'sphinx_copybutton',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    'matplotlib.sphinxext.plot_directive',
    'sphinx_autodoc_typehints',  # needs to be after napoleon
    'git_ref',  # needs to be before scanpydoc.rtd_github_links
    'scanpydoc',  # needs to be before sphinx.ext.linkcode
    'sphinx.ext.linkcode',
    'sphinx_design',
    'sphinxext.opengraph',
    *[p.stem for p in (HERE / 'extensions').glob('*.py') if p.stem not in {'git_ref'}],
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
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True


ogp_site_url = "https://scanpy.readthedocs.io/en/stable/"
ogp_image = "https://scanpy.readthedocs.io/en/stable/_static/Scanpy_Logo_BrightFG.svg"

typehints_defaults = 'braces'

pygments_style = "default"
pygments_dark_style = "native"

intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    cuml=('https://docs.rapids.ai/api/cuml/stable/', None),
    cycler=('https://matplotlib.org/cycler/', None),
    dask=('https://docs.dask.org/en/stable/', None),
    dask_ml=('https://ml.dask.org/', None),
    h5py=('https://docs.h5py.org/en/stable/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    igraph=('https://python.igraph.org/en/stable/api/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/stable/', None),
    networkx=('https://networkx.org/documentation/stable/', None),
    numpy=('https://numpy.org/doc/stable/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    pynndescent=('https://pynndescent.readthedocs.io/en/latest/', None),
    pytest=('https://docs.pytest.org/en/latest/', None),
    python=('https://docs.python.org/3', None),
    rapids_singlecell=('https://rapids-singlecell.readthedocs.io/en/latest/', None),
    scipy=('https://docs.scipy.org/doc/scipy/', None),
    seaborn=('https://seaborn.pydata.org/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
    tutorials=('https://scanpy-tutorials.readthedocs.io/en/latest/', None),
)


# -- Options for HTML output ----------------------------------------------

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
}
html_static_path = ['_static']
html_css_files = ["css/override.css"]
html_show_sphinx = False
html_logo = '_static/img/Scanpy_Logo_BrightFG.svg'
html_title = "scanpy"


def setup(app: Sphinx):
    """App setup hook."""
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
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
    ('py:class', 'scanpy.neighbors.OnFlySymMatrix'),
    # Currently undocumented: https://github.com/mwaskom/seaborn/issues/1810
    ('py:class', 'seaborn.ClusterGrid'),
    # Won’t be documented
    ('py:class', 'scanpy.plotting._utils._AxesSubplot'),
    ('py:class', 'scanpy._utils.Empty'),
    ('py:class', 'numpy.random.mtrand.RandomState'),
    ('py:class', 'scanpy.neighbors._types.KnnTransformerLike'),
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

# extlinks config
extlinks = {
    "issue": ("https://github.com/scverse/scanpy/issues/%s", "issue%s"),
    "pr": ("https://github.com/scverse/scanpy/pull/%s", "pr%s"),
}
