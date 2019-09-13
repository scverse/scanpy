import os
import sys
import warnings
from pathlib import Path
from datetime import datetime

import matplotlib  # noqa

# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import scanpy  # noqa

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=FutureWarning)
    import scanpy.api

on_rtd = os.environ.get('READTHEDOCS') == 'True'

# -- General configuration ------------------------------------------------


nitpicky = True  # Warn about broken links
needs_sphinx = '2.0'  # Nicer param docs
suppress_warnings = ['ref.citation']

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
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    cycler=('https://matplotlib.org/cycler/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/', None),
    networkx=('https://networkx.github.io/documentation/networkx-1.10/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    scvelo=('https://scvelo.readthedocs.io/en/stable/', None),
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
html_logo = '_static/img/Scanpy_Logo_RGB.png'
gh_url = 'https://github.com/{github_user}/{github_repo}'.format_map(
    html_context
)


def setup(app):
    app.warningiserror = on_rtd
    app.add_stylesheet('css/custom.css')
    app.connect('autodoc-process-docstring', insert_function_images)
    app.connect('build-finished', show_param_warnings)
    app.add_role('pr', autolink(f'{gh_url}/pull/{{}}', 'PR {}'))


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


# -- Images for plot functions -------------------------------------------------


def insert_function_images(app, what, name, obj, options, lines):
    path = Path(__file__).parent / 'api' / f'{name}.png'
    if what != 'function' or not path.is_file():
        return
    lines[0:0] = [
        f'.. image:: {path.name}',
        '   :width: 200',
        '   :align: right',
        '',
    ]


# -- GitHub links --------------------------------------------------------------


def autolink(url_template, title_template='{}'):
    from docutils import nodes

    def role(name, rawtext, text, lineno, inliner, options={}, content=[]):
        url = url_template.format(text)
        title = title_template.format(text)
        node = nodes.reference(rawtext, title, refuri=url, **options)
        return [node], []

    return role


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
    lines_raw = list(
        process_return(self._dedent(self._consume_to_next_section()))
    )
    lines = self._format_block(':returns: ', lines_raw)
    if lines and lines[-1]:
        lines.append('')
    return lines


NumpyDocstring._parse_returns_section = scanpy_parse_returns_section


# -- Warn for non-annotated params ---------------------------------------------


_format_docutils_params_orig = NumpyDocstring._format_docutils_params
param_warnings = {}


def scanpy_log_param_types(self, fields, field_role='param', type_role='type'):
    for _name, _type, _desc in fields:
        if not _type:
            continue
        set_item = r"`'[a-z0-9_.-]+'`"
        if re.fullmatch(rf"{{{set_item}(, {set_item})*}}", _type):
            continue
        param_warnings.setdefault((self._name, self._obj), []).append(
            (_name, _type)
        )
    return _format_docutils_params_orig(self, fields, field_role, type_role)


def show_param_warnings(app, exception):
    import inspect

    for (fname, fun), params in param_warnings.items():
        _, line = inspect.getsourcelines(fun)
        file_name = inspect.getsourcefile(fun)
        params_str = '\n'.join(f'\t{n}: {t}' for n, t in params)
        warnings.warn_explicit(
            f'\nParameters in `{fname}` not set-like: {{`elm-1`, `s_el.2`}}.\n'
            'Convert to this format or replace with type annotations:\n'
            + params_str,
            UserWarning,
            file_name,
            line,
        )


NumpyDocstring._format_docutils_params = scanpy_log_param_types


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


# -- Suppress link warnings ----------------------------------------------------

qualname_overrides = {
    "sklearn.neighbors.dist_metrics.DistanceMetric": "sklearn.neighbors.DistanceMetric"
}

nitpick_ignore = [
    # Will probably be documented
    ('py:class', 'scanpy._settings.Verbosity'),
    # Currently undocumented: https://github.com/mwaskom/seaborn/issues/1810
    ('py:class', 'seaborn.ClusterGrid'),
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
