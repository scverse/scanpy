import sys
import inspect
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional

from sphinx.application import Sphinx
from sphinx.ext import autosummary

# remove PyCharm’s old six module
if 'six' in sys.modules:
    print(*sys.path, sep='\n')
    for pypath in list(sys.path):
        if any(p in pypath for p in ['PyCharm', 'pycharm']) and 'helpers' in pypath:
            sys.path.remove(pypath)
    del sys.modules['six']

import matplotlib  # noqa
# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use('agg')

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import scanpy.api  # noqa


logger = logging.getLogger(__name__)

# -- General configuration ------------------------------------------------

needs_sphinx = '1.7'  # autosummary bugfix
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosummary',
    # 'plot_generator',
    # 'plot_directive',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.intersphinx',
    # 'ipython_directive',
    # 'ipython_console_highlighting',
]

# Generate the API documentation when building
autosummary_generate = True
# both of the following two lines don't work
# see falexwolf's issue for numpydoc
# autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = False
napoleon_custom_sections = [('Params', 'Parameters')]

intersphinx_mapping = dict(
    python=('https://docs.python.org/3', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    matplotlib=('https://matplotlib.org/', None),
    anndata=('https://anndata.readthedocs.io/en/latest/', None),
)

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'Scanpy'
author = 'Alex Wolf, Philipp Angerer, Davide Cittaro, Gokcen Eraslan, Fidel Ramirez, Tobias Callies'
copyright = f'{datetime.now():%Y}, {author}'

version = scanpy.__version__.replace('.dirty', '')
release = version
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = dict(
    navigation_depth=2,
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


htmlhelp_basename = 'Scanpydoc'
latex_documents = [
    (master_doc, 'Scanpy.tex', 'Scanpy Documentation',
     'Alex Wolf, Philipp Angerer', 'manual'),
]
man_pages = [
    (master_doc, 'scanpy', 'Scanpy Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'Scanpy', 'Scanpy Documentation',
     author, 'Scanpy', 'One line description of project.',
     'Miscellaneous'),
]


# -- generate_options override ------------------------------------------
# TODO: why?


def process_generate_options(app: Sphinx):
    genfiles = app.config.autosummary_generate

    if genfiles and not hasattr(genfiles, '__len__'):
        env = app.builder.env
        genfiles = [
            env.doc2path(x, base=None)
            for x in env.found_docs
            if Path(env.doc2path(x)).is_file()
        ]
    if not genfiles:
        return

    from sphinx.ext.autosummary.generate import generate_autosummary_docs

    ext = app.config.source_suffix
    genfiles = [
        genfile + (not genfile.endswith(tuple(ext)) and ext[0] or '')
        for genfile in genfiles
    ]

    suffix = autosummary.get_rst_suffix(app)
    if suffix is None:
        return

    generate_autosummary_docs(
        genfiles, builder=app.builder,
        warn=logger.warning, info=logger.info,
        suffix=suffix, base_path=app.srcdir,
        imported_members=True, app=app,
    )


autosummary.process_generate_options = process_generate_options


# -- GitHub URLs for class and method pages ------------------------------------------


def get_obj_module(qualname):
    """Get a module/class/attribute and its original module by qualname"""
    modname = qualname
    classname = None
    attrname = None
    while modname not in sys.modules:
        attrname = classname
        modname, classname = modname.rsplit('.', 1)

    # retrieve object and find original module name
    if classname:
        cls = getattr(sys.modules[modname], classname)
        modname = cls.__module__
        obj = getattr(cls, attrname) if attrname else cls
    else:
        obj = None

    return obj, sys.modules[modname]


def get_linenos(obj):
    """Get an object’s line numbers"""
    try:
        lines, start = inspect.getsourcelines(obj)
    except TypeError:
        return None, None
    else:
        return start, start + len(lines) - 1


project_dir = Path(__file__).parent.parent  # project/docs/conf.py/../.. → project/
github_url1 = 'https://github.com/{github_user}/{github_repo}/tree/{github_version}'.format_map(html_context)
github_url2 = 'https://github.com/theislab/anndata/tree/master'


def modurl(qualname: str) -> str:
    """Get the full GitHub URL for some object’s qualname."""
    obj, module = get_obj_module(qualname)
    github_url = github_url1
    try:
        path = Path(module.__file__).relative_to(project_dir)
    except ValueError:
        # trying to document something from another package
        github_url = github_url2
        path = '/'.join(module.__file__.split('/')[-2:])
    start, end = get_linenos(obj)
    fragment = f'#L{start}-L{end}' if start and end else ''
    return f'{github_url}/{path}{fragment}'


def api_image(qualname: str) -> Optional[str]:
    # I’d like to make this a contextfilter, but the jinja context doesn’t contain the path,
    # so no chance to not hardcode “api/” here.
    path = Path(__file__).parent / 'api' / f'{qualname}.png'
    print(path, path.is_file())
    return f'.. image:: {path.name}\n   :width: 200\n   :align: right' if path.is_file() else ''


# html_context doesn’t apply to autosummary templates ☹
# and there’s no way to insert filters into those templates
# so we have to modify the default filters
from jinja2.defaults import DEFAULT_FILTERS

DEFAULT_FILTERS.update(modurl=modurl, api_image=api_image)


# -- Prettier Param docs --------------------------------------------


from typing import Dict, List, Tuple

from docutils import nodes
from sphinx import addnodes
from sphinx.domains.python import PyTypedField, PyObject
from sphinx.environment import BuildEnvironment


class PrettyTypedField(PyTypedField):
    list_type = nodes.definition_list

    def make_field(
        self,
        types: Dict[str, List[nodes.Node]],
        domain: str,
        items: Tuple[str, List[nodes.inline]],
        env: BuildEnvironment = None
    ) -> nodes.field:
        def makerefs(rolename, name, node):
            return self.make_xrefs(rolename, domain, name, node, env=env)

        def handle_item(fieldarg: str, content: List[nodes.inline]) -> nodes.definition_list_item:
            head = nodes.term()
            head += makerefs(self.rolename, fieldarg, addnodes.literal_strong)
            fieldtype = types.pop(fieldarg, None)
            if fieldtype is not None:
                head += nodes.Text(' : ')
                if len(fieldtype) == 1 and isinstance(fieldtype[0], nodes.Text):
                    typename = ''.join(n.astext() for n in fieldtype)
                    head += makerefs(self.typerolename, typename, addnodes.literal_emphasis)
                else:
                    head += fieldtype

            body_content = nodes.paragraph('', '', *content)
            body = nodes.definition('', body_content)

            return nodes.definition_list_item('', head, body)

        fieldname = nodes.field_name('', self.label)
        if len(items) == 1 and self.can_collapse:
            fieldarg, content = items[0]
            bodynode = handle_item(fieldarg, content)
        else:
            bodynode = self.list_type()
            for fieldarg, content in items:
                bodynode += handle_item(fieldarg, content)
        fieldbody = nodes.field_body('', bodynode)
        return nodes.field('', fieldname, fieldbody)


# replace matching field types with ours
PyObject.doc_field_types = [
    PrettyTypedField(
        ft.name,
        names=ft.names,
        typenames=ft.typenames,
        label=ft.label,
        rolename=ft.rolename,
        typerolename=ft.typerolename,
        can_collapse=ft.can_collapse,
    ) if isinstance(ft, PyTypedField) else ft
    for ft in PyObject.doc_field_types
]
