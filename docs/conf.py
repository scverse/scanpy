import ast
import sys
import inspect
import logging
from pathlib import Path
from datetime import datetime
from typing import List

from sphinx.application import Sphinx
from sphinx.ext import autosummary
from sphinx.ext.autosummary import limited_join

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

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosummary',
    # 'plot_generator',
    # 'plot_directive',
    'numpydoc',
    # 'ipython_directive',
    # 'ipython_console_highlighting',
]

# Generate the API documentation when building
autosummary_generate = True
# both of the following two lines don't work
# see falexwolf's issue for numpydoc
# autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False


templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'Scanpy'
author = 'Alex Wolf, Philipp Angerer, Tobias Callies, Davide Cittaro, Gokcen Eraslan'
copyright = f'{datetime.now():%Y}, {author}'

version = scanpy.__version__.replace('.dirty', '')
release = version
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------


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
        suffix=suffix, base_path=app.srcdir, imported_members=True,
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


def modurl(qualname):
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


# html_context doesn’t apply to autosummary templates ☹
# and there’s no way to insert filters into those templates
# so we have to modify the default filters
from jinja2.defaults import DEFAULT_FILTERS

DEFAULT_FILTERS['modurl'] = modurl


# -- Prettier Autodoc -----------------------------------------------------
# TODO: replace with sphinx.ext.napoleon


def unparse(ast_node: ast.expr, plain: bool = False) -> str:
    if isinstance(ast_node, ast.Attribute):
        if plain:
            return ast_node.attr
        else:
            v = unparse(ast_node.value, plain)
            return f'{v}.{ast_node.attr}'
    elif isinstance(ast_node, ast.Index):
        return unparse(ast_node.value)
    elif isinstance(ast_node, ast.Name):
        return ast_node.id
    elif isinstance(ast_node, ast.Subscript):
        v = unparse(ast_node.value, plain)
        s = unparse(ast_node.slice, plain)
        return f'{v}[{s}]'
    elif isinstance(ast_node, ast.Tuple):
        return ', '.join(unparse(e) for e in ast_node.elts)
    else:
        t = type(ast_node)
        raise NotImplementedError(f'can’t unparse {t}')


def mangle_signature(sig: str, max_chars: int = 30) -> str:
    fn = ast.parse(f'def f{sig}: pass').body[0]

    args_all = [a.arg for a in fn.args.args]
    n_a = len(args_all) - len(fn.args.defaults)
    args: List[str] = args_all[:n_a]
    opts: List[str] = args_all[n_a:]

    # Produce a more compact signature
    s = limited_join(', ', args, max_chars=max_chars - 2)
    if opts:
        if not s:
            opts_str = limited_join(', ', opts, max_chars=max_chars - 4)
            s = f'[{opts_str}]'
        elif len(s) < max_chars - 4 - 2 - 3:
            opts_str = limited_join(', ', opts, max_chars=max_chars - len(sig) - 4 - 2)
            s += f'[, {opts_str}]'

    if False:  # if fn.returns:  # do not show return type in docs
        ret = unparse(fn.returns, plain=True)
        return f'({s}) -> {ret}'
    return f'({s})'


autosummary.mangle_signature = mangle_signature


if __name__ == '__main__':
    print(mangle_signature('(filename: typing.Union[str, pathlib.Path], delim: int=0) -> anndata.base.AnnData'))
    print(mangle_signature('(a, *, b=1) -> int'))
    print(mangle_signature('(a, b=1, *c) -> Union[str, pathlib.Path]'))
    print(mangle_signature('(a, b=1, *c, d=1)'))
