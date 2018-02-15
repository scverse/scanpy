#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Scanpy documentation build configuration file, created by
# sphinx-quickstart on Sun Aug 20 00:29:31 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import ast
import os
import sys
import time
import inspect
from pathlib import Path
sys.path.insert(0, os.path.abspath(os.path.pardir))
import logging

logger = logging.getLogger(__name__)
from sphinx.ext import autosummary, autodoc
from sphinx.ext.autosummary import limited_join

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
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
autodoc_mock_imports = ['_tkinter']
# both of the following two lines don't work
# see falexwolf's issue for numpydoc
# autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False

def process_generate_options(app):
    # type: (Sphinx) -> None
    genfiles = app.config.autosummary_generate
    if genfiles and not hasattr(genfiles, '__len__'):
        env = app.builder.env
        genfiles = [env.doc2path(x, base=None) for x in env.found_docs
                    if os.path.isfile(env.doc2path(x))]
    if not genfiles:
        return
    from sphinx.ext.autosummary.generate import generate_autosummary_docs
    ext = app.config.source_suffix
    genfiles = [genfile + (not genfile.endswith(tuple(ext)) and ext[0] or '')
                for genfile in genfiles]
    suffix = autosummary.get_rst_suffix(app)
    if suffix is None:
        return
    generate_autosummary_docs(genfiles, builder=app.builder,
                              warn=logger.warning, info=logger.info,
                              suffix=suffix, base_path=app.srcdir, imported_members=True)

autosummary.process_generate_options = process_generate_options

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'Scanpy'
author = 'Alex Wolf, Philipp Angerer, Tobias Callies, Davide Cittaro'
copyright = '{}, {}'.format(time.strftime("%Y"), author)

import scanpy
version = scanpy.__version__.replace('.dirty', '')
release = version
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
if html_theme == 'sphinx_rtd_theme':
    html_theme_options = {
        'navigation_depth': 2,
    }
    html_context = {
        "display_github": True,  # Integrate GitHub
        "github_user": "theislab",  # Username
        "github_repo": "scanpy",  # Repo name
        "github_version": "master",  # Version
        "conf_py_path": "/docs/",  # Path in the checkout to the docs root
    }
elif html_theme == 'bootstrap':
    import sphinx_bootstrap_theme
    html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
    html_theme_options = {
        'navbar_site_name': "Site",
        'navbar_pagenav_name': "Page",
        'source_link_position': "footer",
        'bootswatch_theme': "paper",
        'navbar_pagenav': False,
        'navbar_sidebarrel': False,
        'bootstrap_version': "3",
        'navbar_links': [
            ("API", "api"),
        ],
    }

html_static_path = ['_static']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Scanpydoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Scanpy.tex', 'Scanpy Documentation',
     'Alex Wolf, Philipp Angerer', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'scanpy', 'Scanpy Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Scanpy', 'Scanpy Documentation',
     author, 'Scanpy', 'One line description of project.',
     'Miscellaneous'),
]


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
    fragment = '#L{}-L{}'.format(start, end) if start and end else ''
    return '{}/{}{}'.format(github_url, path, fragment)


# html_context doesn’t apply to autosummary templates ☹
# and there’s no way to insert filters into those templates
# so we have to modify the default filters
from jinja2.defaults import DEFAULT_FILTERS

DEFAULT_FILTERS['modurl'] = modurl


# -- Prettier Autodoc -----------------------------------------------------


def f(string):
    frame = sys._getframe(1)
    return string.format_map(frame.f_locals)


def unparse(ast_node: ast.expr, plain: bool=False) -> str:
    if isinstance(ast_node, ast.Attribute):
        if plain:
            return ast_node.attr
        else:
            v = unparse(ast_node.value, plain)
            return f('{v}.{ast_node.attr}')
    elif isinstance(ast_node, ast.Index):
        return unparse(ast_node.value)
    elif isinstance(ast_node, ast.Name):
        return ast_node.id
    elif isinstance(ast_node, ast.Subscript):
        v = unparse(ast_node.value, plain)
        s = unparse(ast_node.slice, plain)
        return f('{v}[{s}]')
    elif isinstance(ast_node, ast.Tuple):
        return ', '.join(unparse(e) for e in ast_node.elts)
    else:
        t = type(ast_node)
        raise NotImplementedError(f('can’t unparse {t}'))


def mangle_signature(sig: str, max_chars: int=30) -> str:
    fn = ast.parse(f('def f{sig}: pass')).body[0]

    args_all = [a.arg for a in fn.args.args]
    n_a = len(args_all) - len(fn.args.defaults)
    args = args_all[:n_a]  # type: List[str]
    opts = args_all[n_a:]  # type: List[str]

    # Produce a more compact signature
    s = limited_join(', ', args, max_chars=max_chars - 2)
    if opts:
        if not s:
            opts_str = limited_join(', ', opts, max_chars=max_chars - 4)
            s = f('[{opts_str}]')
        elif len(s) < max_chars - 4 - 2 - 3:
            opts_str = limited_join(', ', opts, max_chars=max_chars - len(sig) - 4 - 2)
            s += f('[, {opts_str}]')

    if False:  # fn.returns:  # do not show return type in docs
        ret = unparse(fn.returns, plain=True)
        return f('({s}) -> {ret}')
    return f('({s})')


autosummary.mangle_signature = mangle_signature

# TODO: also replace those for individual function pages:
# autodoc.formatargspec
# autodoc.format_annotation


if __name__ == '__main__':
    print(mangle_signature('(filename: typing.Union[str, pathlib.Path], delim: int=0) -> anndata.base.AnnData'))
    print(mangle_signature('(a, *, b=1) -> int'))
    print(mangle_signature('(a, b=1, *c) -> Union[str, pathlib.Path]'))
    print(mangle_signature('(a, b=1, *c, d=1)'))
