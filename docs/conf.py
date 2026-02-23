"""Configuration for Scanpy’s Sphinx documentation."""

from __future__ import annotations

import os
import shutil
import sys
from datetime import datetime
from functools import partial
from importlib.metadata import version as get_version
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

import matplotlib  # noqa
from docutils import nodes
from packaging.version import Version
from sphinxcontrib.katex import NODEJS_BINARY

# Don’t use tkinter agg when importing scanpy → … → matplotlib
matplotlib.use("agg")

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]
os.environ["SPHINX_RUNNING"] = "1"  # for scanpy._singleton

if TYPE_CHECKING:
    from sphinx.application import Sphinx


# -- General configuration ------------------------------------------------

nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = "4.0"  # Nicer param docs

# General information
project = "Scanpy"
author = "Scanpy development team"
repository_url = "https://github.com/scverse/scanpy"
copyright = f"{datetime.now():%Y}, scverse"
version = get_version("scanpy").replace(".dirty", "")

# Bumping the version updates all docs, so don't do that
if Version(version).is_devrelease:
    parsed = Version(version)
    version = f"{parsed.major}.{parsed.minor}.{parsed.micro}.dev"

release = version

# Bibliography settings
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

# default settings
templates_path = ["_templates"]
master_doc = "index"
default_role = "literal"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    # exclude all 0.x.y.md files, but not index.md
    "release-notes/[!i]*.md",
]

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.katex",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "git_ref",  # needs to be before scanpydoc.rtd_github_links
    "scanpydoc",  # needs to be before sphinx.ext.linkcode
    "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinx_issues",
    "sphinxext.opengraph",
    *[p.stem for p in (HERE / "extensions").glob("*.py") if p.stem not in {"git_ref"}],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_typehints = "none"
autodoc_member_order = "bysource"
autodoc_default_options = {
    # Don’t show members in addition to the autosummary table added by `_templates/class.rst`
    "members": False,
    # show “Bases: SomeClass” at the top of class docs
    "show-inheritance": True,
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
api_dir = HERE / "api"  # function_images
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto", "ftp")
myst_heading_anchors = 3
myst_ignore_mime_types = [  # from custom extension patch_myst_nb
    "application/vnd.microsoft.datawrangler.viewer.v0+json",
]
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True


ogp_site_url = "https://scanpy.readthedocs.io/en/stable/"
ogp_image = "https://scanpy.readthedocs.io/en/stable/_static/Scanpy_Logo_BrightFG.svg"

typehints_defaults = "braces"
always_use_bars_union = True  # Don’t use `Union` even when building with Python ≤3.14

pygments_style = "default"
pygments_dark_style = "native"

katex_prerender = shutil.which(NODEJS_BINARY) is not None

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    bbknn=("https://bbknn.readthedocs.io/en/latest/", None),
    cuml=("https://docs.rapids.ai/api/cuml/stable/", None),
    cycler=("https://matplotlib.org/cycler/", None),
    dask=("https://docs.dask.org/en/stable/", None),
    dask_ml=("https://ml.dask.org/", None),
    decoupler=("https://decoupler.readthedocs.io/en/stable/", None),
    fast_array_utils=(
        "https://icb-fast-array-utils.readthedocs-hosted.com/en/stable/",
        None,
    ),
    h5py=("https://docs.h5py.org/en/stable/", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    ipython=("https://ipython.readthedocs.io/en/stable/", None),
    igraph=("https://python.igraph.org/en/stable/api/", None),
    leidenalg=("https://leidenalg.readthedocs.io/en/latest/", None),
    louvain=("https://louvain-igraph.readthedocs.io/en/latest/", None),
    matplotlib=("https://matplotlib.org/stable/", None),
    networkx=("https://networkx.org/documentation/stable/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    pydeseq2=("https://pydeseq2.readthedocs.io/en/stable/", None),
    pynndescent=("https://pynndescent.readthedocs.io/en/latest/", None),
    pytest=("https://docs.pytest.org/en/latest/", None),
    python=("https://docs.python.org/3", None),
    rapids_singlecell=("https://rapids-singlecell.readthedocs.io/en/latest/", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    session_info2=("https://session-info2.readthedocs.io/en/stable/", None),
    squidpy=("https://squidpy.readthedocs.io/en/stable/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
)


array_support: dict[str, tuple[list[str], list[str]]] = {
    "experimental.pp.highly_variable_genes": (["np", "sp"], []),
    "get.aggregate": (["np", "sp", "da"], []),
    "pp.calculate_qc_metrics": (["np", "sp", "da"], []),
    "pp.combat": (["np"], []),
    "pp.downsample_counts": (["np", "sp[csr]"], []),
    "pp.filter_cells": (["np", "sp", "da"], []),
    "pp.filter_genes": (["np", "sp", "da"], []),
    "pp.highly_variable_genes": (["np", "sp", "da"], ["da[sp[csc]]"]),
    "pp.log1p": (["np", "sp", "da"], []),
    "pp.neighbors": (["np", "sp"], []),
    "pp.normalize_total": (["np", "sp[csr]", "da"], []),
    "pp.pca": (["np", "sp", "da"], ["da[sp[csc]]"]),
    "pp.regress_out": (["np"], []),
    "pp.sample": (["np", "sp", "da"], []),
    "pp.scale": (["np", "sp", "da"], []),
    "pp.scrublet": (["np", "sp"], []),
    "pp.scrublet_simulate_doublets": (["np", "sp"], []),
    "tl.dendrogram": (["np", "sp"], []),
    "tl.diffmap": (["np", "sp"], []),
    "tl.dpt": (["np", "sp"], []),
    "tl.draw_graph": (["np", "sp"], []),  # only uses graph in obsp
    "tl.embedding_density": (["np"], []),
    "tl.ingest": (["np", "sp"], []),
    "tl.leiden": (["np", "sp"], []),  # only uses graph in obsp
    "tl.louvain": (["np", "sp"], []),  # only uses graph in obsp
    "tl.paga": (["np", "sp"], []),
    "tl.rank_genes_groups": (["np", "sp"], []),
    "tl.tsne": (["np", "sp"], []),
    "tl.umap": (["np", "sp"], []),
}


# -- Options for HTML output ----------------------------------------------

# The theme is sphinx-book-theme, with patches for readthedocs-sphinx-search
html_theme = "scanpydoc"
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False
html_logo = "_static/img/Scanpy_Logo_BrightFG.svg"
html_title = "scanpy"


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,  # noqa: FBT003
    )


# -- Options for other output formats ------------------------------------------

htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]


# -- Suppress link warnings ----------------------------------------------------

qualname_overrides = {
    "pathlib._local.Path": "pathlib.Path",
    "sklearn.neighbors._dist_metrics.DistanceMetric": "sklearn.metrics.DistanceMetric",
    "scanpy.plotting._matrixplot.MatrixPlot": "scanpy.pl.MatrixPlot",
    "scanpy.plotting._dotplot.DotPlot": "scanpy.pl.DotPlot",
    "scanpy.plotting._stacked_violin.StackedViolin": "scanpy.pl.StackedViolin",
    "scanpy._param_sets.HVGFlavor": "tuple",
    "scanpy._param_sets.FilterCellsCutoffs": "tuple",
    "scanpy._param_sets.FilterGenesCutoffs": "tuple",
    "pandas.core.series.Series": "pandas.Series",
    # https://github.com/pandas-dev/pandas/issues/63810
    "pandas.api.typing.aliases.AnyArrayLike": ("doc", "pandas:reference/aliases"),
    "numpy.bool_": "numpy.bool",  # Since numpy 2, numpy.bool is the canonical dtype
    "numpy.typing.ArrayLike": ("py:data", "numpy.typing.ArrayLike"),
}

nitpick_ignore = [
    # Technical issues
    ("py:class", "numpy.int64"),  # documented as “attribute”
    ("py:class", "numpy._typing._dtype_like._SupportsDType"),
    ("py:class", "numpy._typing._dtype_like._DTypeDict"),
    # Will probably be documented
    ("py:class", "scanpy._settings.Verbosity"),
    ("py:class", "scanpy.neighbors.OnFlySymMatrix"),
    ("py:class", "scanpy.plotting._baseplot_class.BasePlot"),
    # Currently undocumented
    # https://github.com/mwaskom/seaborn/issues/1810
    ("py:class", "seaborn.matrix.ClusterGrid"),
    ("py:class", "samalg.SAM"),
    # Won’t be documented
    ("py:class", "scanpy.plotting._utils._AxesSubplot"),
    ("py:class", "scanpy._utils.Empty"),
    ("py:class", "numpy.random.mtrand.RandomState"),
    ("py:class", "scanpy.neighbors._types.KnnTransformerLike"),
]

# Options for plot examples

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root

# link config
issues_github_path = "scverse/scanpy"
rtd_links_prefix = PurePosixPath("src")
