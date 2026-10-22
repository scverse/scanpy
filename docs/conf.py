"""Configuration for Scanpy’s Sphinx documentation."""

from __future__ import annotations

import os
import shutil
import sys
from datetime import datetime
from functools import partial
from importlib.metadata import version as get_version
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib  # noqa
from docutils import nodes
from packaging.version import Version
from scverse_doc.config import DEFAULTS
from sphinxcontrib.katex import NODEJS_BINARY

if TYPE_CHECKING:
    from sphinx.application import Sphinx

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]


# -- General configuration ------------------------------------------------

needs_sphinx = "4.0"  # Nicer param docs

# General information
project = "Scanpy"
author = "Scanpy development team"
source_repository = "https://github.com/scverse/scanpy"
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
exclude_patterns = [
    *DEFAULTS.exclude_patterns,
    # exclude all 0.x.y.md files, but not index.md
    "release-notes/[!i]*.md",
]

extensions = [
    "scverse_doc",  # sets up napoleon and so needs to be before s-a-t
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx_exec_jupyter",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.katex",
    "sphinx_autodoc_typehints",
    "scverse_misc.sphinx_ext",
    "sphinx.ext.linkcode",
    "sphinx_issues",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

# Generate the API documentation when building
autodoc_typehints = "none"
autodoc_member_order = "bysource"
autodoc_default_options = {
    # Don’t show members in addition to the autosummary table added by `_templates/class.rst`
    "members": False,
    # show “Bases: SomeClass” at the top of class docs
    "show-inheritance": True,
}
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
api_dir = HERE / "api"  # function_images
exec_jupyter_code = """
# setup notebook backend
import matplotlib
matplotlib.use("module://matplotlib_inline.backend_inline")
# import all slow optional imports before running code
import scanpy, umap, seaborn, sklearn.metrics, pynndescent, networkx
del scanpy, umap, seaborn, sklearn, pynndescent, networkx, matplotlib
"""
holoviews_backends = ["bokeh", "matplotlib", "plotly"]
myst_url_schemes = ("http", "https", "mailto", "ftp")
myst_heading_anchors = 3
myst_ignore_mime_types = [  # from custom extension patch_myst_nb
    "application/vnd.microsoft.datawrangler.viewer.v0+json",
]
nb_execution_mode = "cache"
nb_execution_excludepatterns = [
    f"{d}{'/*' * n}" for d in ["tutorials", "how-to"] for n in (1, 2, 3)
]
nb_execution_show_tb = bool(os.environ.get("READTHEDOCS"))

ogp_site_url = "https://scanpy.scverse.org/en/stable/"
ogp_image = f"{ogp_site_url}_static/Scanpy_Logo_BrightFG.svg"

typehints_defaults = "braces"
always_use_bars_union = True  # Don’t use `Union` even when building with Python ≤3.14

katex_prerender = shutil.which(NODEJS_BINARY) is not None

intersphinx_mapping = dict(
    anndata=("https://anndata.scverse.org/en/stable/", None),
    bbknn=("https://bbknn.readthedocs.io/en/latest/", None),
    cellbrowser=("https://cellbrowser.readthedocs.io/en/master/", None),
    cellrank=("https://cellrank.readthedocs.io/en/stable/", None),
    cuml=("https://docs.rapids.ai/api/cuml/stable/", None),
    cycler=("https://matplotlib.org/cycler/", None),
    dask=("https://docs.dask.org/en/stable/", None),
    dask_ml=("https://ml.dask.org/", None),
    decoupler=("https://decoupler.scverse.org/en/stable/", None),
    dynamo=("https://dynamo-release.readthedocs.io/en/stable/", None),
    fast_array_utils=("https://fast-array-utils.scverse.org/en/stable/", None),
    flit=("https://flit.pypa.io/en/stable/", None),
    glasbey=("https://glasbey.readthedocs.io/en/stable/", None),
    h5py=("https://docs.h5py.org/en/stable/", None),
    holoviews=("https://holoviews.org/", None),
    hv_anndata=("https://hv-anndata.readthedocs.io/en/latest/", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    ipython=("https://ipython.readthedocs.io/en/stable/", None),
    igraph=("https://python.igraph.org/en/stable/api/", None),
    leidenalg=("https://leidenalg.readthedocs.io/en/latest/", None),
    louvain=("https://louvain-igraph.readthedocs.io/en/latest/", None),
    marsilea=("https://marsilea.readthedocs.io/en/stable/", None),
    matplotlib=("https://matplotlib.org/stable/", None),
    mudata=("https://mudata.scverse.org/stable/", None),
    napoleon=("https://sphinxcontrib-napoleon.readthedocs.io/en/latest/", None),
    networkx=("https://networkx.org/documentation/stable/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    palantir=("https://palantir.readthedocs.io/en/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    phate=("https://phate.readthedocs.io/en/stable/", None),
    pydeseq2=("https://pydeseq2.readthedocs.io/en/stable/", None),
    pynndescent=("https://pynndescent.readthedocs.io/en/latest/", None),
    pytest=("https://docs.pytest.org/en/latest/", None),
    python=("https://docs.python.org/3", None),
    rapids_singlecell=("https://rapids-singlecell.scverse.org/en/stable/", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    scverse_misc=("https://scverse-misc.readthedocs.io/stable/", None),
    scverse_tutorials=("https://scverse-tutorials.readthedocs.io/en/latest/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    session_info2=("https://session-info2.readthedocs.io/en/stable/", None),
    spatialdata=("https://spatialdata.scverse.org/en/stable/", None),
    squidpy=("https://squidpy.readthedocs.io/en/stable/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    towncrier=("https://towncrier.readthedocs.io/en/stable/", None),
)


array_support: dict[str, tuple[list[str], list[str]]] = {
    "experimental.pp.highly_variable_genes": (["np", "sp"], []),
    "get.aggregate": (["np", "sp", "da"], []),
    "pp.bbknn": (["np", "sp"], []),
    "pp.calculate_qc_metrics": (["np", "sp", "da"], []),
    "pp.combat": (["np"], []),
    "pp.downsample_counts": (["np", "sp[csr]"], []),
    "pp.filter_cells": (["np", "sp", "da"], []),
    "pp.filter_genes": (["np", "sp", "da"], []),
    "pp.harmony_integrate": (["np"], []),
    "pp.hashsolo": (["np", "sp"], []),
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
    "pp.wnn": (["np"], []),
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
    "tl.score_genes": (["np", "sp"], []),
    "tl.tsne": (["np", "sp"], []),
    "tl.umap": (["np", "sp"], []),
}


# -- Options for HTML output ----------------------------------------------

html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False
html_logo = "_static/img/Scanpy_Logo_BrightFG.svg"
html_title = "scanpy"


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))


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

xref_overrides = {
    "numpy.typing.DTypeLike": ("py:data", "numpy.typing.DTypeLike"),
    "numpy._typing._array_like.NDArray": ("py:data", "numpy.typing.NDArray"),
}

nitpick_ignore = [
    # Technical issues
    ("py:class", "numpy.int64"),  # documented as “attribute”
    ("py:class", "numpy._typing._dtype_like._SupportsDType"),
    ("py:obj", "numpy._typing._dtype_like._SupportsDType"),
    ("py:class", "numpy._typing._dtype_like._DTypeDict"),
    # Will probably be documented
    ("py:class", "scanpy.neighbors.OnFlySymMatrix"),
    ("py:class", "scanpy.pl.BasePlot"),
    # Currently undocumented
    # https://github.com/mwaskom/seaborn/issues/1810
    ("py:class", "seaborn.matrix.ClusterGrid"),
    ("py:class", "samalg.SAM"),
    # Won’t be documented
    ("py:class", "scanpy.pl._AxesSubplot"),
    ("py:class", "scanpy._utils.Empty"),
    ("py:class", "numpy.random.mtrand.RandomState"),
    ("py:class", "scanpy.neighbors._types.KnnTransformerLike"),
    ("py:class", "scanpy._settings.presets.Default"),
]

# Options for plot examples

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root

# link config
issues_github_path = "scverse/scanpy"
