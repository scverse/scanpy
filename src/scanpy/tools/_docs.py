"""Shared docstrings for tools function parameters."""

from __future__ import annotations

doc_adata = """\
adata
    The annotated data matrix.\
"""

doc_restrict_to = """\
restrict_to
    Restrict the clustering to the categories within the key for sample
    annotation, tuple needs to contain `(obs_key, list_of_categories)`.\
"""

doc_adjacency = """\
adjacency
    Sparse adjacency matrix of the graph, defaults to neighbors connectivities.\
"""

doc_neighbors_key = """\
neighbors_key
    Use neighbors connectivities as adjacency.
    If not specified, {method} looks at .obsp['connectivities'] for connectivities
    (default storage place for pp.neighbors).
    If specified, {method} looks at
    .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.\
"""

doc_obsp = """\
obsp
    Use .obsp[obsp] as adjacency. You can't specify both
    `obsp` and `neighbors_key` at the same time.\
"""


def doc_use_neighbors(*, meta: bool = False, extra: str = "") -> str:
    """Docs for a `use` parameter choosing a neighborhood graph.

    With `meta`, the function also uses what :func:`~scanpy.pp.neighbors` stored alongside it.
    `extra` is appended as a separate paragraph.
    """
    meta_doc = (
        "\n    Distances and parameters are used from where :func:`~scanpy.pp.neighbors` "
        "stored them alongside it, if they exist."
        if meta
        else ""
    )
    return f"""\
use
    The neighborhood graph (connectivities) to use, e.g. `A.obsp['k']`
    after `sc.pp.neighbors(adata, out=A.obsp['k'])`, or any other graph.
    :class:`str`\\ s are :meth:`~anndata.acc.AdAcc.resolve`\\ d, e.g. `'obsp.k'`.
    If :data:`None`, the one :func:`~scanpy.pp.neighbors` wrote by default.{meta_doc}{extra}\
"""
