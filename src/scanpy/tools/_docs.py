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


def doc_use_neighbors(*, graph: bool = False) -> str:
    """Docs for a `use` parameter finding what :func:`~scanpy.pp.neighbors` wrote.

    With `graph`, the function only needs the connectivities, so `use` can also be an `.obsp` accessor.
    """
    graph_doc = (
        "\n    Or any graph, given by an accessor like `A.obsp['my_graph']`."
        if graph
        else ""
    )
    return f"""\
use
    Where :func:`~scanpy.pp.neighbors` wrote its results, i.e. the `out` passed to it:
    `None` for the default places, `sc.Prefix('p')` for the `p_`-prefixed ones,
    or a mapping like `{{'connectivities': A.obsp['conns'], 'distances': A.obsp['dists']}}`.{graph_doc}\
"""
