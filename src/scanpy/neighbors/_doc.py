from __future__ import annotations

from .._docs import doc_out_multi

doc_neighbors_key = """\
neighbors_key
    Where to look for neighbors connectivities.
    If not specified, this retrieves ``.obsp['connectivities']`` for connectivities
    (default storage place for :func:`~scanpy.pp.neighbors`).
    If specified, this retrieves
    ``.obsp[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
"""

doc_use = r"""use
    Use the indicated representation:
    a :class:`~anndata.acc.LayerAcc` (e.g. `A.X`, `A.layers[...]`) or
    :class:`~anndata.acc.MultiAcc` (e.g. `A.obsm[...]`, `A.varm[...]`).
    :class:`str`\ s are :meth:`~anndata.acc.AdAcc.resolve`\ d, e.g. `'obsm.pca'`.

    If `None`, the representation is chosen automatically:
    For `.n_vars` < :attr:`~scanpy.settings.N_PCS` (default: 50), `.X` is used, otherwise the PCA
    representation (`.obsm['X_pca']`, or `.obsm['pca']` if it was computed under
    :attr:`~scanpy.Preset.ScanpyV2Preview`).
    If it is not present, it’s computed with default parameters or `n_pcs` if present."""

doc_use_rep = f"""\
{doc_use}
use_rep
    `'X'` or a key of :attr:`~anndata.AnnData.obsm`,
    i.e. `use_rep='k'` is `use=A.obsm['k']` and `use_rep='X'` is `use=A.X`."""

doc_n_pcs = """\
n_pcs
    Use this many PCs. If `n_pcs==0` use `.X` if `use is None`.\
"""

_doc_outs = dict(
    connectivities=(
        "The weighted adjacency matrix. "
        "By default written to `A.obsp['connectivities']`, "
        "with `sc.Prefix('p')` to `A.obsp['p_connectivities']`."
    ),
    distances=(
        "The distance matrix. "
        "By default written to `A.obsp['distances']`, "
        "with `sc.Prefix('p')` to `A.obsp['p_distances']`."
    ),
)

doc_out = doc_out_multi(
    _doc_outs,
    example="{'connectivities': A.obsp['conns']}",
    returns=":class:`~scanpy.pp.NeighborsOut`",
    uns="neighbors",
)

doc_out_key_added = doc_out_multi(
    _doc_outs,
    example="{'connectivities': A.obsp['conns']}",
    returns=":class:`~scanpy.pp.NeighborsOut`",
    uns="neighbors",
    key_added=(
        "`key_added='k'` wrote `.obsp['k_connectivities']`, `.obsp['k_distances']`, "
        "and `.uns['k']`, so it is `out=sc.Prefix('k')`, "
        "except that the latter names the `.uns` entry `'k_neighbors'`."
    ),
)
