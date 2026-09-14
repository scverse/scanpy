from __future__ import annotations

doc_neighbors_key = """\
neighbors_key
    Where to look for neighbors connectivities.
    If not specified, this retrieves ``.obsp['connectivities']`` for connectivities
    (default storage place for :func:`~scanpy.pp.neighbors`).
    If specified, this retrieves
    ``.obsp[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
"""

doc_use_rep = r"""use_rep
    Use the indicated representation:
    a :class:`~anndata.acc.LayerAcc` (e.g. `A.X`, `A.layers[...]`) or
    :class:`~anndata.acc.MultiAcc` (e.g. `A.obsm[...]`, `A.varm[...]`).
    A :class:`str` is :meth:`~anndata.acc.AdAcc.resolve`\ d to one of those
    if :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`,
    otherwise interpreted as `'X'` or a key of `.obsm`.

    If `None`, the representation is chosen automatically:
    For `.n_vars` < :attr:`~scanpy.settings.N_PCS` (default: 50), `.X` is used, otherwise the PCA
    representation (`.obsm['X_pca']`, or `.obsm['pca']` if it was computed under
    :attr:`~scanpy.Preset.ScanpyV2Preview`).
    If it is not present, it’s computed with default parameters or `n_pcs` if present."""

doc_n_pcs = """\
n_pcs
    Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\
"""
