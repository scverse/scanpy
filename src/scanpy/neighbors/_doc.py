from __future__ import annotations

doc_use_rep = """\
use_rep
    Use the indicated representation. `'X'` or any key for `.obsm` is valid.
    If `None`, the representation is chosen automatically:
    For `.n_vars` < :attr:`~scanpy.settings.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
    If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.\
"""

doc_n_pcs = """\
n_pcs
    Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\
"""
