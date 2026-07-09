from __future__ import annotations

from typing import TYPE_CHECKING

from hv_anndata import A, register

if TYPE_CHECKING:
    from anndata.acc import AdAcc


__all__ = ["hv_init"]


def hv_init() -> AdAcc:
    """Shortcut to initialize :mod:`hv_anndata`.

    Simply calls :func:`hv_anndata.register` and returns :class:`hv_anndata.A`.

    Examples
    --------
    >>> import scanpy as sc
    >>> A = sc.pl.hv_init()
    """
    register()
    return A
