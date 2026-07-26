from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import holoviews as hv
from hv_anndata import A, register

from ..._utils import get_literal_vals

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import TypeGuard

    from anndata.acc import AdAcc


__all__ = ["hv_init"]

_Backend = Literal["bokeh", "matplotlib", "plotly"]


def hv_init(*backends: _Backend | None) -> AdAcc:
    """Shortcut to initialize :mod:`hv_anndata`.

    #. Calls :func:`hv_anndata.register`
    #. Initializes a :func:`holoviews.extension` with the first available backend.
    #. returns :data:`hv_anndata.A`.

    Parameters
    ----------
    backends
        Unless only `None` is passed, these backends are attempted to be initialized in this order.

    Examples
    --------
    >>> import scanpy as sc
    >>> A = sc.pl.hv_init()
    """
    register()
    if backends != (None,):
        if not _all_match(backends):
            msg = "`*backends` needs to be `None` or 0 or more of 'bokeh' | 'matplotlib' | 'plotly'."
            raise ValueError(msg)
        hv.extension(*(backends or hv.Store.renderers.keys() or ("bokeh",)))
    return A


def _all_match(backends: Iterable[object]) -> TypeGuard[Iterable[_Backend]]:
    return set(backends) <= get_literal_vals(_Backend)
