"""External preprocessing functions."""

from __future__ import annotations

from ..._compat import deprecated
from ._bbknn import bbknn
from ._dca import dca
from ._hashsolo import hashsolo
from ._magic import magic
from ._mnn_correct import mnn_correct
from ._scanorama_integrate import scanorama_integrate

__all__ = [
    "bbknn",
    "dca",
    "hashsolo",
    "magic",
    "mnn_correct",
    "scanorama_integrate",
]


@deprecated("Import from sc.pp instead")
def harmony_integrate(*args, **kwargs):
    from ...preprocessing import harmony_integrate

    return harmony_integrate(*args, **kwargs)


@deprecated("Import from sc.pp instead")
def scrublet(*args, **kwargs):
    from ...preprocessing import scrublet

    return scrublet(*args, **kwargs)


@deprecated("Import from sc.pp instead")
def scrublet_simulate_doublets(*args, **kwargs):
    from ...preprocessing import scrublet_simulate_doublets

    return scrublet_simulate_doublets(*args, **kwargs)
