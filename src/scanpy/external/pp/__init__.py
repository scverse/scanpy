"""External preprocessing functions."""

from __future__ import annotations

from sklearn.utils import deprecated

from ...preprocessing import _scrublet
from ._bbknn import bbknn
from ._dca import dca
from ._harmony_integrate import harmony_integrate
from ._hashsolo import hashsolo
from ._magic import magic
from ._mnn_correct import mnn_correct
from ._scanorama_integrate import scanorama_integrate

scrublet = deprecated("Import from sc.pp instead")(_scrublet.scrublet)
scrublet_simulate_doublets = deprecated("Import from sc.pp instead")(
    _scrublet.scrublet_simulate_doublets
)

__all__ = [
    "bbknn",
    "dca",
    "harmony_integrate",
    "hashsolo",
    "magic",
    "mnn_correct",
    "scanorama_integrate",
]
