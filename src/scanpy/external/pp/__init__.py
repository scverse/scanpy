"""External preprocessing functions."""

from __future__ import annotations

from scverse_misc import Deprecation, deprecated

from ._bbknn import bbknn
from ._harmony_integrate import harmony_integrate
from ._hashsolo import hashsolo
from ._magic import magic
from ._mnn_correct import mnn_correct
from ._scanorama_integrate import scanorama_integrate

__all__ = [
    "bbknn",
    "harmony_integrate",
    "hashsolo",
    "magic",
    "mnn_correct",
    "scanorama_integrate",
]


@deprecated(Deprecation("1.10.0", "Import from sc.pp instead."))
def scrublet(*args, **kwargs):
    from ...preprocessing import scrublet

    return scrublet(*args, **kwargs)


@deprecated(Deprecation("1.10.0", "Import from sc.pp instead."))
def scrublet_simulate_doublets(*args, **kwargs):
    from ...preprocessing import scrublet_simulate_doublets

    return scrublet_simulate_doublets(*args, **kwargs)
