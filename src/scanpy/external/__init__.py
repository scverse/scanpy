"""External preprocessing and analysis tools and their plotting."""

from __future__ import annotations

import sys
import warnings

from .. import _utils
from . import exporting, pl, pp, tl

__all__: list[str] = ["exporting", "pl", "pp", "tl"]

_utils.annotate_doc_types(sys.modules[__name__], "scanpy")

msg = "The `scanpy.external` module is deprecated and will be removed in a future version."
warnings.warn(msg, DeprecationWarning, stacklevel=2)

del msg, sys, _utils, warnings
