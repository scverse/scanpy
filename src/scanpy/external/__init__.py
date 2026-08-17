"""External preprocessing and analysis tools and their plotting."""

from __future__ import annotations

import sys

from .. import _utils
from . import _exporting as exporting
from . import pl, pp, tl

sys.modules[f"{__name__}.exporting"] = exporting

_utils.annotate_doc_types(sys.modules[__name__], "scanpy")
del sys, _utils

__all__ = ["pl", "pp", "tl"]
