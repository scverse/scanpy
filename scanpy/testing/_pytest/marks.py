from __future__ import annotations
from importlib.util import find_spec

import pytest


def make_skip_mark(mod: str, dist: str | None = None):
    """
    Returns a pytest skip marker evaluated at module import.

    This allows us to see the amount of skipped tests at the start of a test run.
    :func:`pytest.importorskip` skips tests after they started running.
    """
    reason = f"needs module `{mod}`"
    if dist:
        reason = f"{reason} (`pip install {dist}`)"
    return pytest.mark.skipif(not find_spec(mod), reason=reason)


needs_leidenalg = make_skip_mark("leidenalg")
needs_louvain = make_skip_mark("louvain")
needs_skmisc = make_skip_mark("skmisc", "scikit-misc")
needs_scanorama = make_skip_mark("scanorama")
needs_scrublet = make_skip_mark("scrublet")
needs_fa2 = make_skip_mark("fa2")
needs_igraph = make_skip_mark("igraph", "python-igraph")
