from __future__ import annotations
from importlib.util import find_spec
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pytest


# Mapping from module name to PyPI name
KNOWN = dict(
    leidenalg='leidenalg',
    louvain='louvain',
    skmisc='scikit-misc',
    fa2='fa2',
    igraph='igraph',
    dask='dask',
    zarr='zarr',
    zappy='zappy',
    pybiomart='pybiomart',
    gprofiler='gprofiler-official',
    # external
    bbknn='bbknn',
    harmony='harmonyTS',
    harmonypy='harmonypy',
    magic='magic-impute',
    palantir='palantir',
    phate='phate',
    phenograph='PhenoGraph',
    pypairs='pypairs',
    samalg='sam-algorithm',
    scanorama='scanorama',
    scrublet='scrublet',
    trimap='trimap',
    wishbone='wishbone-dev',
)


def needs(mod: str) -> pytest.MarkDecorator:
    """
    Returns a pytest skip marker evaluated at module import.

    This allows us to see the amount of skipped tests at the start of a test run.
    :func:`pytest.importorskip` skips tests after they started running.
    """

    try:
        import pytest
    except ImportError:
        return lambda func: func

    try:
        dist = KNOWN[mod]
    except KeyError:
        raise ValueError(
            f"Unknown import {mod}. Please add to KNOWN in this file."
        ) from None
    reason = f"needs module `{mod}`"
    if mod != dist.lower().replace("-", "_"):
        reason = f"{reason} (`pip install {dist}`)"
    return pytest.mark.skipif(not find_spec(mod), reason=reason)
